"""
sbml2kin_batch.py -- SBML (L2/L3) -> Kintecus deck + .kintpack (v3), built around probed Kintecus V2027 behavior.

Encoding rules (each verified on the keyless V2027 build under Wine):
  * mass action             k<TAB>A+B==>C                          rate = k*[A]*[B]
  * any other rate law      1<TAB>A+S[USER;cA,cK:(expr)]==>C        rate = expr*[A]  (column 1 ignored)
  * zero order              1 source+S[USER;...:(expr)]==>C         ('source' = constant species 1.0)
  * nothing produced        A+S[USER;...]==>sink                    ('sink' = constant species 1.0; NULL is not used)
  * negative rate terms are sent through a reverse line (Kintecus stalls on negative rates)
  * consumption only via the left side (negative product coefficients are not parsed)
  * USER expression <= 780 chars and <= 50 variables per line; long laws split into additive terms
  * species never changed by any line are pinned: SPECIES.DAT column 7 = value (unpinned modifiers stall)
  * SBML time      -> pseudo-species TIME_ fed by 1 source==>TIME_
  * rate rules     -> USER production/consumption lines on the target (pseudo-species for parameters)
  * time events    -> Constant File profiles: pinned piecewise profile for parameters/constant species,
                      '-SET' / '-ADD' profile for reacting species when the event times allow it
  * function definitions are inlined into the USER expressions
Anything else (state-triggered events, algebraic rules, piecewise, delay, ...) is flagged and the rest converted.
"""
import json, math, re, sys, os
from fractions import Fraction
import libsbml
import sympy as sp

# Kintecus build 2026-11-03 12:39: the reaction-equation field holds 897 characters (898 -> FATAL #24); the MODEL.DAT
# row itself may be longer (a 2,100-character comment runs). Builds up to 2026-10-30 22:51 limited the whole row to 897
# (reaction field 890, USER expression 840). KIN_FIELD_MAX overrides for newer builds.
REACTION_FIELD_MAX = int(os.environ.get("KIN_FIELD_MAX", 2153))
USER_MAX_LEN = REACTION_FIELD_MAX
MODEL_LINE_MAX = REACTION_FIELD_MAX + len('1\t0\t0\t') + 1
LINE_MAX = REACTION_FIELD_MAX   # reaction field of a USER line; comments go on a '#' line above
USER_SAFE = 840         # was 236 for builds before 2026-10-30, whose parsers mis-evaluated (FATAL #999 or silently wrong) USER expressions past
                        # ~248 characters when they contain '^' (fewer with many '^'); len + count('^') is kept <= this
MAX_REACTANTS = 3       # Kintecus V2027: a 4th reactant entry on one line gives FATAL #39
USER_MAX_VARS = 92     # Kintecus build 2026-10-30 22:51: 92 variables per USER line run; 93 -> "missing a variable list" error
NAME_MAX = 26          # Kintecus >= 2026-10-26 15:30: USER variable 'c'+name up to 27 characters (was 10 -> 9)
RESERVED = {"M", "S", "T", "P", "R", "NULL", "HV", "END", "E", "SOURCE", "SINK"}
# J. Ianni 2026-10-30: no NULL in reactions. Nothing consumed -> reactant "source", nothing produced -> product "sink";
# both are SPECIES.DAT entries pinned at 1.0 (column 7), so rate = k*[source]*... = k*... and sink never changes.
# SBML species called Source/Sink are renamed (Source1, Sink1) because these two names are taken.
SOURCE, SINK = "source", "sink"
# Kintecus builds before 2026-10-26 14:44 recognized function/constant names INSIDE USER variable names
# (pi, ln, exp, log, sin, cos, tan, atn, sqrt, abs, floor) - fixed in that build, so no renaming is needed:
BAD_SUB = ()
BAD_SUB_OLD_BUILDS = ("PI", "LN", "EXP", "LOG", "SIN", "COS", "TAN", "ATN", "SQRT", "ABS", "FLOOR")
NAME_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")


class Unsupported(Exception):
    pass


# ----------------------------------------------------------------------------- math -> sympy
STEP_EPS = sp.Float("1e-300")
DIV_EPS = sp.Float(os.environ.get("DIV_EPS", "1e-30"))
DIV_MODE = os.environ.get("DIV_MODE", "split_signed")   # split_signed (2026-11-01) | eps: [A]+1e-30 (old) | abs
DIV_REL = float(os.environ.get("DIV_REL", "1e-10"))
SIGNED_NONPROP = os.environ.get("SIGNED_NONPROP", "0") == "1"   # non-proportional terms: signed source lines
ARULE_SWITCH = os.environ.get("ARULE_SWITCH", "1") == "1"   # time conditions in assignment rules -> switch species
PROF_MAX = 20001       # Kintecus V2027 (2026-10-30 22:51) reads at most 20,002 values from a profile file (FATAL #94)
PROFILE_FORMAT = os.environ.get("KIN_PROFILE_FORMAT", "pairs")   # "pairs": Time/value pairs (needs build 2026-11-05+);
                                                                  # "spacing": negative spacing + one value per step
PAIR_MAX = 1000        # pairs per Constant File (dimension 1052; largest file in the collection: 7 pairs)
SPLIT_LONG = os.environ.get("KIN_SPLIT", "1") == "1"   # rate too long for one line -> several lines
UNROLL = os.environ.get("KIN_UNROLL", "1") == "1"   # time-determined events -> time events
PAIR_FILES_MAX = 95    # pair files per model (J. Ianni 11-05: TSX1(95,1052), TSY1(95,1052))
FRAC_MAX = 4.0          # non-integer reactant coefficients up to this are written as Kintecus coefficients (rate ~ [A]^n)


def step(z):
    """Heaviside step built from abs(): 1 for z>0, 0 for z<0 (Kintecus USER has no conditionals)."""
    if z.is_Number:
        return sp.Integer(1) if z > 0 else sp.Integer(0)
    return (z / (sp.Abs(z) + STEP_EPS) + 1) / 2


class MathConv:
    used_step = False

    def __init__(self, resolve, switch=None):
        self.resolve = resolve          # name -> sympy expr (raises Unsupported)
        self.switch = switch            # (t0, up) -> sympy symbol of a 0/1 pinned profile species

    def __call__(self, ast):
        if ast is None:
            raise Unsupported("missing math")
        return self.conv(ast)

    def conv(self, n):
        t = n.getType()
        T = libsbml
        lazy = (T.AST_FUNCTION_PIECEWISE, T.AST_RELATIONAL_GT, T.AST_RELATIONAL_GEQ, T.AST_RELATIONAL_LT,
                T.AST_RELATIONAL_LEQ, T.AST_RELATIONAL_EQ, T.AST_RELATIONAL_NEQ, T.AST_LOGICAL_AND,
                T.AST_LOGICAL_OR, T.AST_LOGICAL_NOT, T.AST_LOGICAL_XOR)
        k = None if t in lazy else [self.conv(n.getChild(i)) for i in range(n.getNumChildren())]
        T = libsbml
        if t == T.AST_PLUS:
            return sp.Add(*k) if k else sp.Integer(0)
        if t == T.AST_MINUS:
            return -k[0] if len(k) == 1 else k[0] - k[1]
        if t == T.AST_TIMES:
            return sp.Mul(*k) if k else sp.Integer(1)
        if t == T.AST_DIVIDE:
            return k[0] / k[1]
        if t in (T.AST_POWER, T.AST_FUNCTION_POWER):
            return k[0] ** k[1]
        if t == T.AST_INTEGER:
            return sp.Integer(n.getInteger())
        if t in (T.AST_REAL, T.AST_REAL_E):
            v = n.getReal()
            if math.isinf(v) or math.isnan(v):
                raise Unsupported("inf/nan literal")
            return sp.Float(v, 17)
        if t == T.AST_RATIONAL:
            return sp.Rational(n.getNumerator(), n.getDenominator())
        if t == T.AST_NAME:
            return self.resolve(n.getName())
        if t == T.AST_NAME_TIME:
            return self.resolve("__time__")
        if t == T.AST_NAME_AVOGADRO:
            return sp.Float(6.02214076e23, 17)
        if t == T.AST_CONSTANT_E:
            return sp.E
        if t == T.AST_CONSTANT_PI:
            return sp.pi
        if t == T.AST_FUNCTION_EXP:
            return sp.exp(k[0])
        if t == T.AST_FUNCTION_LN:
            return sp.log(k[0])
        if t == T.AST_FUNCTION_LOG:
            return sp.log(k[1], k[0]) if len(k) == 2 else sp.log(k[0], 10)
        if t == T.AST_FUNCTION_ROOT:
            return k[1] ** (1 / k[0]) if len(k) == 2 else sp.sqrt(k[0])
        if t == T.AST_FUNCTION_ABS:
            return sp.Abs(k[0])
        fn = {T.AST_FUNCTION_SIN: sp.sin, T.AST_FUNCTION_COS: sp.cos, T.AST_FUNCTION_TAN: sp.tan,
              T.AST_FUNCTION_ARCSIN: sp.asin, T.AST_FUNCTION_ARCCOS: sp.acos, T.AST_FUNCTION_ARCTAN: sp.atan,
              T.AST_FUNCTION_SINH: sp.sinh, T.AST_FUNCTION_COSH: sp.cosh, T.AST_FUNCTION_TANH: sp.tanh,
              T.AST_FUNCTION_FLOOR: sp.floor}
        if t in fn:
            return fn[t](k[0])
        if t == T.AST_FUNCTION_CEILING:
            return sp.ceiling(k[0])
        if t == T.AST_LAMBDA:
            raise Unsupported("lambda")
        if t == T.AST_FUNCTION_PIECEWISE:
            kids = [n.getChild(i) for i in range(n.getNumChildren())]
            pairs = [(kids[i], kids[i + 1]) for i in range(0, len(kids) - 1, 2)]
            other = self.conv(kids[-1]) if len(kids) % 2 == 1 else sp.Integer(0)
            out = other
            for val, cond in reversed(pairs):
                c = self.conv(cond)
                out = self.conv(val) * c + (1 - c) * out
            self.used_step = True
            return out
        rel = {T.AST_RELATIONAL_GT: 1, T.AST_RELATIONAL_GEQ: 1, T.AST_RELATIONAL_LT: -1, T.AST_RELATIONAL_LEQ: -1}
        if t in rel:
            kids = [self.conv(n.getChild(i)) for i in range(n.getNumChildren())]
            tsym = sp.Symbol("cTIME_")
            if len(kids) == 2 and self.switch is not None and (
                    (kids[0] == tsym and kids[1].is_Number) or (kids[1] == tsym and kids[0].is_Number)):
                t0 = float(kids[1] if kids[0] == tsym else kids[0])
                up = (rel[t] > 0) == (kids[0] == tsym)        # time > t0  -> 0 before, 1 after
                return self.switch(t0, up)
            out = sp.Integer(1)
            for a, b in zip(kids, kids[1:]):
                out *= step(rel[t] * (a - b))
            self.used_step = True
            return out
        if t in (T.AST_RELATIONAL_EQ, T.AST_RELATIONAL_NEQ):
            a, b = self.conv(n.getChild(0)), self.conv(n.getChild(1))
            ne = step(sp.Abs(a - b))
            self.used_step = True
            return ne if t == T.AST_RELATIONAL_NEQ else 1 - ne
        if t == T.AST_LOGICAL_AND:
            return sp.Mul(*[self.conv(n.getChild(i)) for i in range(n.getNumChildren())])
        if t == T.AST_LOGICAL_OR:
            out = sp.Integer(0)
            for i in range(n.getNumChildren()):
                c = self.conv(n.getChild(i)); out = out + c - out * c
            return out
        if t == T.AST_LOGICAL_NOT:
            return 1 - self.conv(n.getChild(0))
        if t == T.AST_LOGICAL_XOR:
            a, b = self.conv(n.getChild(0)), self.conv(n.getChild(1))
            return a + b - 2 * a * b
        if t == T.AST_CONSTANT_TRUE:
            return sp.Integer(1)
        if t == T.AST_CONSTANT_FALSE:
            return sp.Integer(0)
        if t == T.AST_FUNCTION_DELAY:
            raise Unsupported("delay")
        if t == T.AST_FUNCTION:
            raise Unsupported(f"function '{n.getName()}' not inlined")
        raise Unsupported(f"MathML node '{n.getName() or t}'")


# ----------------------------------------------------------------------------- printing for Kintecus USER
def num(v, digits=17):
    v = float(v)
    if not math.isfinite(v):
        raise Unsupported("non-finite number")
    if v == int(v) and abs(v) < 1e15:
        return str(int(v))
    # default: exact round-trip (shortest repr) - rounding can put a value exactly on a singularity
    # (e.g. 2.99999999999998 -> 3); shorter forms are used only for long USER expressions (see emit)
    out = repr(v) if digits >= 17 else f"{v:.{digits}g}"
    if "e" in out:
        m, e = out.split("e")
        out = f"{m}e{int(e)}"                # 1e-30, 1.5e16 (no '+', no leading zeros)
    return out


from sympy.printing.str import StrPrinter
from sympy.printing.precedence import precedence


OFFSET = sp.Function("OFFSET")     # OFFSET(cX_o, shift) = cX_o - shift; opaque so expansion never splits it


class KPrinter(StrPrinter):
    """Minimal-parenthesis printer in Kintecus USER syntax: ^ powers, (-1*x) negative exponents, log = ln, atn."""
    digits = 17
    nopow = False      # print powers without '^' (Kintecus V2027 mis-evaluates long USER expressions with '^')

    def _print_Float(self, e):
        v = float(e)
        if math.isinf(v) and e.is_finite:
            # finite but beyond double range (sympy keeps e.g. exp(3800) exactly): print the largest double, which
            # gives what Kintecus (IEEE double) computes - m/(1+exp(3800)) -> 0
            v = math.copysign(1e308, v)
        return num(v, self.digits)

    def _print_Integer(self, e):
        return str(int(e))

    def _print_Rational(self, e):
        return num(float(e), self.digits)

    def _print_Pow(self, e, rational=False):
        b, x = e.base, e.exp
        if x.is_Number and float(x) == 1.0:      # e.g. Pow(cB^12.5033, 1.0) left unsimplified by sympy
            # keep the base's parentheses: the caller placed this node as a power (tight binding), so a bare
            # sum here would lose them - (a)*(b-c)^1.0 must not print as a*b-c
            return self.parenthesize(b, precedence(e), strict=True)
        if x == sp.Rational(1, 2):
            return f"sqrt({self._print(b)})"
        if x.is_Number and x < 0:
            if x == -1:
                return f"1/{self.parenthesize(b, precedence(e))}"
            return f"1/{self._print(b ** (-x)) if (b ** (-x)).is_Symbol else '(' + self._print(b ** (-x)) + ')'}"
        if self.nopow:
            atom = lambda z: z.is_Symbol or (z.is_Number and z >= 0)
            if x.is_Integer and 1 < x <= 6:
                bs = self._print(b) if atom(b) else f"({self._print(b)})"
                return "(" + "*".join([bs] * int(x)) + ")"
            xs = self._print(x)
            xs = xs if atom(x) else f"({xs})"
            if b.is_Number and b > 0:
                return f"exp({num(math.log(float(b)), self.digits)}*{xs})"
            return f"exp({xs}*log({self._print(b)}))"
        bs = self.parenthesize(b, precedence(e), strict=True)
        if b.is_Pow and not bs.startswith("("):  # Kintecus cannot read a^b^c: write (a^b)^c
            bs = f"({bs})"
        xs = self._print(x)
        if not (x.is_Symbol or (x.is_Number and x > 0)):
            xs = f"({xs})"
        return f"{bs}^{xs}"

    def _print_Mul(self, e):
        c, rest = e.as_coeff_Mul()
        if c < 0:
            return f"(-1*{self._print(-e)})"
        return super()._print_Mul(e)

    def _print_Add(self, e, order=None):
        parts = []
        for i, a in enumerate(e.as_ordered_terms()):
            c, _ = a.as_coeff_Mul()
            t = self._print(-a) if c < 0 else self._print(a)
            parts.append(("-" if c < 0 else ("+" if i else "")) + t)
        return "".join(parts)

    def _print_Function(self, e):
        if e.func.__name__ == "OFFSET":          # shifted state: (c<state> - shift), kept as one factor
            return f"({self._print(e.args[0])}-{self._print(e.args[1])})"
        fmap = {"exp": "exp", "log": "log", "sin": "sin", "cos": "cos", "tan": "tan", "asin": "asin",
                "acos": "acos", "atan": "atn", "sinh": "sinh", "cosh": "cosh", "tanh": "tanh", "Abs": "abs",
                "floor": "floor"}
        nm = type(e).__name__
        if nm == "log" and len(e.args) == 2:
            return f"log({self._print(e.args[0])})/log({self._print(e.args[1])})"
        if nm == "ceiling":
            return f"(-1*floor(-1*{self._print(e.args[0])}))"
        if nm not in fmap:
            raise Unsupported(f"function {nm} not available in Kintecus USER")
        return f"{fmap[nm]}({self._print(e.args[0])})"

    _print_exp = _print_log = _print_Abs = _print_floor = _print_ceiling = _print_Function

    def _print_Exp1(self, e):
        return "exp(1)"

    def _print_Pi(self, e):
        return repr(math.pi)

    def emptyPrinter(self, e):
        raise Unsupported(f"cannot print {type(e).__name__}")


def ulen(body):
    return len(body) + body.count("^")


def kprint(e, digits=17, nopow=False):
    pr = KPrinter(); pr.digits = digits; pr.nopow = nopow
    s = pr.doprint(e)
    if "**" in s or "Piecewise" in s or re.search(r"(?<![A-Za-z0-9_])(zoo|nan|NaN|inf|oo)(?![A-Za-z0-9_])", s):
        raise Unsupported("expression contains constructs Kintecus USER cannot evaluate")
    return s


# ----------------------------------------------------------------------------- the converter
# ---------------------------------------------------------------------------- time-determined events
class _TD(Exception):
    pass


def _feq(a, b):
    return abs(a - b) <= 1e-9 * max(1.0, abs(a), abs(b))


def safe_cancel(u):
    """sp.cancel (removes common factors, e.g. a reactant concentration a rate is divided by) - but cancel writes
    the result as expanded polynomial / expanded polynomial, which turns factored rate laws such as
    (1+a/Ka)(1+b/Kb)... into thousands of characters. When the cancelled form would not fit on a Kintecus line,
    the expression is kept as written."""
    c = sp.cancel(u)
    if len(str(c)) > USER_MAX_LEN and len(str(c)) > 2 * len(str(u)):
        return u                 # keep the SBML's own (factored) form; sp.factor can run for hours on these
    return c


def split_rate(e):
    """terms whose sum is e: the terms of a sum, or of the largest sum inside a product (the other factors
    multiplied into every term); [e] when e is not a sum."""
    if e.is_Add:
        out = []
        for t in e.args:
            sub = split_rate(t) if not t.is_Add else list(t.args)
            out.extend(sub)
        return out
    if e.is_Mul:
        adds = [f for f in e.args if f.is_Add]
        if adds:
            big = max(adds, key=lambda f: len(str(f)))
            rest = sp.Mul(*[f for f in e.args if f is not big], evaluate=False)
            out = []
            for t in big.args:
                out.extend(split_rate(sp.Mul(rest, t, evaluate=False)) if t.is_Mul or t.is_Add else [sp.Mul(rest, t, evaluate=False)])
            return out
    return [e]


class EventUnroller:
    """Exact discrete simulation of events whose triggers depend only on time, constants and variables changed
    only by these events (time windows, periodic doses written with % or floor, self-resetting clocks, time == c).
    Between firings every trigger is piecewise-affine in time, so its changes are found exactly: breakpoints of
    floor / ceil / % with affine arguments, then the root of each relation inside each affine segment."""
    T = libsbml

    def __init__(self, m, init, t_end, dynamic):
        self.m, self.t_end = m, float(t_end)
        self.events = list(m.getListOfEvents())
        self.consts = {}
        for coll in (m.getListOfParameters(), m.getListOfCompartments(), m.getListOfSpecies()):
            for x in coll:
                v = init.get(x.getId())
                if v is None:
                    v = x.getValue() if hasattr(x, "getValue") and x.isSetValue() else \
                        x.getSize() if hasattr(x, "getSize") and x.isSetSize() else \
                        x.getInitialConcentration() if hasattr(x, "isSetInitialConcentration") and x.isSetInitialConcentration() else \
                        x.getInitialAmount() if hasattr(x, "isSetInitialAmount") and x.isSetInitialAmount() else None
                if v is not None:
                    self.consts[x.getId()] = float(v)
        targets = {a.getVariable() for e in self.events for a in e.getListOfEventAssignments()}
        self.P = {v for v in targets if v not in dynamic}          # discrete variables (changed only by events)
        self.dynamic = dynamic
        for e in self.events:
            if e.isSetDelay():
                raise _TD("event with delay")
            if not e.isSetTrigger() or e.getTrigger().getMath() is None:
                raise _TD("event without trigger")
            for s in self.symbols(e.getTrigger().getMath()):
                if s != "<time>" and s not in self.P and (s in dynamic or s not in self.consts):
                    raise _TD(f"trigger uses '{s}', which changes continuously")

    def symbols(self, ast, out=None):
        out = set() if out is None else out
        if ast is None:
            return out
        if ast.getType() == self.T.AST_NAME:
            out.add(ast.getName())
        elif ast.getType() == self.T.AST_NAME_TIME:
            out.add("<time>")
        for i in range(ast.getNumChildren()):
            self.symbols(ast.getChild(i), out)
        return out

    # evaluation ------------------------------------------------------------------------
    def ev(self, a, t, P):
        T, ty = self.T, a.getType()
        ch = [a.getChild(i) for i in range(a.getNumChildren())]
        E = lambda x: self.ev(x, t, P)
        if ty == T.AST_NAME_TIME:
            return t
        if ty == T.AST_NAME:
            n = a.getName()
            if n in P:
                return P[n]
            if n in self.consts:
                return self.consts[n]
            raise _TD(f"unknown symbol '{n}'")
        if ty in (T.AST_INTEGER, T.AST_REAL, T.AST_REAL_E, T.AST_RATIONAL):
            return float(a.getValue())
        if ty == T.AST_CONSTANT_TRUE: return True
        if ty == T.AST_CONSTANT_FALSE: return False
        if ty == T.AST_CONSTANT_PI: return math.pi
        if ty == T.AST_CONSTANT_E: return math.e
        if ty == T.AST_PLUS: return sum(E(c) for c in ch)
        if ty == T.AST_MINUS: return -E(ch[0]) if len(ch) == 1 else E(ch[0]) - E(ch[1])
        if ty == T.AST_TIMES:
            v = 1.0
            for c in ch: v *= E(c)
            return v
        if ty == T.AST_DIVIDE: return E(ch[0]) / E(ch[1])
        if ty in (T.AST_POWER, T.AST_FUNCTION_POWER): return E(ch[0]) ** E(ch[1])
        if ty == T.AST_FUNCTION_FLOOR: return math.floor(E(ch[0]) + 1e-12 * max(1.0, abs(E(ch[0]))))
        if ty == T.AST_FUNCTION_CEILING: return math.ceil(E(ch[0]) - 1e-12 * max(1.0, abs(E(ch[0]))))
        if ty == T.AST_FUNCTION_ABS: return abs(E(ch[0]))
        if ty == T.AST_FUNCTION_EXP: return math.exp(E(ch[0]))
        if ty == T.AST_FUNCTION_LN: return math.log(E(ch[0]))
        if ty == getattr(T, "AST_FUNCTION_REM", -1):
            x, y = E(ch[0]), E(ch[1]); r = math.fmod(x, y)
            return 0.0 if _feq(r, 0) or _feq(abs(r), abs(y)) else r
        if ty == getattr(T, "AST_FUNCTION_QUOTIENT", -1):
            return float(math.trunc(E(ch[0]) / E(ch[1]) + 1e-12))
        if ty == T.AST_LOGICAL_AND: return all(E(c) for c in ch)
        if ty == T.AST_LOGICAL_OR: return any(E(c) for c in ch)
        if ty == T.AST_LOGICAL_NOT: return not E(ch[0])
        if ty == T.AST_LOGICAL_XOR: return sum(bool(E(c)) for c in ch) % 2 == 1
        if ty in (T.AST_RELATIONAL_EQ, T.AST_RELATIONAL_NEQ, T.AST_RELATIONAL_LT, T.AST_RELATIONAL_LEQ,
                  T.AST_RELATIONAL_GT, T.AST_RELATIONAL_GEQ):
            x, y = E(ch[0]), E(ch[1])
            eq = _feq(float(x), float(y))
            return {T.AST_RELATIONAL_EQ: eq, T.AST_RELATIONAL_NEQ: not eq, T.AST_RELATIONAL_LT: x < y and not eq,
                    T.AST_RELATIONAL_LEQ: x < y or eq, T.AST_RELATIONAL_GT: x > y and not eq,
                    T.AST_RELATIONAL_GEQ: x > y or eq}[ty]
        if ty == T.AST_FUNCTION_PIECEWISE:
            for i in range(0, len(ch) - 1, 2):
                if E(ch[i + 1]):
                    return E(ch[i])
            return E(ch[-1]) if len(ch) % 2 else float("nan")
        raise _TD(f"unsupported function in trigger ({libsbml.formulaToL3String(a)[:40]})")

    # critical times -------------------------------------------------------------------
    def affine(self, a, P):
        """(slope, intercept) if a is affine in time over [0, t_end], else None."""
        ts = [0.0, 0.37 * self.t_end, self.t_end]
        try:
            v = [float(self.ev(a, t, P)) for t in ts]
        except (ZeroDivisionError, ValueError, OverflowError):
            return None
        s = (v[2] - v[0]) / self.t_end
        if not _feq(v[1], v[0] + s * ts[1]):
            return None
        return s, v[0]

    def critical(self, P, lo):
        """sorted times in (lo, t_end] at which some trigger may change value."""
        T = self.T
        steps = (T.AST_FUNCTION_FLOOR, T.AST_FUNCTION_CEILING, getattr(T, "AST_FUNCTION_REM", -1),
                 getattr(T, "AST_FUNCTION_QUOTIENT", -1))
        rel = (T.AST_RELATIONAL_EQ, T.AST_RELATIONAL_NEQ, T.AST_RELATIONAL_LT, T.AST_RELATIONAL_LEQ,
               T.AST_RELATIONAL_GT, T.AST_RELATIONAL_GEQ)
        brk, rels = set(), []

        def walk(a):
            ty = a.getType()
            if ty in steps:
                arg = a.getChild(0) if ty in (T.AST_FUNCTION_FLOOR, T.AST_FUNCTION_CEILING) else None
                if arg is None:                                 # x % y, quotient(x, y): breaks where x/y is integer
                    num, den = self.affine(a.getChild(0), P), self.affine(a.getChild(1), P)
                    if num is None or den is None or den[0] != 0 or den[1] == 0:
                        raise _TD("% or quotient with a non-affine argument")
                    s, b = num[0] / den[1], num[1] / den[1]
                else:
                    af = self.affine(arg, P)
                    if af is None:
                        raise _TD("floor/ceiling with a non-affine argument")
                    s, b = af
                if s != 0:
                    k0, k1 = sorted(((s * lo + b), (s * self.t_end + b)))
                    n = math.floor(k1) - math.ceil(k0) + 1
                    if n > 200000:
                        raise _TD("too many breakpoints")
                    for k in range(math.ceil(k0), math.floor(k1) + 1):
                        brk.add((k - b) / s)
            if ty in rel:
                rels.append(a)
            for i in range(a.getNumChildren()):
                walk(a.getChild(i))
        for e in self.events:
            walk(e.getTrigger().getMath())
        pts = sorted({lo, self.t_end} | {x for x in brk if lo < x < self.t_end})
        roots = set()
        for x0, x1 in zip(pts, pts[1:]):                       # inside a segment each relation is affine in time
            if x1 - x0 <= 1e-12 * max(1.0, self.t_end):
                continue
            ta, tb = x0 + 0.25 * (x1 - x0), x0 + 0.75 * (x1 - x0)
            for r in rels:
                try:
                    ga = float(self.ev(r.getChild(0), ta, P)) - float(self.ev(r.getChild(1), ta, P))
                    gb = float(self.ev(r.getChild(0), tb, P)) - float(self.ev(r.getChild(1), tb, P))
                except (ZeroDivisionError, ValueError, OverflowError, TypeError):
                    continue
                if ga != gb:
                    tr = ta - ga * (tb - ta) / (gb - ga)
                    if x0 - 1e-12 * self.t_end <= tr <= x1 + 1e-12 * self.t_end:
                        roots.add(min(max(tr, x0), x1))
        return sorted({x for x in set(pts) | brk | roots if lo < x <= self.t_end})

    # simulation ------------------------------------------------------------------------
    def run(self, max_fire=20000):
        P = {v: self.consts.get(v, 0.0) for v in self.P}
        prev = [bool(e.getTrigger().getInitialValue()) if e.getTrigger().isSetInitialValue() else True
                for e in self.events]
        trig_vars = set()
        for e in self.events:
            trig_vars |= self.symbols(e.getTrigger().getMath())
        firings = []

        def step(t, probe):
            """evaluate every trigger at `probe`; fire (at time t) those that go false -> true, cascading over
            assignments made at the same time. Returns True if a variable used in a trigger changed."""
            nonlocal P
            moved = False
            for _ in range(100):
                fired = False
                for i, e in enumerate(self.events):
                    now = bool(self.ev(e.getTrigger().getMath(), probe, P))
                    if now and not prev[i]:
                        vals = {a.getVariable(): self.value(a, probe, P) for a in e.getListOfEventAssignments()}
                        firings.append((t, i, vals))
                        if len(firings) > max_fire:
                            raise _TD(f"more than {max_fire:,} event firings")
                        newP = dict(P)
                        for k, (kind, v) in vals.items():
                            if k in self.P:
                                newP[k] = v
                        if any(newP[k] != P[k] for k in P if k in trig_vars):
                            moved = True
                        P = newP
                        fired = True
                    prev[i] = now
                if not fired:
                    return moved
            raise _TD("events keep re-firing at the same time")

        # the initial values come from libRoadRunner after its t = 0 events have fired, so events that fire at
        # t = 0 are already applied: start from the trigger values at t = 0
        prev[:] = [bool(self.ev(e.getTrigger().getMath(), 0.0, P)) for e in self.events]
        t = 0.0
        while True:
            cand = self.critical(P, t)
            restart = False
            for k, tc in enumerate(cand):
                nxt = cand[k + 1] if k + 1 < len(cand) else tc
                moved = step(tc, tc)                              # the instant tc
                if not moved and nxt > tc:                        # the open interval after tc (fires 'at' tc)
                    moved = step(tc, 0.5 * (tc + nxt))
                if moved:
                    t, restart = tc, True
                    break
            if not restart:
                return firings

    def value(self, a, t, P):
        """('set', v) or ('add', d) for one assignment, evaluated at the firing time."""
        tgt, ast = a.getVariable(), a.getMath()
        syms = self.symbols(ast)
        if tgt in self.dynamic:
            if tgt in syms:
                T = self.T
                if ast.getType() in (T.AST_PLUS, T.AST_MINUS) and ast.getNumChildren() == 2:
                    a0, a1 = ast.getChild(0), ast.getChild(1)
                    if a0.getType() == T.AST_NAME and a0.getName() == tgt and tgt not in self.symbols(a1):
                        d = float(self.ev(a1, t, P))
                        return ("add", d if ast.getType() == T.AST_PLUS else -d)
                    if ast.getType() == T.AST_PLUS and a1.getType() == T.AST_NAME and a1.getName() == tgt \
                            and tgt not in self.symbols(a0):
                        return ("add", float(self.ev(a0, t, P)))
                raise _TD(f"value assigned to '{tgt}' depends on its own state")
            if any(s in self.dynamic for s in syms):
                raise _TD(f"value assigned to '{tgt}' depends on a continuously changing variable")
        return ("set", float(self.ev(ast, t, P)))


class Converter:
    def __init__(self, sbml_path, rr_init=None, alias=None, offsets=None):
        self.path = sbml_path
        self.offsets_req = dict(offsets or {})   # sid -> shift: variables the SBML takes below zero (2026-11-01)
        self.offset = {}                         # sid -> shift actually applied
        self.outname = {}                        # sid -> name of its OUTC[] output species
        self.alias = alias or {}   # sid -> short Kintecus name (used only when USER expressions are too long)
        self.long_lines = []       # (comment, ulen, variables) of USER expressions over USER_SAFE
        self.overlong = []         # (comment, variables) of USER lines refused for length (rate term or whole line)
        self.flags = []          # (severity, text)  severity: 'info' | 'partial' | 'blocker'
        self.lines = []          # (k, reaction, comment)
        self.ofiles = []         # dict(filename, rows)
        rd = libsbml.SBMLReader()
        self.doc = rd.readSBMLFromFile(sbml_path)
        if self.doc.getModel() is None:
            raise Unsupported("no <model>")
        cp = libsbml.ConversionProperties(); cp.addOption("expandFunctionDefinitions", True)
        if self.doc.convert(cp) != libsbml.LIBSBML_OPERATION_SUCCESS:
            self.flag("partial", "function definitions could not all be inlined")
        self.m = self.doc.getModel()
        self.init = rr_init or {}          # values at t=0 from roadrunner: id -> value
        self.time_used = False
        self.switches = {}                 # (t0, up) -> kintecus name

    def flag(self, sev, text):
        if (sev, text) not in self.flags:
            self.flags.append((sev, text))

    # -- setup ---------------------------------------------------------------
    def classify(self):
        m = self.m
        self.comp = {c.getId(): c for c in m.getListOfCompartments()}
        self.spec = {s.getId(): s for s in m.getListOfSpecies()}
        self.par = {p.getId(): p for p in m.getListOfParameters()}
        self.arule, self.rrule = {}, {}
        for r in m.getListOfRules():
            if r.isAlgebraic():
                self.flag("blocker", "algebraic rule not convertible (Kintecus has no algebraic constraints)")
            elif r.isAssignment():
                self.arule[r.getVariable()] = r
            elif r.isRate():
                self.rrule[r.getVariable()] = r
        for v in list(self.rrule):
            if v in self.comp:
                self.flag("blocker", f"rate rule on compartment '{v}' (variable volume) not convertible")
                del self.rrule[v]
        for c in self.comp:
            if c in self.arule:
                self.flag("partial", f"compartment '{c}' has an assignment rule; its t=0 size is used")
        if m.getNumConstraints():
            self.flag("info", "constraints ignored")
        for r in m.getListOfReactions():
            if r.isSetFast() and r.getFast():
                self.flag("info", f"fast reaction '{r.getId()}' integrated as normal kinetics")
        for s in m.getListOfSpecies():
            if s.isSetConversionFactor() or m.isSetConversionFactor():
                self.flag("partial", "conversion factors ignored")
                break

    def vol(self, cid):
        c = self.comp[cid]
        v = self.init.get(cid, c.getSize() if c.isSetSize() else 1.0)
        if v is None or not math.isfinite(v) or v == 0:
            v = 1.0
        return float(v)

    # names --------------------------------------------------------------------
    def make_names(self):
        self.kname, used = {}, set()

        def assign(sid, prefix=""):
            # Kintecus USER variables ('c' + species name) longer than 10 characters give FATAL #999 (probed on V2027),
            # so every Kintecus species name is kept to 9 characters.
            if sid in self.alias:
                used.add(self.alias[sid].upper()); self.kname[sid] = self.alias[sid]
                return self.alias[sid]
            base = prefix + sid
            clean = lambda nm: not any(b in ("C" + nm).upper() for b in BAD_SUB)
            ok = (NAME_RE.match(base) and len(base) <= NAME_MAX and base.upper() not in RESERVED
                  and base.upper() not in used and clean(base))
            name = base if ok else None
            if name is None:
                stem = re.sub(r"[^A-Za-z0-9]", "", base)
                guard = 0
                while any(x in ("C" + stem).upper() for x in BAD_SUB) and guard < 20:
                    full = "C" + stem
                    for x in BAD_SUB:
                        i = full.upper().find(x)
                        if i >= 0:
                            full = full[:i + 1] + "x" + full[i + 1:]
                            break
                    stem = full[1:] if full[0] == "C" and full[1] != "x" else "X" + full[2:]
                    guard += 1
                if not stem or not stem[0].isalpha():
                    stem = "X" + stem
                mt = re.match(r"^(.*?[A-Za-z])0*(\d+)$", stem)       # keep the id's trailing number
                if mt and len(mt.group(2)) < NAME_MAX:
                    tail = mt.group(2)
                    cand = mt.group(1)[:NAME_MAX - len(tail)] + tail
                    if cand.upper() not in used and cand.upper() not in RESERVED and clean(cand):
                        name = cand
                i = 1
                while name is None:
                    suf = str(i)
                    cand = stem[:NAME_MAX - len(suf)] + suf
                    if cand.upper() not in used and cand.upper() not in RESERVED and clean(cand):
                        name = cand
                    i += 1
            used.add(name.upper())
            self.kname[sid] = name
            return name
        for sid in self.spec:
            assign(sid)
        self._assign = assign

    def sym(self, sid):
        if sid in self.offset:                   # shifted state: the model's value is (state - shift)
            return OFFSET(sp.Symbol("c" + self.kname[sid]), sp.Float(self.offset[sid], 17))
        return sp.Symbol("c" + self.kname[sid])

    def apply_offsets(self):
        """A variable that the SBML takes below zero (a membrane potential, a deviation from basal, a species of a
        defective model) is integrated as state = value + shift >= 0; rates use (c<state> - shift) and an
        output-only species OUTC[state;-source*shift] gives the value itself (OUTC[] is not used in any rate).
        Kintecus stalls when it has to integrate a variable through zero, so this keeps every state positive."""
        if not self.offsets_req:
            return
        used = {n.upper() for n in self.kname.values()}
        for sid, shift in self.offsets_req.items():
            if sid not in self.kname or not (shift and shift > 0):
                continue
            if sid in self.arule:
                continue                     # defined by an assignment rule, not integrated
            if sid in self.changed_by_event:
                self.flag("info", f"'{sid}' goes below zero in the SBML but is changed by events; not shifted")
                continue
            if sid in self.spec and sid not in self.rrule:
                s = self.spec[sid]
                if s.getBoundaryCondition() or s.getConstant():
                    continue
            if sid in self.par and sid not in self.pseudo:
                continue
            out = self.kname[sid]
            base = out[:NAME_MAX - 2] + "_o"
            i = 1
            while base.upper() in used:
                base = out[:NAME_MAX - 2 - len(str(i))] + f"{i}_o"; i += 1
            used.add(base.upper())
            self.outname[sid] = out
            self.kname[sid] = base
            if sid in self.pseudo:
                self.pseudo[sid] = base
            self.offset[sid] = float(shift)
            self.flag("info", f"'{sid}' goes below zero in the SBML: integrated as {base} = value + {num(shift)}, "
                              f"output as {out} = OUTC[{base};-{num(shift)}*source]")

    # state determination --------------------------------------------------------
    def states(self):
        """Decide which SBML species / parameters become Kintecus species."""
        self.pseudo = {}
        self.events = self.parse_events()
        self.changed_by_event = {a for ev in self.events for a in ev["assign"]}
        self.pseudo = {}       # sbml id (parameter) -> kintecus name
        for pid in self.par:
            if pid in self.rrule or pid in self.changed_by_event:
                self.pseudo[pid] = self._assign(pid, "p_")
        # species in Kintecus = all species without assignment rule
        self.kspec = [s for s in self.spec if s not in self.arule]
        for s in self.spec:
            if s in self.arule:
                self.flag("info", f"species '{s}' is defined by an assignment rule; not a Kintecus state (not output)")

    # symbol resolution --------------------------------------------------------
    def resolver(self, local=None, depth=0):
        local = local or {}

        def res(name):
            if name == "__time__":
                self.time_used = True
                return sp.Symbol("cTIME_")
            if name in local:
                return sp.Float(local[name], 17)
            if name in self.arule:
                if depth > 50:
                    raise Unsupported("assignment rule recursion")
                # time conditions inside assignment rules become exact 0/1 switch species too (not abs() steps)
                return MathConv(self.resolver(local, depth + 1), self.switch if ARULE_SWITCH else None)(
                    self.arule[name].getMath())
            if name in self.spec:
                s = self.spec[name]
                varying = getattr(self, "varying", None)
                if varying is not None and name not in varying:
                    v0 = self.init.get(name)
                    if v0 is None:
                        v0 = (s.getInitialConcentration() if s.isSetInitialConcentration() else
                              s.getInitialAmount() / self.vol(s.getCompartment()) if s.isSetInitialAmount() else 0.0)
                    v0 = float(v0)
                    return sp.Float(v0 * self.vol(s.getCompartment()) if s.getHasOnlySubstanceUnits() else v0, 17)
                c = self.sym(name)
                return c * self.vol(s.getCompartment()) if s.getHasOnlySubstanceUnits() else c
            if name in self.pseudo:
                return self.sym(name)
            if name in self.par:
                v = self.init.get(name, self.par[name].getValue())
                if v is None or not math.isfinite(float(v)):
                    raise Unsupported(f"parameter '{name}' has no finite value")
                return sp.Float(float(v), 17)
            if name in self.comp:
                return sp.Float(self.vol(name), 17)
            for r in self.m.getListOfReactions():          # reaction id used as a symbol = its rate (kinetic law)
                if r.getId() == name:
                    kl = r.getKineticLaw()
                    if kl is None or kl.getMath() is None:
                        raise Unsupported(f"reaction id '{name}' used as a variable (reaction has no kinetic law)")
                    if depth > 50:
                        raise Unsupported("reaction rate recursion")
                    loc = {}
                    for lp in list(kl.getListOfParameters()) + list(getattr(kl, "getListOfLocalParameters", lambda: [])()):
                        if lp.isSetValue():
                            loc[lp.getId()] = lp.getValue()
                    return MathConv(self.resolver(loc, depth + 1))(kl.getMath())
            raise Unsupported(f"unknown symbol '{name}'")
        return res

    # events --------------------------------------------------------------------
    def parse_events(self):
        if UNROLL and self.m.getNumEvents():
            plain = True
            for ev in self.m.getListOfEvents():
                try:
                    self.time_trigger(ev.getTrigger().getMath() if ev.isSetTrigger() else None)
                except Unsupported:
                    plain = False
            if not plain:
                try:
                    return self.unroll_events()
                except _TD as ex:
                    self.flag("info", f"events not unrolled as time events ({ex})")
        out = []
        for ev in self.m.getListOfEvents():
            eid = ev.getId() or "event"
            try:
                t0 = self.time_trigger(ev.getTrigger().getMath() if ev.isSetTrigger() else None)
                if ev.isSetDelay():
                    d = MathConv(self.resolver())(ev.getDelay().getMath())
                    if not d.is_Number:
                        raise Unsupported("state-dependent delay")
                    t0 += float(d)
                assign = {}
                for a in ev.getListOfEventAssignments():
                    assign[a.getVariable()] = a.getMath()
                out.append({"id": eid, "t": t0, "assign": assign})
            except Unsupported as ex:
                self.flag("blocker", f"event '{eid}' not converted ({ex})")
        return out

    def unroll_events(self):
        """All events time-determined (see EventUnroller): their firings become plain time events."""
        dyn = {r.getVariable() for r in self.m.getListOfRules() if r.isRate() or r.isAssignment()}
        for r in self.m.getListOfReactions():
            for sr in list(r.getListOfReactants()) + list(r.getListOfProducts()):
                sp_ = self.m.getSpecies(sr.getSpecies())
                if sp_ is not None and not sp_.getBoundaryCondition() and not sp_.getConstant():
                    dyn.add(sr.getSpecies())
        U = EventUnroller(self.m, self.init, self.t_end, dyn)
        fir = U.run()
        out = []
        for t, i, vals in fir:
            eid = self.m.getEvent(i).getId() or "event"
            assign = {}
            for tgt, (kind, v) in vals.items():
                assign[tgt] = libsbml.parseL3Formula(repr(v) if kind == "set" else f"{tgt} + ({repr(v)})")
            out.append({"id": eid, "t": t, "assign": assign})
        ne = self.m.getNumEvents()
        self.flag("info", f"{ne} event(s) depend only on time (and on variables only they change): unrolled into "
                          f"{len(out)} firing(s) up to t = {num(self.t_end)}")
        return out

    def time_trigger(self, ast):
        """Accept triggers of the form time >= t0 / time > t0 / t0 <= time (optionally inside and(...) with true)."""
        if ast is None:
            raise Unsupported("no trigger")
        T = libsbml
        t = ast.getType()
        if t in (T.AST_RELATIONAL_GEQ, T.AST_RELATIONAL_GT, T.AST_RELATIONAL_LEQ, T.AST_RELATIONAL_LT) \
                and ast.getNumChildren() == 2:
            a, b = ast.getChild(0), ast.getChild(1)
            lhs_time = a.getType() == T.AST_NAME_TIME
            rhs_time = b.getType() == T.AST_NAME_TIME
            if lhs_time ^ rhs_time:
                other = b if lhs_time else a
                ev_t = {x.getVariable() for e in self.m.getListOfEvents() for x in e.getListOfEventAssignments()}
                names = set()
                stack = [other]
                while stack:
                    n_ = stack.pop()
                    if n_.getType() == T.AST_NAME:
                        names.add(n_.getName())
                    stack.extend(n_.getChild(i) for i in range(n_.getNumChildren()))
                if names & ev_t:      # e.g. time >= t_int*cups + t_lag with cups = cups + 1: re-arms, not one firing
                    raise Unsupported("trigger threshold is changed by events")
                v = MathConv(self.resolver())(other)
                if not v.is_Number:
                    raise Unsupported("trigger threshold is not constant")
                forward = (t in (T.AST_RELATIONAL_GEQ, T.AST_RELATIONAL_GT)) == lhs_time
                if not forward:
                    raise Unsupported("trigger fires while time is below a threshold")
                return float(v)
        raise Unsupported("state-triggered (not a time >= t0 trigger)")

    def switch(self, t0, up):
        key = (round(t0, 12), up)
        if key not in self.switches:
            self.switches[key] = f"SW{len(self.switches) + 1}"
        self.time_used = True                       # keeps the time pseudo-species out of the way of names
        return sp.Symbol("c" + self.switches[key])

    # reactions -----------------------------------------------------------------
    def split_terms(self, v):
        if sp.count_ops(v) > 150:                  # expansion of large nested laws can explode
            return list(sp.Add.make_args(v)) if len(sp.Add.make_args(v)) <= 40 else [v]
        try:
            ex = sp.expand_mul(v)
        except Exception:
            ex = v
        terms = list(sp.Add.make_args(ex))
        if len(terms) > 40:
            return [v]
        if max(len(str(t)) for t in terms) > USER_MAX_LEN and len(str(ex)) > 2 * len(str(v)):
            return [v]          # expand_mul multiplied a factored law out into a term too long for a line
        return terms

    def emit(self, term, lhs, rhs, comment, allowed=frozenset()):
        """lhs/rhs: gross (kintecus name, stoich) lists. term: rate contribution in conc/time, already divided by
        the group volume. Species in `allowed` that multiply the rate are moved onto both sides so Kintecus's
        Jacobian sees the dependence (a zero-order line hides it, which breaks stiff semi-implicit steps)."""
        c, rest = term.as_coeff_Mul()
        if c == 0:
            return
        if c < 0:
            lhs, rhs, term = rhs, lhs, -term
        offs = {self.kname[s] for s in self.offset}
        if offs and any(nm in offs for nm, _ in list(lhs) + list(rhs)):
            # shifted states change only through signed source lines (their 'concentration' is value + shift)
            nf = lambda n: sp.Integer(1) if abs(float(n) - 1) < 1e-12 else sp.Float(float(n), 17)
            for nm, st in lhs:
                if nm in offs:
                    self._line(-nf(st) * term, [], nm, f"{comment} (shifted {nm})")
            for nm, st in rhs:
                if nm in offs:
                    self._line(nf(st) * term, [], nm, f"{comment} (shifted {nm})")
            lhs = [x for x in lhs if x[0] not in offs]
            rhs = [x for x in rhs if x[0] not in offs]
            if not lhs and not rhs:
                return
            if not lhs:                          # only products left: zero order from source
                R = "+".join(nm if abs(float(n) - 1) < 1e-12 else f"{coef(float(n))}{nm}" for nm, n in rhs)
                self._line(term, [], R, comment)
                return
        side_l, side_r = {}, {}
        frac = {}                                # non-integer reactant coefficients, written as Kintecus coefficients
        denom = sp.Integer(1)
        tot = {}
        for nm, st in lhs:
            tot[nm] = tot.get(nm, 0.0) + float(st)
        for nm, n in tot.items():
            if abs(n - round(n)) <= 1e-12:
                side_l[nm] = int(round(n))
                denom *= sp.Symbol("c" + nm) ** int(round(n))
            elif n <= FRAC_MAX:
                # Kintecus: a line 'nA+...' has rate = (k or USER) * [A]^n and consumes n*rate of A (probed with
                # 0.25, 0.5, 0.7, 1.5, 4.33 and 0.0009967), so the rate is divided by [A]^n like an integer order
                frac[nm] = n
                denom *= sp.Symbol("c" + nm) ** sp.Float(n, 17)
                self.flag("info", f"{comment}: non-integer reactant coefficient {num(n)} on {nm} written as a "
                                  f"Kintecus coefficient")
            else:
                # large non-integer coefficient: [A]^n in the rate would over/underflow; consume A on its own line
                # at n*rate (A + S[USER; n*rate/cA] ==> sink) and keep A out of the main line
                self.flag("info", f"{comment}: non-integer reactant coefficient {num(n)} on {nm} consumed on a "
                                  f"separate line (n > {FRAC_MAX})")
                self.emit(term * sp.Float(n, 17), [(nm, 1)], [], f"{comment} consumption of {nm}", allowed)
        for nm, st in rhs:
            side_r[nm] = side_r.get(nm, 0) + float(st)
        user = term / denom                      # Mul flattening cancels monomial factors automatically
        inden = lambda u: any(sp.Symbol("c" + n) in u.as_numer_denom()[1].free_symbols for n in list(side_l) + list(frac))
        if (side_l or frac) and inden(user):
            try:
                if len(str(user)) < 600:
                    user = safe_cancel(user)
            except Exception:
                pass
            if inden(user) and DIV_MODE.startswith("split") and self._singular(user, list(side_l) + list(frac)):
                self._decoupled(term, side_l, frac, side_r, comment)
                return
            if inden(user) and SIGNED_NONPROP and self._singular(user, list(side_l) + list(frac)):
                self._signed(term, side_l, frac, side_r, comment)
                return
            if inden(user):
                num_, den_ = user.as_numer_denom()
                syms = {sp.Symbol("c" + n) for n in list(side_l) + list(frac)} & den_.free_symbols
                if DIV_MODE == "abs" and self._singular(user, list(side_l) + list(frac)):
                    # smooth on/off: rate = f*[A]/(|[A]|+eps_A), eps_A = DIV_REL x the species' scale; no pole for
                    # [A] < 0 (it turns into a small restoring production), and a switch width the integrator can
                    # resolve (1e-30 is below double-precision resolution of any concentration -> a step function)
                    sing = {x for x in syms if den_.subs(x, 0) == 0}
                    user = num_ / den_.subs({x: (sp.Abs(x) + self._div_eps(x.name[1:])) if x in sing else x + DIV_EPS
                                             for x in syms})
                else:
                    user = num_ / den_.subs({x: x + DIV_EPS for x in syms})
                self.flag("info", f"{comment}: rate is not proportional to its reactants; divided form used with "
                                  f"1/(c+{num(DIV_EPS)}) so an empty species is not consumed")
        # move multiplicative modifier factors onto both sides
        moved = []
        if not user.is_Number:
            for f in sp.Mul.make_args(user):
                b, e = (f.base, f.exp) if f.is_Pow else (f, sp.Integer(1))
                if b.is_Symbol and b.name.startswith("c") and b.name[1:] in allowed and e.is_Integer and 0 < e <= 4:
                    nm = b.name[1:]
                    side_l[nm] = side_l.get(nm, 0) + int(e)
                    side_r[nm] = side_r.get(nm, 0) + int(e)
                    user = user / f
                    moved.append((nm, int(e), f))
        # Kintecus allows at most MAX_REACTANTS reactant entries per line (FATAL #39 otherwise): first put moved
        # modifiers back into the rate expression, then split the consumed reactants over several lines
        while sum(side_l.values()) + len(frac) > MAX_REACTANTS and moved:
            nm, e, f = moved.pop()
            side_l[nm] -= e; side_r[nm] -= e
            if side_l[nm] == 0: del side_l[nm]
            user = user * f
        R = "+".join((nm if abs(n - 1) < 1e-12 else f"{coef(n)}{nm}") for nm, n in side_r.items() if n > 1e-12) or SINK
        # reactant entries: (text, species, order); an integer coefficient is written as repeated names (A+A)
        ents = [(nm, nm, 1) for nm, n in side_l.items() for _ in range(n)] + \
               [(f"{coef(n)}{nm}", nm, n) for nm, n in frac.items()]
        reac = [e[0] for e in ents]
        if len(reac) <= MAX_REACTANTS:
            self._line(user, reac, R, comment)
            return
        chunks = [ents[i:i + MAX_REACTANTS] for i in range(0, len(ents), MAX_REACTANTS)]
        self.flag("info", f"{comment}: {len(reac)} reactants split over {len(chunks)} lines "
                          f"(Kintecus allows {MAX_REACTANTS} per line)")
        for j, ch in enumerate(chunks):
            rest = list(ents)
            for e in ch:
                rest.remove(e)
            u = user
            for _, nm, o in rest:
                u = u * (sp.Symbol("c" + nm) if o == 1 else sp.Symbol("c" + nm) ** sp.Float(o, 17))
            self._line(u, [e[0] for e in ch], R if j == 0 else SINK, f"{comment} part {j+1}/{len(chunks)}")

    def _div_eps(self, nm):
        """switch width for the divided form: DIV_REL x |initial value| (or x the median nonzero initial value)."""
        if not hasattr(self, "_eps_cache"):
            self._eps_cache = {}
            vals = [abs(float(v)) for v in self.init.values() if v and math.isfinite(float(v))]
            vals.sort()
            self._eps_med = vals[len(vals) // 2] if vals else 1.0
            self._k2sid = {n: sid for sid, n in self.kname.items()}
        if nm not in self._eps_cache:
            v = self.init.get(self._k2sid.get(nm, ""), 0.0) or 0.0
            try:
                v = abs(float(v))
            except (TypeError, ValueError):
                v = 0.0
            base = v if v > 0 and math.isfinite(v) else self._eps_med
            self._eps_cache[nm] = sp.Float(float(f'{max(1e-30, DIV_REL * base):.3g}'))
        return self._eps_cache[nm]

    @staticmethod
    def _singular(user, names):
        """True if the USER factor blows up at [A] = 0 for a reactant A (the rate is not proportional to [A]);
        a denominator like (A + Km) is not singular."""
        _, den_ = user.as_numer_denom()
        for nm in names:
            s_ = sp.Symbol("c" + nm)
            if s_ in den_.free_symbols:
                try:
                    if den_.subs(s_, 0) == 0:            # plain substitution: sympy folds c*(...) and sums of c-terms to 0
                        return True
                except Exception:
                    return True
        return False

    def _decoupled(self, term, side_l, frac, side_r, comment):
        """The rate's order in some reactant differs from its stoichiometry (e.g. rate k*[II_f]*[LIPID] with 99 LIPID
        consumed, or a zero-order drain). Kintecus ties a reactant's order to its coefficient, so such reactants are
        consumed on their own lines; the others stay on the main line with the products:
          main line:  B+...+S[USER;rate/[B]^nB...]==>P      (reactants whose order matches; 'source' if none)
          A, rate ~ [A]:          A+S[USER;nA*rate/cA]==>sink        (first order: no division singularity)
          A, rate not ~ [A]:      source+S[USER;-nA*rate]==>A       (zero order: signed, exactly as the SBML)"""
        def canc(u):
            try:
                return safe_cancel(u) if len(str(u)) < 600 else u
            except Exception:
                return u
        reac = [(nm, n, True) for nm, n in side_l.items()] + [(nm, n, False) for nm, n in frac.items()]
        keep, sep, u = [], [], term
        for nm, n, is_int in reac:
            if nm in side_r:                    # catalyst-like species on both sides: keep the old handling
                sep.append((nm, n)); continue
            c_ = sp.Symbol("c" + nm)
            cand = canc(u / (c_ ** (int(n) if is_int else sp.Float(n, 17))))
            if not self._singular(cand, [nm]):
                keep.append((nm, n, is_int)); u = cand
            else:
                sep.append((nm, n))
        nf = lambda n: sp.Integer(1) if abs(n - 1) < 1e-12 else sp.Float(n, 17)
        net_r = dict(side_r)
        for nm, n in sep:                        # a separated species that is also a product: net it here
            if nm in net_r:
                net_r[nm] -= n
        R = "+".join(nm if abs(n - 1) < 1e-12 else f"{coef(n)}{nm}" for nm, n in net_r.items() if n > 1e-12)
        sep = [(nm, n - side_r[nm] if nm in side_r else n) for nm, n in sep]
        sep = [(nm, n) for nm, n in sep if n > 1e-12]
        ents = [nm for nm, n, i in keep if i for _ in range(int(n))] + [f"{coef(n)}{nm}" for nm, n, i in keep if not i]
        if R or ents:
            if len(ents) <= MAX_REACTANTS:
                self._line(u, ents, R or SINK, comment)
            else:                                # too many entries: fall back to the products on a source line
                self._line(term, [], R or SINK, comment)
                sep += [(nm, n) for nm, n, i in keep]
        kinds = []
        for nm, n in sep:
            c_ = sp.Symbol("c" + nm)
            u1 = canc(term / c_)
            if not self._singular(u1, [nm]):
                self._line(nf(n) * u1, [nm], SINK, f"{comment} (consumes {nm})"); kinds.append(f"{nm}: first order")
            else:
                self._line(-nf(n) * term, [], nm, f"{comment} (consumes {nm})"); kinds.append(f"{nm}: zero order, signed")
        self.flag("info", f"{comment}: rate order differs from the stoichiometry; consumed on separate line(s) - "
                          + "; ".join(kinds))

    def _any_singular(self, terms, lhs, rhs):
        """number of terms that are not proportional to the species they consume (all of them: lin-log kinetics)"""
        n_sing = 0
        for tm in terms:
            c, _ = tm.as_coeff_Mul()
            side = rhs if c < 0 else lhs
            den = sp.Integer(1)
            for nm, st in side:
                den *= sp.Symbol("c" + nm) ** (int(round(float(st))) if abs(float(st) - round(float(st))) < 1e-12
                                                else sp.Float(float(st), 17))
            u = tm / den
            try:
                u = safe_cancel(u) if len(str(u)) < 600 else u
            except Exception:
                pass
            if self._singular(u, [nm for nm, _ in side]):
                n_sing += 1
        return n_sing

    def _whole_signed(self, rate, lhs, rhs, comment):
        net = {}
        for nm, st in lhs:
            net[nm] = net.get(nm, 0.0) - float(st)
        for nm, st in rhs:
            net[nm] = net.get(nm, 0.0) + float(st)
        prods = [(nm, n) for nm, n in net.items() if n > 1e-12]
        cons = [(nm, -n) for nm, n in net.items() if n < -1e-12]
        if prods:
            self._line(rate, [], "+".join(nm if abs(n - 1) < 1e-12 else f"{coef(n)}{nm}" for nm, n in prods), comment)
        for nm, n in cons:
            self._line(-rate * (sp.Integer(1) if abs(n - 1) < 1e-12 else sp.Float(n, 17)), [], nm,
                       f"{comment} (consumes {nm})")
        self.flag("info", f"{comment}: rate not proportional to its reactants; written once as signed lines "
                          f"(source+S[USER;...]) for " + ", ".join(nm for nm, _ in cons))

    def _signed(self, term, side_l, frac, side_r, comment):
        """A rate term that is not proportional to its reactants, written exactly as in the SBML: the products are
        made by 'source+S[USER;term]==>P', and each consumed species A by 'source+S[USER;-n*term]==>A' (Kintecus takes
        negative rates). No division by [A], so no near-discontinuity at [A] ~ 0, and A may go below zero exactly
        where the SBML does."""
        net = {}
        for nm, n in side_l.items():
            net[nm] = net.get(nm, 0) - n
        for nm, n in frac.items():
            net[nm] = net.get(nm, 0) - n
        for nm, n in side_r.items():
            net[nm] = net.get(nm, 0) + n
        prods = [(nm, n) for nm, n in net.items() if n > 1e-12]
        cons = [(nm, -n) for nm, n in net.items() if n < -1e-12]
        self.flag("info", f"{comment}: rate is not proportional to its reactants; written as signed zero-order lines "
                          f"(source+S[USER;...]) for " + ", ".join(nm for nm, _ in cons))
        if prods:
            R = "+".join(nm if abs(n - 1) < 1e-12 else f"{coef(n)}{nm}" for nm, n in prods)
            self._line(term, [], R, comment)
        for nm, n in cons:
            self._line(-term * (sp.Integer(1) if abs(n - 1) < 1e-12 else sp.Float(n, 17)), [], nm,
                       f"{comment} (consumes {nm})")

    def _line(self, user, reac, R, comment):
        L = "+".join(reac) or SOURCE
        if user.is_Number or (not user.free_symbols and user.is_number):
            # a constant rate, also when written with functions (e.g. 0.4*log(2)): Kintecus needs a number here -
            # a USER expression with no variables (S[USER;:(...)]) stops Kintecus with FATAL #999
            v = float(sp.N(user, 20)) if not user.is_Number else float(user)
            if not math.isfinite(v):
                raise Unsupported(f"rate constant is not finite ({user})")
            if v == 0 or (v < 0 and reac):
                return
            self.lines.append((num(v), f"{L}==>{R}", comment))     # v < 0 only on a signed source line
            return
        vars_ = sorted({sym.name for sym in user.free_symbols})
        if len(vars_) > USER_MAX_VARS:
            raise Unsupported(f"rate term needs {len(vars_)} variables (> {USER_MAX_VARS})")
        body = kprint(user)
        if ulen(body) > USER_SAFE:
            # long expression: shortest form that fits; otherwise the shortest '^'-free form (long '^' expressions
            # are mis-evaluated silently; long '^'-free ones either evaluate correctly or stop with FATAL #999)
            # exact digits first ('^'-free forms are read correctly at any length, and the per-line evaluation check
            # catches the rest); rounded constants only when nothing exact fits - rounding can move a near-cancelling
            # denominator onto zero
            pref = [kprint(user, 17, True), kprint(user, 10, True), kprint(user, 10), kprint(user, 7, True),
                    kprint(user, 7)]
            for cand in pref:
                if ("^" not in cand or ulen(cand) <= USER_SAFE) and len(cand) <= USER_MAX_LEN:
                    body = cand
                    break
            else:
                body = min([body] + pref, key=len)
            if ulen(body) > USER_SAFE:
                self.long_lines.append((comment, ulen(body), vars_))
                self.flag("info", f"{comment}: USER expression is {len(body)} characters; Kintecus V2027 mis-evaluates "
                                  f"USER expressions longer than about {USER_SAFE} characters (see KINTECUS_ISSUES)")
        if len(body) > USER_MAX_LEN:
            parts = split_rate(user)
            if SPLIT_LONG and len(parts) > 1:
                # rates add: the same reactants and products on several lines, each carrying part of the rate
                chunks, cur = [], []
                for t in parts:
                    trial = cur + [t]
                    if cur and len(kprint(sp.Add(*trial, evaluate=False))) > USER_MAX_LEN - 40:
                        chunks.append(sp.Add(*cur, evaluate=False)); cur = [t]
                    else:
                        cur = trial
                if cur:
                    chunks.append(sp.Add(*cur, evaluate=False))
                if len(chunks) > 1:
                    n0, nl0 = len(self.lines), len(self.long_lines)
                    try:
                        for k, ch in enumerate(chunks, 1):
                            self._line(ch, reac, R, f"{comment} (part {k}/{len(chunks)})")
                    except Unsupported:
                        del self.lines[n0:]; del self.long_lines[nl0:]      # all or nothing
                        raise
                    self.flag("info", f"{comment}: rate is {len(body)} characters (> {USER_MAX_LEN}); written as "
                                      f"{len(chunks)} lines with the same reactants and products, each with part "
                                      f"of the rate (the parts add up to the rate)")
                    return
            self.overlong.append((comment, vars_))     # convert() retries with short names
            raise Unsupported(f"rate term is {len(body)} characters (> {USER_MAX_LEN})")
        line = f"{L}+S[USER;{','.join(vars_)}:({body})]==>{R}"
        if len(line) > LINE_MAX:
            self.overlong.append((comment, vars_))
            raise Unsupported(f"reaction line is {len(line)} characters")
        self.lines.append(("1", line, comment))

    def stoich(self, r):
        """Gross stoichiometry {sid: n} for reactants and products."""
        out = []
        for lst in (r.getListOfReactants(), r.getListOfProducts()):
            d = {}
            for sr in lst:
                sid = sr.getSpecies()
                if sr.isSetStoichiometryMath():
                    e = MathConv(self.resolver())(sr.getStoichiometryMath().getMath())
                    if not e.is_Number:
                        raise Unsupported("non-constant stoichiometryMath")
                    st = float(e)
                else:
                    st = sr.getStoichiometry() if (sr.isSetStoichiometry() or self.m.getLevel() < 3) else 1.0
                    if sr.isSetId() and (sr.getId() in self.arule or sr.getId() in self.rrule):
                        raise Unsupported("variable stoichiometry")
                d[sid] = d.get(sid, 0.0) + st
            out.append(d)
        return out

    def is_dyn(self, sid):
        return sid in self.kspec and not self.spec[sid].getBoundaryCondition() and not self.spec[sid].getConstant()

    def convert_reaction(self, r):
        rid = r.getId()
        if not r.isSetKineticLaw() or not r.getKineticLaw().isSetMath():
            raise Unsupported("no kinetic law")
        kl = r.getKineticLaw()
        local = {}
        for p in list(kl.getListOfParameters()) + list(getattr(kl, "getListOfLocalParameters", lambda: [])()):
            local[p.getId()] = p.getValue()
        mc = MathConv(self.resolver(local), self.switch)
        v = mc(kl.getMath())            # substance / time
        if mc.used_step:
            self.flag("info", f"{rid}: piecewise/conditional rate encoded with abs()-based step functions")
        reac, prod = self.stoich(r)
        reac = {s: n for s, n in reac.items() if self.is_dyn(s)}
        prod = {s: n for s, n in prod.items() if self.is_dyn(s)}
        net = {s: prod.get(s, 0) - reac.get(s, 0) for s in set(reac) | set(prod)}
        if not any(abs(n) > 0 for n in net.values()):
            return
        vol = lambda s: self.vol(self.spec[s].getCompartment())
        groups = {}
        for s in sorted(set(reac) | set(prod)):          # sorted: the same deck on every run (no hash order)
            groups.setdefault(vol(s), []).append(s)
        terms = self.split_terms(v)
        for V, members in groups.items():
            if not any(abs(net[s]) > 0 for s in members):
                continue
            lhs = [(self.kname[s], reac[s]) for s in members if s in reac]
            rhs = [(self.kname[s], prod[s]) for s in members if s in prod]
            allowed = frozenset(self.kname[s] for s in self.reacting_ids if vol(s) == V)
            if DIV_MODE.startswith("split") and len(terms) > 1 and self._any_singular(terms, lhs, rhs) == len(terms):
                # e.g. lin-log kinetics J*(1 + e1*log(S/S0) + ...): no term is proportional to the species it
                # consumes; write the whole rate once, signed (1 + number of species lines, not one per term)
                self._whole_signed(v / sp.Float(V, 17), lhs, rhs,
                                   f"{rid}" + (f" (volume {num(V)})" if len(groups) > 1 else ""))
                continue
            for i, tm in enumerate(terms):
                self.emit(tm / sp.Float(V, 17), lhs, rhs,
                          f"{rid}" + (f" term {i+1}/{len(terms)}" if len(terms) > 1 else "") +
                          (f" (volume {num(V)})" if len(groups) > 1 else ""), allowed)

    def convert_rate_rules(self):
        for var, rule in self.rrule.items():
            try:
                f = MathConv(self.resolver(), self.switch)(rule.getMath())
                if var in self.spec:
                    s = self.spec[var]
                    if s.getHasOnlySubstanceUnits():
                        f = f / self.vol(s.getCompartment())
                    name = self.kname[var]
                elif var in self.pseudo:
                    name = self.kname[var]
                else:
                    raise Unsupported("rate rule target unknown")
                vol_ok = frozenset(self.kname[x] for x in self.reacting_ids)
                for i, tm in enumerate(self.split_terms(f)):
                    self.emit(tm, [], [(name, 1)], f"rate rule d{var}/dt term {i+1}", vol_ok)
            except Unsupported as ex:
                self.flag("blocker", f"rate rule for '{var}' not converted ({ex})")

    # profiles --------------------------------------------------------------------
    @staticmethod
    def grid(times):
        fr = [Fraction(t).limit_denominator(10 ** 6) for t in times if t > 0]
        if not fr:
            return None
        g = fr[0]
        for f in fr[1:]:
            g = Fraction(math.gcd(g.numerator * f.denominator, f.numerator * g.denominator),
                         g.denominator * f.denominator)
        return float(g)

    @staticmethod
    def pair_rows(pairs):
        """Constant File as Time/value pairs (Kintecus build 2026-11-03): header line with 'Time', pairs, END."""
        return [["Time(s)", "Concentrations"]] + [[num(t), num(v)] for t, v in pairs] + [["END", ""]]

    @staticmethod
    def is_pair_file(o):
        return bool(o["rows"]) and "TIME" in str(o["rows"][0][0]).upper()

    def use_pairs(self):
        """time/value pairs for the next profile: pair format chosen and fewer than PAIR_FILES_MAX pair files so far
        (Kintecus holds 50 pair files); later profiles fall back to the spacing format."""
        return PROFILE_FORMAT == "pairs" and sum(1 for o in self.ofiles if self.is_pair_file(o)) < PAIR_FILES_MAX

    def pair_profile(self, tgt, evs, reacting, incs, sets, name, fn, t_end, col7):
        """Time-event changes of one target as Time/value pairs: -ADD (doses on a reacting species) or -SET
        (absolute values; also every non-reacting species/parameter, whose initial value stays in SPECIES.DAT)."""
        evs = [e for e in evs if e[0] <= t_end]
        if reacting and incs:
            acc = {}
            for t, d in incs:
                if t <= t_end:
                    acc[t] = acc.get(t, 0.0) + d
            pairs, kind = [(t, d) for t, d in sorted(acc.items()) if d != 0], "-ADD"
        elif reacting:
            last = {}
            for t, v in sorted(sets):
                if t <= t_end:
                    last[t] = v
            pairs, kind = sorted(last.items()), "-SET"
        else:
            cur = self.init.get(tgt)
            last = {}
            for t, ast, _ in evs:
                e = MathConv(self.resolver())(ast)
                if not e.is_Number:
                    e = e.subs(sp.Symbol("c" + name), cur)
                    if not e.is_Number:
                        raise Unsupported("assigned value depends on dynamic state")
                cur = float(e)
                if tgt in self.spec and self.spec[tgt].getHasOnlySubstanceUnits():
                    cur /= self.vol(self.spec[tgt].getCompartment())
                last[t] = cur
            pairs, kind = sorted(last.items()), "-SET"
        if not pairs:
            self.flag("info", f"events leave '{tgt}' unchanged; no profile written")
            return True
        if len(pairs) > PAIR_MAX:
            return False                    # more than Kintecus holds: the caller writes the spacing format
        # a pair at t = 0 is applied at the start (time > 0 events); 11-03 build: a pair between 0 and the starting
        # integration time (1e-12, say) is never applied
        pairs = [(max(t, 0.0), v) for t, v in pairs]
        self.ofiles.append({"filename": fn, "rows": self.pair_rows(pairs)})
        col7[tgt] = fn + kind
        self.flag("info", f"time event(s) on '{tgt}' encoded as time/value pairs in {fn} ({kind[1:]})")
        return True

    def build_events(self, t_end):
        """Return {sbml_target: (filename_on_disk, col7_reference)} and write profile rows."""
        col7 = {}
        by_target = {}
        for ev in self.events:
            for tgt, ast in ev["assign"].items():
                by_target.setdefault(tgt, []).append((ev["t"], ast, ev["id"]))
        for tgt, evs in by_target.items():
            try:
                if tgt in self.comp:
                    raise Unsupported("assignment to a compartment")
                if tgt not in self.spec and tgt not in self.pseudo:
                    raise Unsupported("target is not a species or parameter")
                evs.sort(key=lambda x: x[0])
                times = [t for t, _, _ in evs if t <= t_end]
                if not times:
                    continue
                # rate-rule variables are Kintecus species integrated by their own lines: dosed like reacting ones
                reacting = (tgt in self.spec and tgt in self.reacting) or tgt in self.rrule
                name = self.kname[tgt]
                fn = f"EV_{name}"[:20].upper() + ".TXT"
                if reacting:
                    sym = sp.Symbol("c" + name)
                    incs, sets = [], []
                    for t, ast, _ in evs:
                        e = MathConv(self.resolver())(ast)
                        if tgt in self.spec and self.spec[tgt].getHasOnlySubstanceUnits():
                            e = e / self.vol(self.spec[tgt].getCompartment())
                        d = sp.simplify(e - sym)
                        if not d.is_Number and d.free_symbols == {sym}:
                            # X + dose written through volumes: (dose/V + X/V)*V leaves X*(1 - 1.0000000000000002)
                            dd = sp.expand(d)
                            k1 = dd.coeff(sym, 1)
                            if k1.is_Number and abs(float(k1)) < 1e-12 and not dd.coeff(sym, 0).free_symbols:
                                d = dd.coeff(sym, 0)
                        if d.is_Number:
                            incs.append((t, float(d)))
                        elif e.is_Number:
                            sets.append((t, float(e)))
                        else:
                            raise Unsupported("assigned value depends on other dynamic species")
                    if sets and incs:
                        raise Unsupported("mixes absolute and incremental assignments")
                    if incs and all(d == 0 for _, d in incs):
                        self.flag("info", f"events add 0 to '{tgt}'; no profile written")
                        continue
                if self.use_pairs() and self.pair_profile(tgt, evs, reacting, incs if reacting else [],
                                                          sets if reacting else [], name, fn, t_end, col7):
                    continue
                dt = self.grid(times + [t_end])
                if dt is not None and t_end / dt + 2 <= PROF_MAX:
                    n = int(round(t_end / dt)) + 1
                else:
                    # the profile runs once (negative spacing) and then holds its last value, so it only has to
                    # reach just past the last event: spacing from the event times alone (pinned and -ADD profiles;
                    # -ADD ends with a 0 entry, so the re-applied last entry adds nothing)
                    dt = self.grid(times)
                    if dt is None or max(times) / dt + 3 > PROF_MAX:
                        raise Unsupported(f"event times need more than {PROF_MAX:,} profile points")
                    n = int(round(max(times) / dt)) + 2
                if not reacting:
                    # pinned piecewise-constant profile: value k applies on ((k-1)dt, k dt]
                    cur = self.init.get(tgt)
                    vals, ei = [], 0
                    for k in range(1, n + 1):
                        lo = (k - 1) * dt
                        while ei < len(evs) and evs[ei][0] <= lo + 1e-12 * max(1, t_end):
                            e = MathConv(self.resolver())(evs[ei][1])
                            if not e.is_Number:
                                e = e.subs(sp.Symbol("c" + name), cur)
                                if not e.is_Number:
                                    raise Unsupported("assigned value depends on dynamic state")
                            cur = float(e)
                            if tgt in self.spec and self.spec[tgt].getHasOnlySubstanceUnits():
                                cur /= self.vol(self.spec[tgt].getCompartment())
                            ei += 1
                        vals.append(cur)
                    self.ofiles.append({"filename": fn, "rows": [[num(-dt)]] + [[num(v)] for v in vals]})
                    col7[tgt] = fn
                else:
                    # reacting species: -ADD for increments (X + c), -SET when every grid point is an event
                    if incs:
                        vals = [0.0] * n
                        for t, d in incs:
                            vals[int(round(t / dt))] += d
                        if vals[-1] != 0:
                            vals.append(0.0)       # Kintecus keeps re-applying the last -ADD entry after the profile ends
                        self.ofiles.append({"filename": fn, "rows": [[num(-dt)]] + [[num(v)] for v in vals]})
                        col7[tgt] = fn + "-ADD"
                    else:
                        sets.sort()
                        if len(sets) == 1 and sets[0][0] > 0:
                            g, vals = sets[0][0], [self.init.get(tgt, 0.0), sets[0][1]]
                        else:
                            g = self.grid([t for t, _ in sets])
                            ks = [int(round(t / g)) for t, _ in sets]
                            if g is None or ks != list(range(1, len(ks) + 1)):
                                raise Unsupported("absolute assignments at irregular times on a reacting species")
                            vals = [self.init.get(tgt, 0.0)] + [v for _, v in sets]
                        self.ofiles.append({"filename": fn, "rows": [[num(-g)]] + [[num(v)] for v in vals]})
                        col7[tgt] = fn + "-SET"
                        if len(vals) * g <= t_end:
                            self.flag("partial", f"'{tgt}': Kintecus V2027 re-applies the last -SET value every "
                                                 f"{num(g)} time units after the profile ends (from t={num(len(vals) * g)}); "
                                                 f"the SBML sets it once")
                self.flag("info", f"time event(s) on '{tgt}' encoded as Constant File {fn}")
            except Unsupported as ex:
                self.flag("blocker", f"event assignment to '{tgt}' not converted ({ex})")
        if len(self.ofiles) > 10:
            self.flag("partial", f"{len(self.ofiles)} profile files exceed the GUI's 10 O_ file limit")
        return col7

    # driver ----------------------------------------------------------------------
    def build(self, t_end):
        self.t_end = float(t_end)
        self.classify()
        self.make_names()
        self.states()
        self.apply_offsets()
        self.reacting_ids = set()
        for r in self.m.getListOfReactions():
            try:
                reac, prod = self.stoich(r)
            except Unsupported:
                continue
            for sid in set(reac) | set(prod):
                if self.is_dyn(sid) and abs(prod.get(sid, 0) - reac.get(sid, 0)) > 0:
                    self.reacting_ids.add(sid)
        self.reacting_ids |= {v for v in self.rrule if v in self.spec}
        self.varying = set(self.reacting_ids) | {t for t in self.changed_by_event if t in self.spec}
        for r in self.m.getListOfReactions():
            try:
                self.convert_reaction(r)
            except Unsupported as ex:
                self.flag("blocker", f"reaction '{r.getId()}' not converted ({ex})")
        self.convert_rate_rules()
        if any("cTIME_" in rx for _, rx, _ in self.lines):
            self.lines.append(("1", f"{SOURCE}==>TIME_", "SBML time as a pseudo-species (d TIME_/dt = 1)"))
        self.reacting = set(self.reacting_ids)
        col7 = self.build_events(t_end)
        rows = []
        for sid in self.kspec:
            s = self.spec[sid]
            c0 = self.init.get(sid)
            if c0 is None:
                c0 = (s.getInitialConcentration() if s.isSetInitialConcentration() else
                      s.getInitialAmount() / self.vol(s.getCompartment()) if s.isSetInitialAmount() else 0.0)
            if c0 is None or not math.isfinite(float(c0)):
                self.flag("partial", f"species '{sid}' has no finite initial value; 0 used")
                c0 = 0.0
            pin = "No"
            if sid in col7:
                pin = col7[sid]
            elif sid not in self.rrule and (s.getBoundaryCondition() or s.getConstant() or sid not in self.reacting):
                pin = num(c0)
            if sid in self.offset:
                sh = self.offset[sid]
                rows.append((self.kname[sid], float(c0) + sh, "N", pin, f"{s.getName() or sid} + {num(sh)} (shifted state)",
                             f"__shifted__{sid}"))
                continue
            rows.append((self.kname[sid], c0, "Y", pin, s.getName() or sid, sid))
        for pid, name in self.pseudo.items():
            v0 = self.init.get(pid, self.par[pid].getValue())
            pin = col7.get(pid, "No" if pid in self.rrule else num(v0))
            if pid in self.offset:
                sh = self.offset[pid]
                rows.append((name, float(v0) + sh, "N", pin, f"parameter {pid} + {num(sh)} (shifted state)",
                             f"__shifted__{pid}"))
                continue
            rows.append((name, v0, "Y" if pid in self.rrule else "N", pin, f"parameter {pid}", pid))
        if self.time_used and any("TIME_" in rx for _, rx, _ in self.lines):
            rows.append(("TIME_", 0.0, "N", "No", "time", "__time__"))
        for nm, side_ in ((SOURCE, 0), (SINK, 1)):
            pat = re.compile(r"(^|\+)" + nm + r"(\+|$)", re.I)
            if any(pat.search(rx.split("==>")[side_]) for _, rx, _ in self.lines):
                rows.append((nm, 1.0, "N", "1", f"constant 1.0: {'nothing consumed' if side_ == 0 else 'nothing produced'}",
                             f"__{nm}__"))
        for (t0, up), nm in self.switches.items():
            fn = f"{nm}.TXT"
            if t0 <= 0 or t0 >= t_end:
                v = (1.0 if up else 0.0) if t0 <= 0 else (0.0 if up else 1.0)
                rows.append((nm, v, "N", num(v), f"time switch at {num(t0)}", "__switch__"))
                continue
            if self.use_pairs():
                self.ofiles.append({"filename": fn, "rows": self.pair_rows([(t0, 1.0 if up else 0.0)])})
                rows.append((nm, 0.0 if up else 1.0, "N", fn + "-SET",
                             f"time switch at {num(t0)} ({'off->on' if up else 'on->off'})", "__switch__"))
                self.switch_times = getattr(self, "switch_times", []) + [t0]
                continue
            dt = self.grid([t0, t_end])
            n = int(round(t_end / dt)) + 1
            if n > PROF_MAX:
                # the profile runs once and holds its last value (negative spacing), so it only has to reach just
                # past t0: spacing from t0 alone, two points after the switch
                dt = self.grid([t0])
                n = int(round(t0 / dt)) + 2
            if n > PROF_MAX:
                self.flag("blocker", f"time switch at {num(t0)} needs more than {PROF_MAX:,} profile points")
                continue
            vals = [(1.0 if ((k - 1) * dt >= t0 - 1e-12 * max(1.0, t_end)) == up else 0.0) for k in range(1, n + 1)]
            self.ofiles.append({"filename": fn, "rows": [[num(-dt)]] + [[num(v)] for v in vals]})
            rows.append((nm, 0.0 if up else 1.0, "N", fn, f"time switch at {num(t0)} ({'off->on' if up else 'on->off'})",
                         "__switch__"))
            self.switch_times = getattr(self, "switch_times", []) + [t0]
        if len(self.ofiles) > 10:
            self.flag("partial", f"{len(self.ofiles)} profile files exceed the GUI's 10 O_ file limit")
        npair = sum(1 for o in self.ofiles if self.is_pair_file(o))
        if PROFILE_FORMAT == "pairs" and npair < len(self.ofiles):
            self.flag("info", f"{npair} profiles as time/value pairs (Kintecus limit {PAIR_FILES_MAX} files), "
                              f"{len(self.ofiles) - npair} in the spacing format")
        if self.offset and not any(r[0] == SOURCE for r in rows):
            rows.append((SOURCE, 1.0, "N", "1", "constant 1.0: nothing consumed", f"__{SOURCE}__"))
        for sid, sh in self.offset.items():          # output-only species last (OUTC families sum top to bottom)
            lab = (self.spec[sid].getName() or sid) if sid in self.spec else f"parameter {sid}"
            rows.append((self.outname[sid], 0.0, "Y", "No", f"{lab} (= {self.kname[sid]} - {num(sh)})", sid,
                         f"OUTC[{self.kname[sid]};-{num(sh)}*{SOURCE}]"))   # factor first: '-source*2' (one-character
                                                                               # factor last) gives Kintecus soft error u82n
        self.species_rows = rows
        if not self.lines:
            self.flag("blocker", "no reactions or rate rules could be converted")
        return self


# ----------------------------------------------------------------------------- deck / kintpack writers
MODEL_COLS = ['Rate Const (A)', 'T^m', 'Ea', 'Reaction Equation', 'Comments']
SPECIES_COLS = ['Species', 'Residence (CSTR,s)', 'Initial Conc.', 'Display (Y/N)', 'External Conc.',
                'SSA? (Y/N)', 'Const File/#/No', 'Comments', '']
PARM_COLS = ['Starting Time', 'Max Time', 'Ea Units', 'Conc Units', 'X0', '']
THERM_COLS = ['DB Filename', 'Input Format', 'MAP Fields', 'Special Switches', 'Reservation List', '']


def coef(v):
    """a stoichiometric coefficient in plain decimal (no exponent: '1e-5A' would read as a species name)."""
    out = num(v)
    if "e" in out:
        out = f"{float(v):.20f}".rstrip("0").rstrip(".")
    return out


def num_up(v, digits=10):
    """v printed with `digits` significant digits, never below v. Kintecus applies a Constant File change at the
    first output step at or after the profile time, so k x (sampling interval) must not fall a hair short of an
    event time (1/7 printed exactly gives 7 x 0.14285714285714285 < 1 and the change lands one step late)."""
    v = float(v)
    if v == int(v):
        return num(v)
    out = num(v, digits)
    if float(out) < v:
        step = 10.0 ** (math.floor(math.log10(abs(v))) - digits + 1)
        out = num(math.ceil(v / step) * step, digits)
    return out


def grids(conv, header, t_end, n_points, acc=1e-8):
    model = [["#", "", "", "Model Description SpreadSheet", ""]]
    for h in header:
        model.append(["#", "", "", h[:250], ""])
    model.append(["#", "T^m", "Ea", "Reaction", "Comments"])
    for k, rx, cm in conv.lines:
        if cm:                                   # comment on its own line above the reaction (keeps lines <= 897)
            model.append(["#", "", "", cm[:MODEL_LINE_MAX - 8], ""])
        model.append([k, "0", "0", rx, ""])
    model.append(["END", "", "", "", ""])
    species = [["#", "", "", "", "", "", "", "Species Description Spreadsheet", ""],
               ["# Species", "Residence", "Initial", "Display", "External", "Species Special", "Constant File?",
                "Comments", ""]]
    for row in conv.species_rows:
        nm, c0, disp, pin, label, sid = row[:6]
        ssd = row[6] if len(row) > 6 else "No"
        species.append([nm, "0", num(c0), disp, "0", ssd, pin, f"{label} [{sid}]"[:120], ""])
    species.append(["END"] + [""] * 8)
    samp = t_end / n_points
    # with Constant File profiles, write the interval a hair long (relative 1e-9): Kintecus accumulates it, and
    # 100 x 0.1 = 9.99999999999998 < 10 would apply a change at t=10 one output step late
    samp_w = num_up(samp * (1 + 1e-9)) if conv.ofiles else num_up(samp)
    secs = t_end
    parm = [["#", "Parameter Description SpreadSheet", "", "", "", ""],
            ["# Starting Integration Time", "Maximum Integration Time", "Ea UNITS (Kcal, KJ,J,CAL,K)",
             "Conc. Units(moles/L, moles/cc, molecules/cc)", "X0", ""],
            [num(min(1e-6, samp / 1000)), samp_w, "KCAL", "MOLES/LITER", "0", ""],
            ["#Temperature (K) or Filename", "Pressure (Constant ? (Yes/No) )", "Volume Profile (Filename/No)",
             "External Heat Source", "X0", ""],
            ["298", "NOPE", "NO", "No", "0", ""],
            ["#       DAYS", "Hours", "Minutes", "Seconds", "PicoSeconds", ""],
            ["0", "0", "0", num(secs), "0", ""],
            ["#hv(filename)", "Sampling Interval (s)", "Percent(%)", "Accuracy", "X0", ""],
            ["None1", samp_w, "0", num(acc), "0", ""],
            ["# CSTR/PFR inlet Flow Temperature", "X0", "X0", "X0", "X0", ""],
            ["0", "0", "0", "0", "0", ""],
            ["END", "", "", "", "", ""]]
    return model, species, parm


def write_deck(d, model, species, parm, ofiles):
    """Write files exactly as the Kintecus GUI/VBA does: tab-joined, CRLF; species/parm get a trailing tab."""
    def w(fn, rows, n, trail):
        out = []
        for r in rows:
            cells = (list(r) + [""] * n)[:n]
            if cells[0] == "":
                cells[0] = "# "
            out.append("\t".join(cells) + ("\t" if trail else ""))
        open(os.path.join(d, fn), "w", newline="", encoding="latin-1", errors="replace").write("\r\n".join(out) + "\r\n")
    w("MODEL.DAT", model, 5, False)
    w("SPECIES.DAT", species, 9, True)
    w("PARM.DAT", parm, 6, True)
    for of in ofiles:
        open(os.path.join(d, of["filename"]), "w", newline="").write("\r\n".join("\t".join(r) for r in of["rows"]) + "\r\n")


def kintpack(model, species, parm, ofiles, switches, description, created):
    return {
        "version": 3,
        "created": created,
        "description": description,
        "settings": {"kpath": "", "switches": switches, "nation": False},
        "grids": {
            "model": {"cols": MODEL_COLS, "rows": model},
            "species": {"cols": SPECIES_COLS, "rows": species},
            "parm": {"cols": PARM_COLS, "rows": parm},
            "therm": {"cols": THERM_COLS, "rows": []},
            "fitdata": {"cols": [], "rows": []},
        },
        "ofiles": [{"filename": of["filename"], "key": "O_" + of["filename"], "cols": [f"Col{i + 1}" for i in range(max(len(r) for r in of["rows"]))], "rows": of["rows"]}
                   for of in ofiles],
    }


def convert(path, init, t_end, offsets=None):
    """Build the deck; if some USER expressions are too long for Kintecus, rebuild once with short aliases
    (Q1, Q2, ...) for the species/parameters named in those expressions and keep the better build."""
    conv = Converter(path, init, offsets=offsets).build(t_end)
    if not conv.long_lines and not conv.overlong:
        return conv
    inv = {("c" + n).upper(): sid for sid, n in conv.kname.items()}
    need = []
    for vars_ in [v for _, _, v in conv.long_lines] + [v for _, v in conv.overlong]:
        for v in vars_:
            sid = inv.get(v.upper())
            if sid and len(conv.kname[sid]) > 3 and sid not in need:
                need.append(sid)
    if not need:
        return conv
    taken = {n.upper() for n in conv.kname.values()}
    alias, i = {}, 1
    for sid in need:
        while f"Q{i}" in taken or f"Q{i}" in RESERVED or any(b in f"CQ{i}" for b in BAD_SUB):
            i += 1
        alias[sid] = f"Q{i}"; taken.add(f"Q{i}"); i += 1
    conv2 = Converter(path, init, alias, offsets=offsets).build(t_end)
    # fewer refused lines first (each is a reaction not converted), then less USER text over USER_SAFE
    excess = lambda c: (len(c.overlong), sum(max(0, n - USER_SAFE) for _, n, _ in c.long_lines))
    if excess(conv2) < excess(conv):
        conv2.flag("info", f"{len(alias)} species/parameters given short names (Q1...) so USER expressions fit "
                           f"Kintecus's expression length limit; original ids are in the SPECIES comments column")
        return conv2
    return conv
