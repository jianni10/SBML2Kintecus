#!/usr/bin/env python3
"""
kintpack2sbml.py -- convert a Kintecus model (.kintpack, or a folder / MODEL.DAT of a Kintecus deck) to SBML
Level 3 Version 2, a SED-ML time course, and/or a lossless Kintecus XML file. Anything SBML cannot represent is
exported as far as possible and listed - with the original text - in an XML comment block in the SBML file,
in the SBML model notes, and in a text report.

    kintpack2sbml.py INPUT [options]

INPUT   a .kintpack file, a folder containing MODEL.DAT / SPECIES.DAT / PARM.DAT, or the MODEL.DAT itself

Options
  -o DIR               output folder (default: next to the input)
  --name NAME          base name of the output files (default: the input's name)
  --format FMT         sbml (default) | xml (lossless Kintecus XML: every sheet and file) | both
  --no-sedml           do not write the SED-ML time course (<name>.sedml)
  --no-report          do not write <name>_kintecus_unsupported.txt
  --max-events N       largest number of SBML events generated from Constant File profiles [5000]
  -i                   international inputs: numbers with a decimal comma (4,34) are read as 4.34
  --validate           run the libSBML consistency check on the written SBML and print the result
  -h, --help           this help

Needs: pip install python-libsbml        (only for SBML output; --format xml needs nothing)
"""
import datetime, json, math, os, re, sys
from xml.sax.saxutils import escape, quoteattr

KIN_NS = "http://www.kintecus.com/sbml/kintecus"
R_GAS = {"KCAL": 1.98720425864083e-3, "KJ": 8.31446261815324e-3, "J": 8.31446261815324,
         "CAL": 1.98720425864083, "K": 1.0, "KELVIN": 1.0}
FUNMAP = {"exp": "exp", "log": "ln", "ln": "ln", "log10": "log10", "sqrt": "sqrt", "abs": "abs", "sin": "sin",
          "cos": "cos", "tan": "tan", "asin": "arcsin", "acos": "arccos", "atn": "arctan", "atan": "arctan",
          "sinh": "sinh", "cosh": "cosh", "tanh": "tanh", "floor": "floor", "ceil": "ceiling"}
NUM_COMMA = re.compile(r"^\s*[+-]?(\d*,\d+|\d+,\d*)([eEdD][+-]?\d+)?\s*$")
RESERVED = {"cell", "T", "R_Ea", "M_total", "time", "pi", "exponentiale", "true", "false", "avogadro", "inf",
            "nan", "notanumber", "infinity"}


class Unsupported(Exception):
    pass


# ------------------------------------------------------------------ reading input
def read_dat(path):
    with open(path, encoding="latin-1", errors="replace", newline="") as fh:
        txt = fh.read()
    return [l.split("\t") for l in re.split(r"\r\n|\n|\r", txt.rstrip("\r\n"))]


def load(inp):
    """-> dict(grids={name: rows}, files={filename: rows}, description, switches, source)"""
    if inp.lower().endswith(".kintpack"):
        p = json.load(open(inp, encoding="utf-8"))
        grids = {k.lower(): [[("" if c is None else str(c)) for c in r] for r in g.get("rows", [])]
                 for k, g in p.get("grids", {}).items() if g.get("rows")}
        cols = {k.lower(): g.get("cols", []) for k, g in p.get("grids", {}).items()}
        files = {o["filename"]: [[str(c) for c in r] for r in o.get("rows", [])] for o in p.get("ofiles", [])}
        st = p.get("settings") or {}
        return {"grids": grids, "cols": cols, "files": files, "description": p.get("description", ""),
                "switches": st.get("switches", ""), "nation": bool(st.get("nation")), "source": inp,
                "version": p.get("version"), "created": p.get("created", ""), "settings": st}
    d = inp if os.path.isdir(inp) else os.path.dirname(os.path.abspath(inp))
    names = {f.upper(): f for f in os.listdir(d)}
    grids, files = {}, {}
    for key, fn in (("model", "MODEL.DAT"), ("species", "SPECIES.DAT"), ("parm", "PARM.DAT"),
                    ("therm", "THERM.DAT"), ("transport", "TRANSPORT.DAT"), ("spname", "SPNAME.DAT")):
        if fn in names:
            grids[key] = read_dat(os.path.join(d, names[fn]))
    if "model" not in grids:
        raise SystemExit(f"no MODEL.DAT in {d}")
    for f in sorted(names):                       # profile files referenced by SPECIES column 7, fitdata files
        if f.startswith("FITDATA") and f.endswith(".TXT"):
            files[names[f]] = read_dat(os.path.join(d, names[f]))
    for r in grids.get("species", []):
        if len(r) > 6 and r[6].strip() and not _is_num(r[6]) and r[6].strip().upper() not in ("NO", "N", "YES"):
            fn = re.sub(r"-(SET|ADD)$", "", r[6].strip(), flags=re.I)
            if fn.upper() in names:
                files[names[fn.upper()]] = read_dat(os.path.join(d, names[fn.upper()]))
    for v in parm_fields(grids.get("parm", [])).values():      # temperature / volume / hv profile files
        if v and v.upper() in names and v.upper() not in ("MODEL.DAT", "SPECIES.DAT", "PARM.DAT"):
            files[names[v.upper()]] = read_dat(os.path.join(d, names[v.upper()]))
    return {"grids": grids, "cols": {}, "files": files, "description": "", "switches": "", "nation": False,
            "source": d, "version": None, "created": "", "settings": {}}


def _is_num(s):
    try:
        float(str(s).strip().replace("D", "E").replace("d", "e"))
        return True
    except ValueError:
        return False


class Reader:
    def __init__(self, intl):
        self.intl = intl

    def f(self, s, what="number"):
        s = str(s).strip()
        if self.intl and NUM_COMMA.match(s):
            s = s.replace(",", ".")
        s = s.replace("D", "E").replace("d", "e")
        try:
            return float(s)
        except ValueError:
            raise Unsupported(f"'{s}' is not a {what}")


def data_rows(rows):
    for i, r in enumerate(rows, 1):
        if not r or not r[0].strip() or r[0].lstrip().startswith("#"):
            continue
        if r[0].strip().upper() == "END":
            break
        yield i, r


def label_key(s):
    s = s.strip().lstrip("#").strip()
    s = re.split(r"\(", s)[0]
    s = re.split(r"\bor\b", s)[0]
    return re.sub(r"[^a-z0-9]", "", s.lower())


def parm_fields(rows):
    out, pos = {}, []
    for i, r in enumerate(rows):
        if r and r[0].strip().upper() == "END":
            break
        if r and r[0].lstrip().startswith("#") and i + 1 < len(rows):
            nxt = rows[i + 1]
            if nxt and nxt[0].strip() and not nxt[0].lstrip().startswith("#") and nxt[0].strip().upper() != "END":
                for c, lab in enumerate(r):
                    k = label_key(lab)
                    if k and k != "x0" and c < len(nxt):
                        out.setdefault(k, nxt[c].strip())
        elif r and r[0].strip() and not r[0].lstrip().startswith("#"):
            pos.append(r)
    # positional fallback (standard PARM layout)
    std = [["startingintegrationtime", "maximumintegrationtime", "eaunits", "concunits"],
           ["temperature", "pressure", "volumeprofile", "externalheatsource"],
           ["days", "hours", "minutes", "seconds", "picoseconds"],
           ["hv", "samplinginterval", "percent", "accuracy"]]
    for r, keys in zip(pos, std):
        for c, k in enumerate(keys):
            if k not in out and c < len(r):
                out[k] = r[c].strip()
    return out


# ------------------------------------------------------------------ reaction parsing
def split_brackets(rx):
    """remove X[...] tokens: -> (rest, [(kind, content)])"""
    toks, out, i = [], [], 0
    while i < len(rx):
        m = re.compile(r"([A-Za-z]+)\[").match(rx, i)
        if m and (i == 0 or rx[i - 1] in "+= >"):
            j = rx.find("]", m.end())
            if j < 0:
                raise Unsupported("unbalanced '[' in reaction")
            toks.append((m.group(1).upper(), rx[m.end():j]))
            i = j + 1
            continue
        out.append(rx[i]); i += 1
    rest = "".join(out)
    rest = re.sub(r"\+\s*(?=\+|==>|=|$)", "", rest)          # dangling '+'
    rest = re.sub(r"(^|==>|=)\s*\+", r"\1", rest)
    return rest, toks


def parse_side(s, species_names):
    """'2E2F' is 2 x E2F (Kintecus coefficients are plain decimals: no exponent), '0.5B', 'A+A'."""
    def known(x):
        return x in species_names or x.upper() in species_names
    terms = []
    for t in s.split("+"):
        t = t.strip()
        if not t or re.match(r"^(0\s+)?NULL$", t, re.I):
            continue
        if known(t):
            terms.append((t, 1.0)); continue
        plain = re.match(r"^(\d+\.?\d*|\.\d+)\s*(\S.*)$", t)
        expo = re.match(r"^(\d+\.?\d*(?:[eE][-+]?\d+)|\.\d+[eE][-+]?\d+)\s*(\S.*)$", t)
        if plain and known(plain.group(2).strip()):
            terms.append((plain.group(2).strip(), float(plain.group(1))))
        elif expo and known(expo.group(2).strip()):
            terms.append((expo.group(2).strip(), float(expo.group(1))))
        elif plain and not re.match(r"^\d", plain.group(2)):
            terms.append((plain.group(2).strip(), float(plain.group(1))))
        else:
            terms.append((t, 1.0))
    return terms


def parse_reaction(rx, species_names):
    rest, toks = split_brackets(rx)
    if "==>" in rest:
        lhs, rhs = rest.split("==>", 1); rev = False
    elif "<=>" in rest:
        lhs, rhs = rest.split("<=>", 1); rev = True
    elif "=" in rest:
        lhs, rhs = rest.split("=", 1); rev = True
    else:
        raise Unsupported("no '==>' or '=' in the reaction")
    return parse_side(lhs, species_names), parse_side(rhs, species_names), rev, toks


# ------------------------------------------------------------------ the converter
class Kin2SBML:
    def __init__(self, data, intl=False, max_events=5000):
        self.d = data
        self.rd = Reader(intl or data.get("nation"))
        self.max_events = max_events
        self.unsupported = []        # (where, what, original text)
        self.notes = []              # informational
        self.sid = {}                # kintecus name (upper) -> SBML id
        self.used = set(RESERVED)

    def unsup(self, where, what, text=""):
        self.unsupported.append((where, what, text))

    def make_id(self, name, prefix=""):
        base = re.sub(r"[^A-Za-z0-9_]", "_", prefix + name)
        if not re.match(r"[A-Za-z_]", base):
            base = "s_" + base
        sid, i = base, 1
        while sid in self.used:
            i += 1; sid = f"{base}_{i}"
        self.used.add(sid)
        return sid

    # -------------------------------------------------------------- build
    def build(self):
        import libsbml
        self.L = libsbml
        g = self.d["grids"]
        self.doc = libsbml.SBMLDocument(3, 2)
        m = self.m = self.doc.createModel()
        title = self._title()
        m.setId(self.make_id(re.sub(r"\W+", "_", title)[:40] or "kintecus_model"))
        m.setName(title)
        # ---- PARM
        P = parm_fields(g.get("parm", []))
        self.P = P
        units = (P.get("concunits") or "MOLES/LITER").upper().replace(" ", "")
        ea_units = (P.get("eaunits") or "KCAL").upper()
        self._units(units)
        comp = m.createCompartment(); comp.setId("cell"); comp.setConstant(True); comp.setSize(1.0)
        comp.setSpatialDimensions(3)
        temp = P.get("temperature", "298")
        if _is_num(temp):
            T = self.rd.f(temp)
        else:
            T = 298.0
            self.unsup("PARM.DAT", f"temperature profile file '{temp}' - T held at 298 K in the SBML",
                       self._file_text(temp))
            self.used_files = getattr(self, "used_files", set()) | {temp.strip().upper()}
        self._param("T", T, "Temperature (K) from PARM.DAT")
        if ea_units not in R_GAS:
            self.unsup("PARM.DAT", f"Ea units '{ea_units}' not recognized; KCAL assumed")
            ea_units = "KCAL"
        self.R = R_GAS[ea_units]; self.ea_units = ea_units
        for key, label in (("pressure", "constant pressure"), ("volumeprofile", "volume profile"),
                           ("externalheatsource", "external heat source"), ("hv", "photolysis (hv) file")):
            v = (P.get(key) or "").strip()
            if v and v.upper() not in ("NO", "N", "NOPE", "NONE", "NONE1", "0", ""):
                self.unsup("PARM.DAT", f"{label}: '{v}' (not representable in SBML)", self._file_text(v))
                self.used_files = getattr(self, "used_files", set()) | {v.upper()}
        try:
            days, hours, mins, secs, ps = (self.rd.f(P.get(k) or 0) for k in
                                            ("days", "hours", "minutes", "seconds", "picoseconds"))
            self.t_end = days * 86400 + hours * 3600 + mins * 60 + secs + ps * 1e-12
        except Unsupported:
            self.t_end = 0.0
        try:
            self.samp = self.rd.f(P.get("samplinginterval") or 0)
        except Unsupported:
            self.samp = 0.0
        sw = self.d.get("switches") or ""
        if re.search(r"-THERM|-EQUIL", sw, re.I):
            self.unsup("switches", f"thermodynamic/equilibrium mode ({sw}) is not exported; the SBML is the "
                                   f"isothermal kinetic model")
        # ---- species
        self.species = {}            # upper name -> dict
        rows = g.get("species", [])
        for line, r in data_rows(rows):
            r = (r + [""] * 9)[:9]
            nm = r[0].strip()
            try:
                c0 = self.rd.f(r[2] or 0, "number (initial concentration)")
            except Unsupported as ex:
                self.unsup(f"SPECIES.DAT line {line}", str(ex), "\t".join(r)); c0 = 0.0
            self.species[nm.upper()] = {"name": nm, "init": c0, "res": r[1].strip(), "disp": r[3].strip(),
                                        "ext": r[4].strip(), "ssa": r[5].strip(), "pin": r[6].strip(),
                                        "comment": r[7].strip(), "line": line, "raw": "\t".join(r).rstrip("\t")}
        # species that appear only in MODEL.DAT
        self.lines = list(data_rows(g.get("model", [])))
        mrows = g.get("model", [])
        for j, (line, r) in enumerate(self.lines):      # a '#' comment line right above a reaction is its comment
            prev = mrows[line - 2] if line >= 2 else None
            if prev and prev[0].strip() == "#" and not "".join(prev[1:3]).strip() and len(prev) > 3 \
                    and prev[3].strip() and not (len(r) > 4 and r[4].strip()):
                self.lines[j] = (line, (list(r) + [""] * 5)[:4] + [prev[3].strip()])
        names = set(self.species)
        for line, r in self.lines:
            rx = r[3] if len(r) > 3 else (r[1] if len(r) > 1 else "")
            try:
                lhs, rhs, _, _ = parse_reaction(rx, names)
            except Unsupported:
                continue
            for s, _ in lhs + rhs:
                if s.upper() not in self.species and s.upper() not in ("M", "HV"):
                    self.species[s.upper()] = {"name": s, "init": 0.0, "res": "", "disp": "", "ext": "", "ssa": "",
                                               "pin": "", "comment": "", "line": None, "raw": ""}
                    self.notes.append(f"species {s} is used in MODEL.DAT but not listed in SPECIES.DAT; "
                                      f"initial concentration 0")
        for key, sp in self.species.items():
            self._species(sp)
        self.var = {("c" + sp["name"]).lower(): self.sid[k] for k, sp in self.species.items()}
        # ---- reactions
        self.need_M = False
        for n, (line, r) in enumerate(self.lines, 1):
            try:
                self._reaction(n, line, r)
            except Unsupported as ex:
                self.unsup(f"MODEL.DAT line {line}", f"reaction not exported: {ex}", "\t".join(r).rstrip("\t"))
        if self.need_M:
            p = self._param("M_total", 0.0, "Third-body concentration M = sum of all species concentrations",
                            constant=False)
            ar = m.createAssignmentRule(); ar.setVariable("M_total")
            terms = [self.sid[k] for k, sp in self.species.items()] or ["0"]
            ar.setMath(libsbml.parseL3Formula(" + ".join(terms)))
            self.notes.append("third body M exported as M_total = sum of all species concentrations "
                              "(collision efficiencies M[...] are not exported)")
        # ---- CSTR flow
        for k, sp in self.species.items():
            if sp["res"] and _is_num(sp["res"]) and self.rd.f(sp["res"]) > 0:
                self._cstr(sp)
        # ---- events from Constant File profiles
        self._events()
        # ---- other sheets and files
        for key in g:
            if key not in ("model", "species", "parm"):
                self.unsup(f"{key.upper()} sheet", "not representable in SBML; full content below",
                           "\n".join("\t".join(r).rstrip("\t") for r in g[key]))
        used_files = getattr(self, "used_files", set())
        for fn, rows in self.d["files"].items():
            if fn.upper() not in used_files:
                self.unsup(f"file {fn}", "not used by the SBML model (fit data or other data file); full content below",
                           "\n".join("\t".join(r).rstrip("\t") for r in rows))
        self._notes()
        return self.doc

    def _title(self):
        desc = self.d.get("description") or ""
        if desc.strip():
            return desc.strip().split("\n")[0][:200]
        for r in self.d["grids"].get("model", [])[:4]:
            if r and r[0].lstrip().startswith("#") and len(r) > 3 and r[3].strip() and \
                    "Description SpreadSheet" not in r[3]:
                return r[3].strip()[:200]
        return os.path.splitext(os.path.basename(self.d["source"]))[0]

    def _file_text(self, name):
        for fn, rows in self.d["files"].items():
            if fn.upper() == str(name).strip().upper():
                return "\n".join("\t".join(r) for r in rows)
        return ""

    def _units(self, units):
        L, m = self.L, self.m
        def ud(uid, parts):
            u = m.createUnitDefinition(); u.setId(uid)
            for kind, exp, scale, mult in parts:
                x = u.createUnit(); x.setKind(kind); x.setExponent(exp); x.setScale(scale); x.setMultiplier(mult)
        if units in ("MOLES/LITER", "MOLES/L", "M", "MOL/L", "MOLES/LITRE"):
            sub, vol = L.UNIT_KIND_MOLE, ("litre", [(L.UNIT_KIND_LITRE, 1, 0, 1.0)])
        elif units in ("MOLES/CC", "MOL/CC", "MOLES/CM3"):
            sub, vol = L.UNIT_KIND_MOLE, ("cc", [(L.UNIT_KIND_LITRE, 1, -3, 1.0)])
        elif units in ("MOLECULES/CC", "MOLECULES/CM3"):
            sub, vol = L.UNIT_KIND_ITEM, ("cc", [(L.UNIT_KIND_LITRE, 1, -3, 1.0)])
        else:
            sub, vol = L.UNIT_KIND_MOLE, ("litre", [(L.UNIT_KIND_LITRE, 1, 0, 1.0)])
            self.unsup("PARM.DAT", f"concentration units '{units}' not recognized; mole/litre assumed")
        ud("substance_unit", [(sub, 1, 0, 1.0)])
        ud("volume_unit", vol[1])
        m.setSubstanceUnits("substance_unit"); m.setExtentUnits("substance_unit")
        m.setVolumeUnits("volume_unit"); m.setTimeUnits("second")
        self.units = units

    def _param(self, pid, val, note="", constant=True):
        p = self.m.createParameter(); p.setId(pid); p.setValue(val); p.setConstant(constant)
        self.used.add(pid)
        if note:
            p.setNotes(f'<body xmlns="http://www.w3.org/1999/xhtml"><p>{escape(note)}</p></body>')
        return p

    def _species(self, sp):
        L = self.L
        sid = self.make_id(sp["name"])
        self.sid[sp["name"].upper()] = sid
        s = self.m.createSpecies(); s.setId(sid); s.setName(sp["name"]); s.setCompartment("cell")
        s.setHasOnlySubstanceUnits(False); s.setInitialConcentration(sp["init"])
        s.setBoundaryCondition(False); s.setConstant(False)
        pin = sp["pin"]
        sp["mode"] = None
        if pin and pin.upper() not in ("NO", "N", ""):
            if _is_num(pin) or (self.rd.intl and NUM_COMMA.match(pin)):
                v = self.rd.f(pin)
                s.setBoundaryCondition(True); s.setConstant(True); s.setInitialConcentration(v)
                if sp["line"] and abs(v - sp["init"]) > 1e-12 * max(1.0, abs(v)):
                    self.notes.append(f"species {sp['name']} is held at {v:g} (SPECIES column 7) - its "
                                      f"initial concentration {sp['init']:g} is not used")
            else:
                mode = "SET" if pin.upper().endswith("-SET") else "ADD" if pin.upper().endswith("-ADD") else "PIN"
                fn = re.sub(r"-(SET|ADD)$", "", pin, flags=re.I)
                rows = None
                for f, rr in self.d["files"].items():
                    if f.upper() == fn.upper():
                        rows = rr; self.used_files = getattr(self, "used_files", set()) | {f.upper()}
                if rows is None:
                    self.unsup(f"SPECIES.DAT line {sp['line']}", f"Constant File '{fn}' for species {sp['name']} "
                                                                f"not found in the pack; species exported as free")
                elif rows and rows[0] and "TIME" in str(rows[0][0]).upper():
                    # Time/value pairs (Kintecus build 2026-11-03): header line, "t<TAB>value" lines, END
                    pairs = []
                    try:
                        for r in rows[1:]:
                            if not r or not r[0].strip():
                                continue
                            if r[0].strip().upper() == "END":
                                break
                            pairs.append((self.rd.f(r[0]), self.rd.f(r[1])))
                    except (Unsupported, IndexError) as ex:
                        self.unsup(f"file {fn}", f"time/value pairs not numeric ({ex}); not exported", self._file_text(fn))
                        pairs = []
                    if pairs:
                        sp["mode"], sp["pairs"], sp["profile_file"] = mode, pairs, fn
                        if mode == "PIN":
                            s.setBoundaryCondition(True)
                        for t, v in pairs:
                            if t <= 0:            # applied at the start
                                s.setInitialConcentration(v if mode != "ADD" else sp["init"] + v)
                else:
                    try:
                        vals = [self.rd.f(r[0]) for r in rows if r and r[0].strip()]
                    except Unsupported as ex:
                        self.unsup(f"file {fn}", f"profile not numeric ({ex}); not exported", self._file_text(fn))
                        vals = []
                    if len(vals) >= 2:
                        sp["mode"], sp["profile"], sp["profile_file"] = mode, vals, fn
                        if mode == "PIN":
                            s.setBoundaryCondition(True); s.setInitialConcentration(vals[1])
                        elif mode == "SET":
                            s.setInitialConcentration(vals[1])
                        else:
                            s.setInitialConcentration(sp["init"] + vals[1])
        n = [f"Kintecus SPECIES.DAT: {sp['raw']}" if sp["raw"] else "Kintecus: not in SPECIES.DAT"]
        if sp["comment"]:
            n.append(f"Comment: {sp['comment']}")
        if sp["ssa"] and sp["ssa"].upper() not in ("NO", "N"):
            self.unsup(f"SPECIES.DAT line {sp['line']}", f"species special/SSA setting '{sp['ssa']}' for "
                                                        f"{sp['name']} not exported", sp["raw"])
        s.setNotes('<body xmlns="http://www.w3.org/1999/xhtml">' +
                   "".join(f"<p>{escape(x)}</p>" for x in n) + "</body>")

    def _expr(self, expr, where):
        expr = re.sub(r"(?<![A-Za-z0-9_.])(\d+\.?\d*|\.\d+)[dD]([-+]?\d+)", r"\1e\2", expr)
        if self.rd.intl:
            pass                                     # USER expressions always use '.' decimals and ',' lists

        def rep(mo):
            name, paren = mo.group(1), mo.group(2)
            if paren:
                f = FUNMAP.get(name.lower())
                if not f:
                    raise Unsupported(f"function '{name}' in a USER expression has no SBML equivalent")
                return f + paren
            low = name.lower()
            if low in self.var:
                return self.var[low]
            if low == "t":
                return "T"
            if low == "m":
                self.need_M = True; return "M_total"
            if low == "pi":
                return "pi"
            raise Unsupported(f"unknown variable '{name}' in a USER expression")
        return re.sub(r"(?<![0-9.])\b([A-Za-z_][A-Za-z0-9_]*)\b(\s*\()?", rep, expr)

    def _reaction(self, n, line, r):
        L = self.L
        r = (r + [""] * 5)[:5]
        A, mexp, Ea, rx, cm = r[0].strip(), r[1].strip() or "0", r[2].strip() or "0", r[3].strip(), r[4].strip()
        lhs, rhs, rev, toks = parse_reaction(rx, set(self.species))
        user = None
        for kind, content in toks:
            if kind == "S" and content.upper().startswith("USER;"):
                body = content[5:]
                if ":" not in body:
                    raise Unsupported("USER term without ':'")
                vars_, ex = body.split(":", 1)
                user = ex.strip()
            elif kind == "M":
                raise Unsupported(f"third-body collision efficiencies M[{content}]")
            else:
                raise Unsupported(f"special rate form {kind}[{content[:60]}] (fall-off/Troe/SRI/PLOG/Chebyshev "
                                  f"or other Kintecus special reaction)")
        third = False
        for side in (lhs, rhs):
            for i, (s, c) in enumerate(list(side)):
                if s.upper() == "HV" and s.upper() not in self.species:
                    raise Unsupported("photolysis (hv) reaction")
        if any(s.upper() == "M" and "M" not in self.species for s, _ in lhs):
            third = True; self.need_M = True
        lhs = [(s, c) for s, c in lhs if not (s.upper() == "M" and "M" not in self.species)]
        rhs = [(s, c) for s, c in rhs if not (s.upper() == "M" and "M" not in self.species)]
        rid = self.make_id(f"R{line}")
        rxn = self.m.createReaction(); rxn.setId(rid); rxn.setReversible(bool(rev))
        if cm:
            rxn.setName(cm[:200])
        net = {}
        for s, c in lhs:
            net[s.upper()] = net.get(s.upper(), 0.0) - c
        for s, c in rhs:
            net[s.upper()] = net.get(s.upper(), 0.0) + c
        appear = {s.upper() for s, _ in lhs + rhs}
        for k in appear:
            v = net[k]
            if abs(v) < 1e-15:
                ms = rxn.createModifier(); ms.setSpecies(self.sid[k])
            elif v < 0:
                sr = rxn.createReactant(); sr.setSpecies(self.sid[k]); sr.setStoichiometry(-v); sr.setConstant(True)
            else:
                sr = rxn.createProduct(); sr.setSpecies(self.sid[k]); sr.setStoichiometry(v); sr.setConstant(True)
        # rate = (USER or k) * prod[reactants]  (Kintecus semantics), times the compartment volume
        conc = []
        for s, c in lhs:
            sid = self.sid[s.upper()]
            conc.append(sid if abs(c - 1) < 1e-15 else f"{sid}^{c:.17g}")
        if third:
            conc.append("M_total")
        kl = rxn.createKineticLaw()
        if user is not None:
            core = f"({self._expr(user, f'MODEL.DAT line {line}')})"
        else:
            a = self.rd.f(A, "rate constant"); mm = self.rd.f(mexp, "T^m exponent"); ea = self.rd.f(Ea, "Ea")
            ids = {}
            for pid, v in (("A", a), ("m", mm), ("Ea", ea)):
                lid = "kin_" + pid
                while lid in self.used:
                    lid += "_"
                ids[pid] = lid
                lp = kl.createLocalParameter(); lp.setId(lid); lp.setValue(v)
            if mm == 0 and ea == 0:
                core = ids["A"]
            else:
                core = ids["A"] + ("" if mm == 0 else f" * T^{ids['m']}") + \
                    ("" if ea == 0 else f" * exp(-{ids['Ea']} / (R_Ea * T))")
                if not self.m.getParameter("R_Ea"):
                    self._param("R_Ea", self.R, f"gas constant in Ea units ({self.ea_units}/mol/K)")
        f = " * ".join(["cell", core] + conc)
        ast = L.parseL3FormulaWithSettings(f, self._l3)
        if ast is None:
            raise Unsupported(f"rate expression could not be parsed ({L.getLastParseL3Error()})")
        kl.setMath(ast)
        declared = {self.sid[k] for k in appear}
        species_ids = set(self.sid.values())
        for nm in _ast_names(ast):
            if nm in species_ids and nm not in declared:
                ms = rxn.createModifier(); ms.setSpecies(nm); declared.add(nm)
        ann = (f'<kintecus:line xmlns:kintecus="{KIN_NS}" file="MODEL.DAT" number="{line}">'
               f'{escape(chr(9).join(r).rstrip(chr(9)))}</kintecus:line>')
        rxn.setAnnotation(ann)
        if rev:
            self.unsup(f"MODEL.DAT line {line}", "reversible reaction: Kintecus computes the reverse rate from "
                                                 "thermodynamic data; only the forward rate law is exported",
                       "\t".join(r).rstrip("\t"))
        if cm:
            rxn.setNotes(f'<body xmlns="http://www.w3.org/1999/xhtml"><p>{escape(cm)}</p></body>')

    @property
    def _l3(self):
        if not hasattr(self, "_l3s"):
            s = self.L.L3ParserSettings()
            s.setParseLog(self.L.L3P_PARSE_LOG_AS_LN)
            s.setModel(self.m)
            self._l3s = s
        return self._l3s

    def _cstr(self, sp):
        tau = self.rd.f(sp["res"]); ext = self.rd.f(sp["ext"] or 0)
        sid = self.sid[sp["name"].upper()]
        if self.m.getSpecies(sid).getBoundaryCondition():
            return
        for direction in ("in", "out"):
            if direction == "in" and ext == 0:
                continue
            r = self.m.createReaction(); r.setId(self.make_id(f"flow_{direction}_{sid}")); r.setReversible(False)
            kl = r.createKineticLaw()
            lp = kl.createLocalParameter(); lp.setId("tau"); lp.setValue(tau)
            if direction == "in":
                lp2 = kl.createLocalParameter(); lp2.setId("ext"); lp2.setValue(ext)
                pr = r.createProduct(); pr.setSpecies(sid); pr.setStoichiometry(1); pr.setConstant(True)
                kl.setMath(self.L.parseL3Formula("cell * ext / tau"))
            else:
                re_ = r.createReactant(); re_.setSpecies(sid); re_.setStoichiometry(1); re_.setConstant(True)
                kl.setMath(self.L.parseL3Formula(f"cell * {sid} / tau"))
        self.notes.append(f"species {sp['name']}: CSTR residence time {tau:g} s, external concentration {ext:g} "
                          f"-> flow reactions in = ext/tau, out = [{sp['name']}]/tau")

    def _events(self):
        L = self.L
        sched = {}                                     # time -> list of (sid, formula)
        t_end = self.t_end if self.t_end > 0 else None
        for k, sp in self.species.items():
            if not sp.get("mode"):
                continue
            if sp.get("pairs"):
                sid = self.sid[k]
                for t, val in sp["pairs"]:
                    if t <= 0 or (t_end is not None and t > t_end + 1e-12 * max(1, t_end)):
                        continue
                    if sp["mode"] == "ADD":
                        if val != 0:
                            sched.setdefault(round(t, 12), []).append((sid, f"{sid} + {val!r}"))
                    else:
                        sched.setdefault(round(t, 12), []).append((sid, repr(val)))
                self.notes.append(f"species {sp['name']}: Constant File {sp['profile_file']} ({sp['mode']}, "
                                  f"{len(sp['pairs'])} time/value pairs) exported as SBML events")
                continue
            vals = sp["profile"]; s = vals[0]; v = vals[1:]; sid = self.sid[k]
            step = abs(s)
            if step <= 0:
                self.unsup(f"file {sp['profile_file']}", "profile spacing 0; not exported", self._file_text(sp["profile_file"]))
                continue
            repeat = s > 0
            if repeat and t_end is None:
                self.notes.append(f"profile {sp['profile_file']} repeats (positive spacing) but no run length is "
                                  f"known; exported once")
            n_steps = len(v)
            total = (int(math.ceil(t_end / step)) + 1) if (repeat and t_end) else n_steps
            prev = v[0] if sp["mode"] == "PIN" else None
            for kk in range(1, total):
                t = kk * step
                if t_end is not None and t > t_end + 1e-12 * max(1, t_end):
                    break
                val = v[kk % n_steps] if repeat else (v[kk] if kk < n_steps else None)
                if val is None:
                    break
                if sp["mode"] == "PIN":
                    if val == prev:
                        continue
                    prev = val
                    sched.setdefault(round(t, 12), []).append((sid, repr(val)))
                elif sp["mode"] == "SET":
                    sched.setdefault(round(t, 12), []).append((sid, repr(val)))
                elif val != 0:
                    sched.setdefault(round(t, 12), []).append((sid, f"{sid} + {val!r}"))
            if repeat:
                self.notes.append(f"profile {sp['profile_file']} (positive spacing) repeats; events generated up to "
                                  f"the run length {t_end:g}")
            self.notes.append(f"species {sp['name']}: Constant File {sp['profile_file']} ({sp['mode']}, spacing "
                              f"{step:g}) exported as SBML events")
        times = sorted(sched)
        if len(times) > self.max_events:
            self.unsup("Constant File profiles", f"{len(times)} event times; only the first {self.max_events} are "
                                                 f"exported (raise --max-events)")
            times = times[:self.max_events]
        for j, t in enumerate(times, 1):
            e = self.m.createEvent(); e.setId(self.make_id(f"profile_t{j}")); e.setUseValuesFromTriggerTime(True)
            tr = e.createTrigger(); tr.setPersistent(True); tr.setInitialValue(False)
            tr.setMath(L.parseL3Formula(f"time >= {t!r}"))
            for sid, fml in sched[t]:
                ea = e.createEventAssignment(); ea.setVariable(sid); ea.setMath(L.parseL3Formula(fml))

    def _notes(self):
        d = self.d
        p = [f"Converted from Kintecus by kintpack2sbml.py on {datetime.date.today()} from "
             f"{os.path.basename(d['source'])}.",
             f"Run length {self.t_end:g} s, sampling interval {self.samp:g} s (see the SED-ML file)."
             if self.t_end else "Run length not given in PARM.DAT.",
             f"Kintecus switches stored in the pack: {d.get('switches') or '(none)'}",
             "Rate laws follow Kintecus semantics: rate = (k or USER expression) x product of reactant "
             "concentrations; k = A * T^m * exp(-Ea/(R T)). Compartment 'cell' (volume 1) carries the concentrations."]
        if d.get("description"):
            p.append("Pack description:")
            p += [l for l in d["description"].split("\n")[:60]]
        if self.notes:
            p.append("Conversion notes:")
            p += [f"- {x}" for x in self.notes]
        if self.unsupported:
            p.append(f"NOT REPRESENTED IN SBML ({len(self.unsupported)} items; the XML comment at the top of this file "
                     f"and the _kintecus_unsupported.txt report have the full original text):")
            p += [f"- {w}: {x}" for w, x, _ in self.unsupported]
        self.m.setNotes('<body xmlns="http://www.w3.org/1999/xhtml">' + "".join(f"<p>{escape(x)}</p>" for x in p) +
                        "</body>")


def _ast_names(ast):
    out, stack = set(), [ast]
    while stack:
        n = stack.pop()
        if n is None:
            continue
        if n.isName():
            out.add(n.getName())
        stack.extend(n.getChild(i) for i in range(n.getNumChildren()))
    return out


# ------------------------------------------------------------------ writers
def comment_block(conv):
    if not conv.unsupported:
        return ""
    out = ["", "KINTECUS ITEMS NOT REPRESENTED IN THIS SBML MODEL (exported by kintpack2sbml.py)",
           "=" * 80]
    for w, x, t in conv.unsupported:
        out.append(f"[{w}] {x}")
        if t:
            out += ["    " + l for l in t.split("\n")]
    txt = "\n".join(out) + "\n"
    txt = txt.replace("--", "- -").replace("--", "- -")
    return f"<!--{txt}-->\n"


def write_sbml(conv, path):
    import libsbml
    xml = libsbml.writeSBMLToString(conv.doc)
    blk = comment_block(conv)
    if blk:
        i = xml.find("?>") + 2
        xml = xml[:i] + "\n" + blk + xml[i:].lstrip("\n")
    open(path, "w", encoding="utf-8").write(xml)


def write_report(conv, path):
    L = [f"kintpack2sbml report - {os.path.basename(conv.d['source'])}   ({datetime.date.today()})", ""]
    m = conv.m
    L.append(f"SBML: {m.getNumSpecies()} species, {m.getNumReactions()} reactions, {m.getNumEvents()} events, "
             f"{m.getNumParameters()} parameters, {m.getNumRules()} rules")
    if conv.notes:
        L += ["", "Conversion notes:"] + [f"  - {x}" for x in conv.notes]
    L += ["", f"NOT REPRESENTED IN SBML: {len(conv.unsupported)} item(s)"]
    for w, x, t in conv.unsupported:
        L += ["", f"[{w}] {x}"]
        if t:
            L += ["    " + l for l in t.split("\n")]
    open(path, "w", encoding="utf-8", newline="\r\n").write("\n".join(L) + "\n")


def write_sedml(conv, path, sbml_name):
    if not conv.t_end or conv.t_end <= 0:
        return False
    n = int(round(conv.t_end / conv.samp)) if conv.samp > 0 else 1000
    n = max(1, min(n, 10 ** 6))
    disp = [(k, sp) for k, sp in conv.species.items() if sp["disp"].upper().startswith("Y")] or list(conv.species.items())
    dg, ds = [], []
    for k, sp in disp:
        sid = conv.sid[k]
        dg.append(f'''    <dataGenerator id="dg_{sid}" name={quoteattr(sp["name"])}>
      <listOfVariables><variable id="v_{sid}" taskReference="task1" target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='{sid}']"/></listOfVariables>
      <math xmlns="http://www.w3.org/1998/Math/MathML"><ci>v_{sid}</ci></math>
    </dataGenerator>''')
        ds.append(f'      <dataSet id="ds_{sid}" label={quoteattr(sp["name"])} dataReference="dg_{sid}"/>')
    xml = f'''<?xml version="1.0" encoding="UTF-8"?>
<sedML xmlns="http://sed-ml.org/sed-ml/level1/version3" xmlns:sbml="http://www.sbml.org/sbml/level3/version2/core" level="1" version="3">
  <listOfModels>
    <model id="model1" language="urn:sedml:language:sbml.level-3.version-2" source={quoteattr(sbml_name)}/>
  </listOfModels>
  <listOfSimulations>
    <uniformTimeCourse id="sim1" initialTime="0" outputStartTime="0" outputEndTime="{conv.t_end!r}" numberOfSteps="{n}">
      <algorithm kisaoID="KISAO:0000019"/>
    </uniformTimeCourse>
  </listOfSimulations>
  <listOfTasks>
    <task id="task1" modelReference="model1" simulationReference="sim1"/>
  </listOfTasks>
  <listOfDataGenerators>
    <dataGenerator id="dg_time" name="time">
      <listOfVariables><variable id="v_time" taskReference="task1" symbol="urn:sedml:symbol:time"/></listOfVariables>
      <math xmlns="http://www.w3.org/1998/Math/MathML"><ci>v_time</ci></math>
    </dataGenerator>
{chr(10).join(dg)}
  </listOfDataGenerators>
  <listOfOutputs>
    <report id="report1" name="time course">
      <listOfDataSets>
      <dataSet id="ds_time" label="Time(s)" dataReference="dg_time"/>
{chr(10).join(ds)}
      </listOfDataSets>
    </report>
  </listOfOutputs>
</sedML>
'''
    open(path, "w", encoding="utf-8").write(xml)
    return True


def write_kintecus_xml(data, path):
    """Lossless dump: every sheet (with column headers) and every file, cell by cell."""
    o = ['<?xml version="1.0" encoding="UTF-8"?>',
         f'<kintecus xmlns="{KIN_NS}" version={quoteattr(str(data.get("version") or ""))} '
         f'created={quoteattr(data.get("created") or "")} source={quoteattr(os.path.basename(data["source"]))}>']
    o.append(f"  <description>{escape(data.get('description') or '')}</description>")
    st = data.get("settings") or {}
    o.append("  <settings " + " ".join(f"{k}={quoteattr(str(v))}" for k, v in st.items()) + "/>")
    fnames = {"model": "MODEL.DAT", "species": "SPECIES.DAT", "parm": "PARM.DAT", "therm": "THERM.DAT",
              "transport": "TRANSPORT.DAT", "fitdata": "FITDATA.TXT"}
    for key, rows in data["grids"].items():
        o.append(f'  <sheet name={quoteattr(key)} file={quoteattr(fnames.get(key, key.upper() + ".DAT"))}>')
        cols = (data.get("cols") or {}).get(key) or []
        if cols:
            o.append("    <columns>" + "".join(f"<col>{escape(str(c))}</col>" for c in cols) + "</columns>")
        for i, r in enumerate(rows, 1):
            o.append(f'    <row n="{i}">' + "".join(f"<c>{escape(c)}</c>" for c in r) + "</row>")
        o.append("  </sheet>")
    for fn, rows in data["files"].items():
        o.append(f"  <file name={quoteattr(fn)}>")
        for i, r in enumerate(rows, 1):
            o.append(f'    <row n="{i}">' + "".join(f"<c>{escape(c)}</c>" for c in r) + "</row>")
        o.append("  </file>")
    o.append("</kintecus>")
    open(path, "w", encoding="utf-8").write("\n".join(o) + "\n")


# ------------------------------------------------------------------ main
def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    a = {"in": None, "o": None, "name": None, "fmt": "sbml", "sedml": True, "report": True, "max": 5000,
         "i": False, "validate": False}
    i = 0
    while i < len(argv):
        x = argv[i]
        nxt = lambda: argv[i + 1] if i + 1 < len(argv) else sys.exit(f"{x} needs a value")
        if x in ("-h", "--help", "/?"):
            print(__doc__); return 0
        elif x == "-o": a["o"] = nxt(); i += 1
        elif x == "--name": a["name"] = nxt(); i += 1
        elif x == "--format": a["fmt"] = nxt().lower(); i += 1
        elif x == "--no-sedml": a["sedml"] = False
        elif x == "--no-report": a["report"] = False
        elif x == "--max-events": a["max"] = int(nxt()); i += 1
        elif x == "-i": a["i"] = True
        elif x == "--validate": a["validate"] = True
        elif x.startswith("-"): print(f"kintpack2sbml: unknown option {x} (see --help)", file=sys.stderr); return 2
        elif a["in"] is None: a["in"] = x
        else: print(f"kintpack2sbml: unexpected argument {x}", file=sys.stderr); return 2
        i += 1
    if not a["in"]:
        print(__doc__); return 2
    if a["fmt"] not in ("sbml", "xml", "both"):
        print("kintpack2sbml: --format must be sbml, xml or both", file=sys.stderr); return 2
    inp = os.path.abspath(a["in"])
    if not os.path.exists(inp):
        print(f"kintpack2sbml: {inp} not found", file=sys.stderr); return 2
    data = load(inp)
    base = a["name"] or (os.path.splitext(os.path.basename(inp))[0] if not os.path.isdir(inp) else os.path.basename(inp))
    if base.upper() == "MODEL":
        base = os.path.basename(os.path.dirname(inp)) or "model"
    out = os.path.abspath(a["o"] or (os.path.dirname(inp) if not os.path.isdir(inp) else inp))
    os.makedirs(out, exist_ok=True)
    written = []
    if a["fmt"] in ("xml", "both"):
        p = os.path.join(out, base + ".kintecus.xml"); write_kintecus_xml(data, p); written.append(p)
    rc = 0
    if a["fmt"] in ("sbml", "both"):
        try:
            import libsbml  # noqa
        except ImportError:
            print("kintpack2sbml: SBML output needs python-libsbml:  pip install python-libsbml", file=sys.stderr)
            return 2
        conv = Kin2SBML(data, a["i"], a["max"])
        conv.build()
        p = os.path.join(out, base + ".xml"); write_sbml(conv, p); written.append(p)
        if a["sedml"] and write_sedml(conv, os.path.join(out, base + ".sedml"), base + ".xml"):
            written.append(os.path.join(out, base + ".sedml"))
        if a["report"]:
            rp = os.path.join(out, base + "_kintecus_unsupported.txt"); write_report(conv, rp); written.append(rp)
        m = conv.m
        print(f"{m.getName()}\n  SBML L3V2: {m.getNumSpecies()} species, {m.getNumReactions()} reactions, "
              f"{m.getNumEvents()} events, {m.getNumRules()} rules")
        if conv.unsupported:
            print(f"  not representable in SBML: {len(conv.unsupported)} item(s) - see the report / XML comment")
            for w, x, _ in conv.unsupported[:8]:
                print(f"    - {w}: {x[:110]}")
            if len(conv.unsupported) > 8:
                print(f"    ... and {len(conv.unsupported) - 8} more")
        if a["validate"]:
            import libsbml
            doc = libsbml.readSBMLFromFile(p)
            doc.setConsistencyChecks(libsbml.LIBSBML_CAT_UNITS_CONSISTENCY, False)
            doc.checkConsistency()
            errs = [doc.getError(k) for k in range(doc.getNumErrors())
                    if doc.getError(k).getSeverity() >= libsbml.LIBSBML_SEV_ERROR]
            print(f"  libSBML consistency check: {len(errs)} error(s)")
            for e in errs[:10]:
                print(f"    line {e.getLine()}: {e.getShortMessage()}: {e.getMessage().strip()[:160]}")
            rc = 1 if errs else 0
    for w in written:
        print(f"  wrote {w}")
    return rc


if __name__ == "__main__":
    sys.exit(main())
