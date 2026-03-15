#!/usr/bin/env python3
"""
sbml2kintecus.py  —  Convert SBML Level 2 XML files to Kintecus model.dat format.

Usage:
    python sbml2kintecus.py input.xml                  # writes to stdout
    python sbml2kintecus.py input.xml -o model.dat     # writes to file
    python sbml2kintecus.py input.xml -o model.dat -v  # verbose warnings

Output format (tab-delimited):
    <rate_constant>  <LHS==>RHS>  [# comment]

Supported SBML L2 features:
    - Integer stoichiometry (stoichiometry="N" attribute)
    - stoichiometryMath referencing global parameters (resolved to numeric value)
    - Reversible reactions split into TWO lines:
          forward:  LHS==>RHS  with k_forward
          reverse:  RHS==>LHS  with k_reverse
      Rate constants extracted from MathML kinetic law parsed as
          forward_expr - reverse_expr  (top-level <minus/>)
      Each sub-expression is a <times/> block; the <ci> that matches a local
      kineticLaw parameter (and is not a species or compartment) is taken as k.
    - Irreversible reactions written as a single LHS==>RHS line
    - Compartment variables in kinetic law (ignored — Kintecus is compartment-free)
    - Multi-compartment models: species renamed speciesID{compartmentID} throughout
      (single-compartment models are unaffected — no suffix added)
    - Local kineticLaw parameters as rate constants
    - Global model parameters (used to resolve stoichiometryMath)
    - Empty reactant / product sides  (written as NULL)
"""

import sys
import re
import math
import argparse
import xml.etree.ElementTree as ET
from pathlib import Path

# ── Namespace helpers ──────────────────────────────────────────────────────────

MATHML_NS = "http://www.w3.org/1998/Math/MathML"

def write_kintecus_file(path: str, text: str) -> None:
    """Write text to a file using Windows CRLF line endings (required by Kintecus).
    Binary mode is used to guarantee \r\n on all platforms."""
    with open(path, "wb") as f:
        f.write(text.replace("\n", "\r\n").encode("utf-8"))


def detect_sbml_ns(root_tag: str) -> str:
    """Extract namespace URI from a Clark-notation tag like '{uri}sbml'."""
    m = re.match(r"^\{(.+)\}", root_tag)
    return m.group(1) if m else ""

def make_tagger(ns: str):
    """Return a function that prepends namespace to a local tag name."""
    if ns:
        return lambda t: f"{{{ns}}}{t}"
    return lambda t: t

def mtag(t: str) -> str:
    """MathML-namespaced tag."""
    return f"{{{MATHML_NS}}}{t}"


# ── SBML level / version detection ────────────────────────────────────────────

def get_model_level_version(root) -> tuple[int, int]:
    """Return (level, version) from the SBML root element attributes.
    Defaults to (2, 1) if attributes are absent or unparseable."""
    try:
        level   = int(root.get("level",   "2"))
        version = int(root.get("version", "1"))
    except (ValueError, TypeError):
        level, version = 2, 1
    return level, version


# ── SBML L3 package / extension detection ─────────────────────────────────────

# Known L3 packages whose xmlns URIs we recognise.
# Models using these may not convert faithfully.
_L3_PACKAGES: dict[str, str] = {
    "fbc":   "Flux Balance Constraints (FBC) — constraint-based model, not an ODE system",
    "comp":  "Hierarchical Model Composition — sub-models not flattened",
    "qual":  "Qualitative Models — logical/Boolean network, not rate equations",
    "spatial": "Spatial extension — spatial PDE model not supported",
    "dyn":   "Dynamic structures extension — not supported",
    "groups":"Groups extension — purely organisational, conversion may still work",
    "render":"Render / graphical layout — ignored (cosmetic only)",
    "layout":"Layout extension — ignored (cosmetic only)",
}

def check_unsupported_packages(root, xml_path: str | None = None) -> list[str]:
    """
    Detect known SBML L3 package namespaces in the document and return warning strings.

    Two detection strategies (both needed because ET discards raw xmlns declarations):
      1. Walk every element tag in the parsed tree — catches packages that actually
         contribute elements (fbc:GeneProduct, comp:Submodel, qual:Transition …).
      2. Regex the raw XML header for xmlns:prefix="..." declarations — catches
         packages that are declared on the root but contribute no elements
         (render, layout, groups when used only as annotations).

    Returns an empty list if no known packages are found.
    """
    # Collect all namespace URIs seen in the document
    ns_uris: set[str] = set()

    # Strategy 1 — parsed tree
    for el in root.iter():
        m = re.match(r"^\{(.+)\}", el.tag)
        if m:
            ns_uris.add(m.group(1))
        for attr_key in el.attrib:
            m = re.match(r"^\{(.+)\}", attr_key)
            if m:
                ns_uris.add(m.group(1))

    # Strategy 2 — raw XML header (first 8 KB is enough for xmlns declarations)
    if xml_path:
        try:
            with open(xml_path, "rb") as fh:
                header = fh.read(8192).decode("utf-8", errors="replace")
            for uri in re.findall(r'xmlns:\w+\s*=\s*["\']([^"\']+)["\']', header):
                ns_uris.add(uri)
        except OSError:
            pass

    # URIs that are SBML core (not packages) or other standard namespaces.
    # Must NOT use the bare "http://www.sbml.org/sbml" prefix because package
    # URIs like ".../level3/version1/fbc/version2" share that prefix.
    # We match core specifically by the presence of "/core" in the path.
    def _is_core_ns(uri: str) -> bool:
        u = uri.lower()
        if "/core" in u and "sbml" in u:
            return True
        return any(u.startswith(p) for p in (
            "http://www.w3.org/1998/math/mathml",
            "http://www.w3.org/1999/xhtml",
            "http://www.w3.org/2000/xmlns",
            "http://www.w3.org/xml/1998/namespace",
            "http://www.w3.org/2001/xmlschema",
        ))

    seen_pkgs: set[str] = set()
    warnings: list[str] = []
    for ns_uri in ns_uris:
        if _is_core_ns(ns_uri):
            continue
        for pkg_key, pkg_desc in _L3_PACKAGES.items():
            if pkg_key in ns_uri.lower() and pkg_key not in seen_pkgs:
                seen_pkgs.add(pkg_key)
                warnings.append(
                    f"L3 package detected — '{pkg_key}': {pkg_desc}. "
                    f"Conversion will proceed but results may be incomplete."
                )
                break

    return warnings


# ── Constant MathML expression evaluator ──────────────────────────────────────

def _eval_mathml(el, params: dict) -> float | None:
    """
    Recursively evaluate a MathML element that represents a *constant* expression
    (i.e. no time-dependent variables).  Returns a float, or None if the expression
    cannot be fully evaluated (dynamic variable, unsupported operator, parse error).

    Supported:
      Literals : cn (plain, integer, real, e-notation, rational)
      Reference: ci  (looks up name in params)
      Constants: csymbol pi, exponentiale
      Operators: times  divide  plus  minus (unary & binary)
                 power  root  abs  exp  ln  log  floor  ceiling
      Wrapper  : math (descends to first child)
    """
    local_tag = el.tag.split("}")[-1] if "}" in el.tag else el.tag

    # ── <math> wrapper — unwrap ──
    if local_tag == "math":
        children = list(el)
        return _eval_mathml(children[0], params) if children else None

    # ── Literal number ──
    if local_tag == "cn":
        cn_type = el.get("type", "real").lower()
        text = (el.text or "").strip()
        try:
            if cn_type == "e-notation":
                # SBML uses <cn type="e-notation"> mantissa <sep/> exponent </cn>
                sep = el.find(mtag("sep"))
                if sep is not None and sep.tail:
                    mantissa = float(text)
                    exponent = float(sep.tail.strip())
                    return mantissa * (10.0 ** exponent)
                return float(text)
            elif cn_type == "rational":
                # <cn type="rational"> numerator <sep/> denominator </cn>
                sep = el.find(mtag("sep"))
                if sep is not None and sep.tail:
                    numerator   = float(text)
                    denominator = float(sep.tail.strip())
                    return numerator / denominator if denominator != 0 else None
                return float(text)
            else:
                return float(text)
        except (ValueError, TypeError, ZeroDivisionError):
            return None

    # ── Parameter reference ──
    if local_tag == "ci":
        name = (el.text or "").strip()
        if name in params:
            try:
                return float(params[name])
            except (ValueError, TypeError):
                return None
        return None          # unknown symbol — may be time, species, etc.

    # ── Named mathematical constant ──
    if local_tag == "csymbol":
        defn = el.get("definitionURL", "").lower()
        if "pi" in defn:
            return math.pi
        if "exponentiale" in defn:
            return math.e
        return None

    # ── Piecewise — too dynamic in general ──
    if local_tag == "piecewise":
        return None

    # ── Compound expression <apply> ──
    if local_tag == "apply":
        children = list(el)
        if not children:
            return None
        op_tag = children[0].tag.split("}")[-1] if "}" in children[0].tag else children[0].tag

        # Collect qualifier children (logbase, degree) separately from value operands
        qualifier_tags = {"logbase", "degree", "bvar", "uplimit", "lowlimit"}
        operand_els = [c for c in children[1:] if
                       (c.tag.split("}")[-1] if "}" in c.tag else c.tag) not in qualifier_tags]
        qualifier_els = {
            (c.tag.split("}")[-1] if "}" in c.tag else c.tag): c
            for c in children[1:] if
            (c.tag.split("}")[-1] if "}" in c.tag else c.tag) in qualifier_tags
        }

        operands = [_eval_mathml(c, params) for c in operand_els]
        if any(v is None for v in operands):
            return None

        try:
            if op_tag == "times":
                result = 1.0
                for v in operands:
                    result *= v
                return result

            if op_tag == "divide":
                if len(operands) == 2 and operands[1] != 0:
                    return operands[0] / operands[1]
                return None

            if op_tag == "plus":
                return sum(operands)

            if op_tag == "minus":
                if len(operands) == 1:
                    return -operands[0]
                if len(operands) == 2:
                    return operands[0] - operands[1]
                return None

            if op_tag == "power":
                if len(operands) == 2:
                    return operands[0] ** operands[1]
                return None

            if op_tag == "root":
                if len(operands) == 1:
                    base = 2.0
                    if "degree" in qualifier_els:
                        b = _eval_mathml(qualifier_els["degree"], params)
                        if b is None or b == 0:
                            return None
                        base = b
                    val = operands[0]
                    if val < 0 and base == 2.0:
                        return None
                    return val ** (1.0 / base)
                return None

            if op_tag == "abs" and len(operands) == 1:
                return abs(operands[0])

            if op_tag == "exp" and len(operands) == 1:
                return math.exp(operands[0])

            if op_tag == "ln" and len(operands) == 1:
                return math.log(operands[0]) if operands[0] > 0 else None

            if op_tag == "log" and len(operands) == 1:
                if "logbase" in qualifier_els:
                    base_val = _eval_mathml(qualifier_els["logbase"], params)
                    if base_val is None or base_val <= 0 or base_val == 1:
                        return None
                    return math.log(operands[0], base_val) if operands[0] > 0 else None
                return math.log10(operands[0]) if operands[0] > 0 else None

            if op_tag == "floor" and len(operands) == 1:
                return float(math.floor(operands[0]))

            if op_tag == "ceiling" and len(operands) == 1:
                return float(math.ceil(operands[0]))

            if op_tag == "factorial" and len(operands) == 1:
                n = int(operands[0])
                return float(math.factorial(n)) if n >= 0 else None

        except (ValueError, OverflowError, ZeroDivisionError):
            return None

    return None   # unrecognised tag


# ── Initial assignments (SBML L2v2+ / L3) ────────────────────────────────────

def parse_initial_assignments(model_el, stag,
                              global_params: dict) -> tuple[dict[str, str], list[str]]:
    """
    Parse <listOfInitialAssignments> and return ({symbol: value_str}, [warnings]).

    Each <initialAssignment symbol="id"> contains a MathML <math> block.
    We evaluate constant expressions; non-constant or unparseable ones are skipped
    with a warning so the XML attribute value (if any) remains in effect.

    The symbol may refer to a species, compartment, or parameter — callers decide
    which symbols are relevant (typically, only species IDs are applied to ICs).
    """
    assignments: dict[str, str] = {}
    warnings: list[str] = []

    loia = model_el.find(stag("listOfInitialAssignments"))
    if loia is None:
        return assignments, warnings

    for ia in loia.findall(stag("initialAssignment")):
        symbol = ia.get("symbol", "").strip()
        if not symbol:
            warnings.append("initialAssignment with no 'symbol' attribute — skipped")
            continue

        math_el = ia.find(mtag("math"))
        if math_el is None:
            warnings.append(
                f"initialAssignment for '{symbol}': no <math> element — "
                f"XML attribute value (if any) will be used"
            )
            continue

        val = _eval_mathml(math_el, global_params)
        if val is None:
            warnings.append(
                f"initialAssignment for '{symbol}': expression could not be evaluated "
                f"statically (may reference species, time, or unsupported operators) — "
                f"XML attribute value (if any) will be used"
            )
        else:
            # Use %g to avoid ugly scientific notation for normal concentrations
            # but preserve it for very small/large values.
            assignments[symbol] = f"{val:g}"

    if assignments:
        count = len(assignments)
        warnings.append(
            f"Applied {count} initialAssignment(s) to override XML attribute values: "
            + ", ".join(f"{k}={v}" for k, v in assignments.items())
        )

    return assignments, warnings


# ── Stoichiometry helpers ──────────────────────────────────────────────────────

def resolve_stoichiometry(spec_ref_el, global_params: dict, stag) -> tuple[float, str | None]:
    """
    Return (stoich_value, warning_or_None) for a speciesReference element.
    Handles: stoichiometryMath/<ci>, stoichiometry attribute, or defaults to 1.
    """
    # 1. stoichiometryMath (SBML L2 — math expression for stoich)
    stoich_math = spec_ref_el.find(stag("stoichiometryMath"))
    if stoich_math is not None:
        math_el = stoich_math.find(mtag("math"))
        if math_el is not None:
            ci = math_el.find(mtag("ci"))
            if ci is not None:
                param_id = ci.text.strip()
                if param_id in global_params:
                    val = float(global_params[param_id])
                    return val, f"stoich resolved from global param {param_id}={val}"
                else:
                    return 1.0, (
                        f"stoichiometryMath references unknown param '{param_id}'; "
                        f"defaulting to 1"
                    )
        return 1.0, "stoichiometryMath present but unreadable; defaulting to 1"

    # 2. Plain numeric attribute
    stoich_attr = spec_ref_el.get("stoichiometry")
    if stoich_attr is not None:
        return float(stoich_attr), None

    # 3. Default
    return 1.0, None


def side_string(spec_refs, global_params: dict, stag,
                species_rename: dict,
                is_reactant: bool = False) -> tuple[str, list[str]]:
    """
    Build a '+'-joined species string (with repetition for stoichiometry > 1)
    from a list of speciesReference elements.
    species_rename maps original species IDs to their Kintecus display name
    (e.g. "C" -> "C{BZ1}" in multi-compartment models; identity otherwise).
    is_reactant: when True and the side is empty, returns "0 NULL" to signal
                 a zero-order reaction to Kintecus (rather than plain "NULL").
    Returns (string, [warnings]).
    """
    parts = []
    warnings = []
    for sr in spec_refs:
        sp_id = sr.get("species")
        display = species_rename.get(sp_id, sp_id)  # renamed or unchanged
        stoich, warn = resolve_stoichiometry(sr, global_params, stag)
        if warn:
            warnings.append(f"  [{sp_id}] {warn}")
        if stoich != int(stoich):
            # Fractional stoichiometry — Kintecus accepts a numeric prefix, e.g. 0.5Br
            parts.append(f"{stoich}{display}")
        else:
            n = int(stoich)
            parts.append("+".join([display] * n))  # e.g. HBrO2+HBrO2 for stoich=2
    if parts:
        return "+".join(parts), warnings
    # Empty side — zero-order source uses "0 NULL", empty product sink uses "NULL"
    return ("0 NULL" if is_reactant else "NULL"), warnings


# ── Rate-constant extraction ───────────────────────────────────────────────────

def _load_local_params(kinetic_law_el, stag) -> dict:
    """Return {param_id: {"value": ..., "name": ...}} from kineticLaw/listOfParameters."""
    local_params = {}
    lop = kinetic_law_el.find(stag("listOfParameters"))
    if lop is not None:
        for p in lop.findall(stag("parameter")):
            pid  = p.get("id")
            pval = p.get("value", "1.0")
            pname = p.get("name", pid)
            local_params[pid] = {"value": pval, "name": pname}
    return local_params


def _pick_rate_from_params(local_params: dict, warnings: list) -> tuple[str, str]:
    """
    Given a dict of local params, choose the one that is the rate constant.
    Priority: single param → only option; else prefer id starting with 'k';
    else fall back to first with a warning.
    Returns (value_str, display_name).
    """
    if not local_params:
        warnings.append("No local kineticLaw parameters found; rate constant set to 1.0")
        return "1.0", "?"

    if len(local_params) == 1:
        pid, pdata = next(iter(local_params.items()))
        return pdata["value"], pdata["name"]

    k_candidates = [pid for pid in local_params if pid.lower().startswith("k")]
    if len(k_candidates) == 1:
        pid = k_candidates[0]
        others = [p for p in local_params if p != pid]
        warnings.append(
            f"Multiple local params; selected '{pid}' as rate constant (others: {others})"
        )
        return local_params[pid]["value"], local_params[pid]["name"]

    pid, pdata = next(iter(local_params.items()))
    warnings.append(
        f"Multiple local params; arbitrarily using first '{pid}={pdata['value']}' — please verify"
    )
    return pdata["value"], pdata["name"]


def _find_rate_in_times_expr(apply_el,
                             local_params: dict,
                             species_ids: set,
                             compartment_ids: set,
                             warnings: list,
                             global_params: dict | None = None) -> tuple[str, str]:
    """
    Scan a MathML <apply><times/>...</apply> block for the rate-constant term.

    Search priority:
      1. A <ci> matching a local kineticLaw parameter.
      2. A <cn> literal number.
      3. A <ci> matching a global model parameter that is neither a species
         nor a compartment — preferring names that start with 'k'.
    Returns (value_str, display_name).
    """
    children = list(apply_el)   # first child is the operator element
    operands = children[1:]     # everything after the operator

    # ── 1. Local params — direct children first, then full subtree ──
    for op in operands:
        local_tag = op.tag.split("}")[-1] if "}" in op.tag else op.tag
        if local_tag == "ci":
            name = op.text.strip()
            if name in local_params:
                return local_params[name]["value"], local_params[name]["name"]
        elif local_tag == "cn":
            val = op.text.strip()
            return val, val

    for ci in apply_el.iter(mtag("ci")):
        name = ci.text.strip()
        if name in local_params:
            return local_params[name]["value"], local_params[name]["name"]

    # ── 2. Global params fallback ──
    if global_params:
        candidates = []
        for ci in apply_el.iter(mtag("ci")):
            name = ci.text.strip()
            if (name in global_params
                    and name not in species_ids
                    and name not in compartment_ids):
                candidates.append(name)

        if candidates:
            # Prefer names starting with 'k' (rate constants vs scalar multipliers)
            k_cands = [n for n in candidates if n.lower().startswith("k")]
            chosen = k_cands[0] if k_cands else candidates[0]
            if not k_cands:
                warnings.append(
                    f"No k-prefixed global param in expression; "
                    f"using '{chosen}' as rate constant — please verify"
                )
            return global_params[chosen], chosen

    # ── 3. Nothing found — last resort ──
    warnings.append(
        "Could not identify rate constant in <times> expression; "
        "falling back to first local param"
    )
    return _pick_rate_from_params(local_params, warnings)


def _unwrap_compartment_times(top_apply):
    """
    Many SBML models wrap the net rate in an outer <times/>(compartment, expr).
    If top_apply is a <times/> and one of its operands is an <apply>, return
    that inner <apply> (the actual rate expression).  Otherwise return top_apply
    unchanged.
    This handles both:
        <times/> <ci>comp</ci> <apply><minus/>...</apply>   (reversible)
        <times/> <ci>comp</ci> <ci>k</ci> ...              (irreversible — unchanged)
    """
    children = list(top_apply)
    if not children:
        return top_apply
    op_tag = children[0].tag.split("}")[-1] if "}" in children[0].tag else children[0].tag
    if op_tag != "times":
        return top_apply
    # Look for a nested <apply> among the operands
    for child in children[1:]:
        ctag = child.tag.split("}")[-1] if "}" in child.tag else child.tag
        if ctag == "apply":
            inner_children = list(child)
            if inner_children:
                inner_op = (inner_children[0].tag.split("}")[-1]
                            if "}" in inner_children[0].tag
                            else inner_children[0].tag)
                # Only unwrap if the inner apply is <minus/> (reversible net rate)
                if inner_op == "minus":
                    return child
    return top_apply


def extract_irreversible_rate(kinetic_law_el, species_ids: set,
                              compartment_ids: set, stag,
                              global_params: dict | None = None) -> tuple[str, str, list]:
    """
    Extract a single rate constant from an irreversible reaction's kineticLaw.

    Strategy:
      1. Local kineticLaw parameters (fast path for simple models).
      2. MathML scan using global model parameters as fallback (for models
         where all parameters are defined at the model level).
    Returns (value_str, display_name, warnings).
    """
    warnings = []
    local_params = _load_local_params(kinetic_law_el, stag)

    # Fast path: local params present — use existing logic
    if local_params:
        val, name = _pick_rate_from_params(local_params, warnings)
        return val, name, warnings

    # No local params — scan the MathML with global params
    if global_params:
        math_el = kinetic_law_el.find(mtag("math"))
        if math_el is not None:
            top_apply = math_el.find(mtag("apply"))
            if top_apply is not None:
                # Use the unwrapped apply (handles outer compartment×rate pattern)
                effective_apply = _unwrap_compartment_times(top_apply)
                val, name = _find_rate_in_times_expr(
                    effective_apply, local_params, species_ids,
                    compartment_ids, warnings, global_params
                )
                return val, name, warnings

    warnings.append("No local kineticLaw parameters and no MathML found; rate set to 1.0")
    return "1.0", "?", warnings


def extract_reversible_rates(kinetic_law_el, species_ids: set,
                             compartment_ids: set, stag,
                             global_params: dict | None = None) -> tuple[str, str, str, str, list]:
    """
    Parse a reversible reaction's kineticLaw MathML.

    Handles two common top-level structures:
      A)  <apply><minus/> fwd_expr rev_expr </apply>
      B)  <apply><times/> <ci>compartment</ci>
                <apply><minus/> fwd_expr rev_expr </apply>
          </apply>
    Structure B is unwrapped to A automatically.

    Rate constants are resolved from local kineticLaw parameters first,
    then from global model parameters (excluding species and compartments,
    preferring k-prefixed names).

    Returns (fwd_val, fwd_name, rev_val, rev_name, warnings).
    """
    warnings = []
    local_params = _load_local_params(kinetic_law_el, stag)

    math_el = kinetic_law_el.find(mtag("math"))
    if math_el is None:
        warnings.append("No <math> element in kineticLaw; cannot split reversible rates")
        val, name = _pick_rate_from_params(local_params, warnings)
        return val, name, "1.0", "?", warnings

    top_apply = math_el.find(mtag("apply"))
    if top_apply is None:
        warnings.append("No top-level <apply> in MathML; cannot split reversible rates")
        val, name = _pick_rate_from_params(local_params, warnings)
        return val, name, "1.0", "?", warnings

    # Unwrap outer <times/>(compartment, <minus/>(...)) if present
    minus_apply = _unwrap_compartment_times(top_apply)
    children = list(minus_apply)

    if not children:
        warnings.append("Empty <apply>; cannot split reversible rates")
        val, name = _pick_rate_from_params(local_params, warnings)
        return val, name, "1.0", "?", warnings

    op_tag = children[0].tag.split("}")[-1] if "}" in children[0].tag else children[0].tag

    if op_tag != "minus":
        warnings.append(
            f"Kinetic law operator is <{op_tag}>, not <minus/>; "
            f"treating entire expression as forward rate and setting reverse to 1.0"
        )
        fwd_apply = minus_apply if op_tag == "times" else top_apply
        fwd_val, fwd_name = _find_rate_in_times_expr(
            fwd_apply, local_params, species_ids, compartment_ids, warnings, global_params
        )
        return fwd_val, fwd_name, "1.0", "?", warnings

    # We have <minus/> fwd_expr rev_expr
    operands = children[1:]
    if len(operands) < 2:
        warnings.append("<minus/> has fewer than 2 operands; cannot extract both rates")
        val, name = _pick_rate_from_params(local_params, warnings)
        return val, name, "1.0", "?", warnings

    fwd_expr, rev_expr = operands[0], operands[1]

    # ── Forward rate ──
    fwd_op_tag = fwd_expr.tag.split("}")[-1] if "}" in fwd_expr.tag else fwd_expr.tag
    if fwd_op_tag == "apply":
        fwd_val, fwd_name = _find_rate_in_times_expr(
            fwd_expr, local_params, species_ids, compartment_ids, warnings, global_params
        )
    elif fwd_op_tag == "ci":
        name = fwd_expr.text.strip()
        if name in local_params:
            fwd_val, fwd_name = local_params[name]["value"], local_params[name]["name"]
        elif global_params and name in global_params:
            fwd_val, fwd_name = global_params[name], name
        else:
            fwd_val, fwd_name = name, name
    else:
        warnings.append(f"Forward expression is <{fwd_op_tag}>; expected <apply>")
        fwd_val, fwd_name = _pick_rate_from_params(local_params, warnings)

    # ── Reverse rate ──
    # Exclude the already-found forward param so the reverse search picks a different one
    rev_local = {k: v for k, v in local_params.items() if v["value"] != fwd_val}
    if not rev_local:
        rev_local = local_params

    # For global params: exclude the forward param name
    rev_global = (
        {k: v for k, v in global_params.items() if k != fwd_name}
        if global_params else None
    )

    rev_op_tag = rev_expr.tag.split("}")[-1] if "}" in rev_expr.tag else rev_expr.tag
    if rev_op_tag == "apply":
        rev_val, rev_name = _find_rate_in_times_expr(
            rev_expr, rev_local, species_ids, compartment_ids, warnings, rev_global
        )
    elif rev_op_tag == "ci":
        name = rev_expr.text.strip()
        if name in local_params:
            rev_val, rev_name = local_params[name]["value"], local_params[name]["name"]
        elif global_params and name in global_params:
            rev_val, rev_name = global_params[name], name
        else:
            rev_val, rev_name = name, name
    else:
        warnings.append(f"Reverse expression is <{rev_op_tag}>; expected <apply>")
        rev_val, rev_name = "1.0", "?"

    return fwd_val, fwd_name, rev_val, rev_name, warnings


# ── MathML → Kintecus USER expression ────────────────────────────────────────
#
# Kintecus V6.51+ accepts a USER-defined rate function in MODEL.DAT column 1:
#   [USER;var1,var2,...:(expression)]
#
# Allowed variable names:
#   T           current temperature (K)
#   P           pressure
#   M           Loschmidt's value (total gas concentration)
#   R           gas constant
#   cSPECIES    concentration of species SPECIES  (e.g. cH2O, cHBrO2)
#
# The USER expression replaces the scalar rate constant k in column 1.
# Because Kintecus still multiplies column-1 by the reactant concentrations,
# the expression stored here is:  full_rate / (cA * cB * ...)
# so that:  USER_value × [A] × [B] = full_rate  ✓
#
# For zero-order reactions (no reactants) the expression IS the full rate.
# Simple mass-action laws (no complex operators, no T/P/M/R) fall back to
# the plain numeric rate-constant format — the USER path is not taken.

_KINTECUS_OP_FUNCS: dict[str, str] = {
    "abs":     "abs",
    "exp":     "exp",
    "ln":      "log",    # Kintecus 'log' = natural logarithm
    "sin":     "sin",
    "cos":     "cos",
    "tan":     "tan",
    "arctan":  "atn",    # Kintecus uses 'atn' (not 'atan')
    "arcsin":  "asin",
    "arccos":  "acos",
    "sinh":    "sinh",
    "cosh":    "cosh",
    "tanh":    "tanh",
    "floor":   "floor",
    "ceiling": "floor",  # ceiling → floor (Kintecus V6 approximation; warns)
    "sqrt":    "sqrt",
}

# Regex that matches operators / functions which require USER syntax.
# Division '/' is the key signal (MM denominators, Hill functions, etc.).
# Power '^' is deliberately excluded — it appears in pure mass-action kinetics
# whenever a reactant has stoichiometry > 1 (e.g. k*[A]^2), and the existing
# rate-constant extractor handles those correctly without USER.
_COMPLEX_EXPR_RE = re.compile(
    r'/'
    r'|(?<![a-zA-Z])'
    r'(log10?|exp|sin|cos|tan|atn|asin|acos|sinh|cosh|tanh|floor|ceil|abs|sqrt)'
    r'\s*\('
)


def _cn_to_str(el) -> str:
    """Format a MathML <cn> element as a Kintecus-compatible number string."""
    cn_type = el.get("type", "real").lower()
    text    = (el.text or "").strip()
    try:
        sep = el.find(mtag("sep"))
        if cn_type == "e-notation" and sep is not None and sep.tail:
            val = float(text) * (10.0 ** float(sep.tail.strip()))
        elif cn_type == "rational" and sep is not None and sep.tail:
            d   = float(sep.tail.strip())
            val = float(text) / d if d != 0 else float(text)
        else:
            val = float(text)
        return f"{val:g}"
    except (ValueError, TypeError):
        return text or "1"


def _mathml_to_kintecus_expr(
    el,
    species_ids:      set,
    compartment_ids:  set,
    params:           dict,    # merged global + local params {id: value_str}
    compartment_sizes: dict,   # {id: size_str}
    species_rename:   dict,
    used_vars:        set,     # OUTPUT — Kintecus vars referenced (T,P,M,R,cXXX)
    substitute_params: bool = True,
) -> str:
    """
    Recursively convert a MathML element to a Kintecus USER expression string.

    Raises ValueError with a descriptive message for any unsupported construct
    (piecewise, 'time' csymbol, unknown symbols, unsupported operators).
    """
    def qtag(e) -> str:
        return e.tag.split("}")[-1] if "}" in e.tag else e.tag

    local_tag = qtag(el)

    # ── <math> wrapper ──
    if local_tag == "math":
        kids = list(el)
        if not kids:
            raise ValueError("empty <math>")
        return _mathml_to_kintecus_expr(
            kids[0], species_ids, compartment_ids, params,
            compartment_sizes, species_rename, used_vars, substitute_params
        )

    # ── Literal number ──
    if local_tag == "cn":
        return _cn_to_str(el)

    # ── Variable / parameter reference ──
    if local_tag == "ci":
        name = (el.text or "").strip()
        # Species concentration → cDISPLAY
        if name in species_ids:
            display = species_rename.get(name, name)
            var = f"c{display}"
            used_vars.add(var)
            return var
        # Compartment → substitute its numeric size
        if name in compartment_ids:
            return compartment_sizes.get(name, "1")
        # Kintecus built-in scalars (checked BEFORE param lookup so T/P/R are
        # recognised even if they also appear as SBML parameters)
        nl = name.lower()
        if nl in ("t", "temperature", "temp"):
            used_vars.add("T")
            return "T"
        if nl in ("p", "pressure"):
            used_vars.add("P")
            return "P"
        if nl in ("r", "rgas", "r_gas", "gasconstant"):
            used_vars.add("R")
            return "R"
        # Named parameter → substitute value
        if name in params:
            if substitute_params:
                try:
                    return f"{float(params[name]):g}"
                except (ValueError, TypeError):
                    return str(params[name])
            return name   # display mode — keep name
        raise ValueError(f"unresolvable symbol '{name}'")

    # ── Named mathematical constant ──
    if local_tag == "csymbol":
        defn = el.get("definitionURL", "").lower()
        if "pi" in defn:
            return "3.14159265358979"
        if "exponentiale" in defn:
            return "2.71828182845905"
        if "avogadro" in defn:
            used_vars.add("M")
            return "M"
        if "time" in defn:
            raise ValueError("'time' csymbol is dynamic — not convertible to a USER constant")
        raise ValueError(f"unknown csymbol '{el.get('definitionURL', '')}'")

    # ── Piecewise — too dynamic ──
    if local_tag == "piecewise":
        raise ValueError("piecewise expressions are not supported in USER")

    # ── Compound <apply> ──
    if local_tag != "apply":
        raise ValueError(f"unsupported MathML element <{local_tag}>")

    children = list(el)
    if not children:
        raise ValueError("empty <apply>")

    op_tag = qtag(children[0])
    qualifier_tags = {"logbase", "degree", "bvar", "uplimit", "lowlimit"}
    qualifiers   = {qtag(c): c for c in children[1:] if qtag(c) in qualifier_tags}
    operand_els  = [c for c in children[1:] if qtag(c) not in qualifier_tags]

    def sub(e) -> str:
        return _mathml_to_kintecus_expr(
            e, species_ids, compartment_ids, params,
            compartment_sizes, species_rename, used_vars, substitute_params
        )

    def paren_add(s: str, trigger=r'[+\-]') -> str:
        """Wrap s in parens if it contains low-precedence operators."""
        return f"({s})" if re.search(trigger, s) else s

    # ── n-ary product ──
    if op_tag == "times":
        parts = [paren_add(sub(o)) for o in operand_els]
        return "*".join(parts) if parts else "1"

    # ── Division ──
    if op_tag == "divide":
        if len(operand_els) != 2:
            raise ValueError(f"<divide/> expects 2 operands, got {len(operand_els)}")
        num = sub(operand_els[0])
        den = sub(operand_els[1])
        return f"{paren_add(num)}/{paren_add(den, r'[+\-/*]')}"

    # ── n-ary sum ──
    if op_tag == "plus":
        return "+".join(sub(o) for o in operand_els)

    # ── Subtraction / unary negation ──
    if op_tag == "minus":
        if len(operand_els) == 1:
            s = sub(operand_els[0])
            return f"-({s})" if re.search(r'[+\-]', s) else f"-{s}"
        if len(operand_els) == 2:
            lhs_s = sub(operand_els[0])
            rhs_s = sub(operand_els[1])
            return f"{lhs_s}-{paren_add(rhs_s)}"
        raise ValueError(f"<minus/> with {len(operand_els)} operands not supported")

    # ── Power — Kintecus requires x^(-1*n) for negative exponents ──
    if op_tag == "power":
        if len(operand_els) != 2:
            raise ValueError(f"<power/> expects 2 operands, got {len(operand_els)}")
        base = sub(operand_els[0])
        exp  = sub(operand_els[1])
        base_s = paren_add(base, r'[+\-/*]')
        # Detect negative literal exponent (e.g. "-0.72") and reformat for Kintecus
        neg_lit = re.match(r'^-([\d.e+\-]+)$', exp.strip())
        if neg_lit:
            return f"{base_s}^(-1*{neg_lit.group(1)})"
        # Detect unary minus wrapping an expression (e.g. "-(a/b)")
        if re.match(r'^-\(', exp.strip()):
            inner = exp.strip()[2:-1]   # strip "-(" and ")"
            return f"{base_s}^(-1*({inner}))"
        return f"{base_s}^({exp})"

    # ── Root ──
    if op_tag == "root":
        if len(operand_els) != 1:
            raise ValueError("<root/> expects 1 operand")
        val = sub(operand_els[0])
        if "degree" in qualifiers:
            deg = sub(qualifiers["degree"])
            return f"{paren_add(val, r'[+\-/*]')}^(1/{deg})"
        return f"sqrt({val})"

    # ── Standard single-arg functions ──
    if op_tag in _KINTECUS_OP_FUNCS and len(operand_els) == 1:
        if op_tag == "ceiling":
            pass   # fall through — warn below after building the string
        fn = _KINTECUS_OP_FUNCS[op_tag]
        return f"{fn}({sub(operand_els[0])})"

    # ── log (base-10 default, or custom base with <logbase>) ──
    if op_tag == "log" and len(operand_els) == 1:
        val = sub(operand_els[0])
        if "logbase" in qualifiers:
            base_s = sub(qualifiers["logbase"])
            return f"log({val})/log({base_s})"
        return f"log10({val})"

    raise ValueError(f"unsupported MathML operator <{op_tag}> with {len(operand_els)} operand(s)")


def _species_product_str(
    r_refs, global_params: dict, stag, species_rename: dict
) -> str:
    """
    Build the Kintecus expression string for the product of reactant concentrations.
    e.g. A (stoich 1) + B (stoich 2) → "cA*cB*cB"
    Returns "" for zero-order reactions (empty reactant list).
    """
    parts: list[str] = []
    for sr in r_refs:
        sp_id = sr.get("species")
        if not sp_id:
            continue
        display = species_rename.get(sp_id, sp_id)
        var     = f"c{display}"
        stoich, _ = resolve_stoichiometry(sr, global_params, stag)
        if stoich == int(stoich):
            parts.extend([var] * int(stoich))
        else:
            parts.append(f"{var}^({stoich:g})")
    return "*".join(parts)


def _split_reversible_math(kl_el, stag):
    """
    Extract the forward and reverse MathML <apply> sub-elements from a
    reversible kineticLaw that uses the standard <minus/> pattern:
        net = fwd_expr - rev_expr
    (with optional outer compartment×rate wrapper stripped first).
    Returns (fwd_apply, rev_apply) or (None, None) if structure is not recognised.
    """
    if kl_el is None:
        return None, None
    math_el = kl_el.find(mtag("math"))
    if math_el is None:
        return None, None
    top = math_el.find(mtag("apply"))
    if top is None:
        return None, None
    minus_el = _unwrap_compartment_times(top)
    kids = list(minus_el)
    if not kids:
        return None, None
    op = kids[0].tag.split("}")[-1] if "}" in kids[0].tag else kids[0].tag
    if op != "minus" or len(kids) < 3:
        return None, None
    return kids[1], kids[2]


def _kinetic_apply_to_user_rate(
    apply_el,
    r_refs:           list,
    species_ids:      set,
    compartment_ids:  set,
    params:           dict,
    compartment_sizes: dict,
    species_rename:   dict,
    global_params:    dict,
    stag,
) -> tuple[str | None, list[str]]:
    """
    Convert a single MathML <apply> element (one direction of a kinetic law)
    into a Kintecus rate-column string using the USER function.

    The returned string is the full_rate_expression divided by the reactant
    concentration product, so that:
        USER_value × [A] × [B] × ... = full_kinetic_rate   ✓

    Returns (rate_col_str, warnings):
      - rate_col_str is None   → fall back to existing simple k extraction
      - "[USER;vars:(expr)]"   → complex rate; use this in MODEL.DAT column 1
      - "numeric_string"       → constant (evaluated at parse time)
    """
    warnings: list[str] = []
    used_vars: set[str] = set()

    try:
        full_expr = _mathml_to_kintecus_expr(
            apply_el, species_ids, compartment_ids, params,
            compartment_sizes, species_rename, used_vars
        )
    except ValueError as e:
        warnings.append(f"USER conversion failed: {e}")
        return None, warnings

    # Simple mass-action check: no T/P/M/R variables and no complex operators
    non_conc_vars = {v for v in used_vars if not v.startswith("c")}
    if not non_conc_vars and not _COMPLEX_EXPR_RE.search(full_expr):
        # Pure mass-action — existing rate-constant extraction handles this better
        return None, warnings

    # Build denominator: product of reactant concentrations
    reactant_prod = _species_product_str(r_refs, global_params, stag, species_rename)

    # USER expression = full_rate / (cA * cB * ...)
    if reactant_prod:
        user_expr = f"({full_expr})/({reactant_prod})"
    else:
        user_expr = full_expr   # zero-order: no division needed

    # If no runtime variables at all, evaluate to a plain constant
    if not used_vars:
        val = _eval_mathml(apply_el, params)
        if val is not None:
            return f"{val:g}", warnings
        return user_expr, warnings

    var_list = ",".join(sorted(used_vars))
    return f"[USER;{var_list}:({user_expr})]", warnings


def get_rate_column(
    kl_el,
    r_refs:           list,
    p_refs:           list,
    reversible:       bool,
    species_ids:      set,
    compartment_ids:  set,
    merged_params:    dict,
    compartment_sizes: dict,
    species_rename:   dict,
    global_params:    dict,
    stag,
) -> tuple:
    """
    Determine the rate-column string(s) for one reaction, preferring USER syntax
    for complex kinetic laws and falling back to simple k extraction otherwise.

    Returns:
        irreversible → ((rate_val, rate_name), None, warnings)
        reversible   → ((fwd_val, fwd_name), (rev_val, rev_name), warnings)
    """
    all_warnings: list[str] = []

    def _try_user(apply_el, refs):
        return _kinetic_apply_to_user_rate(
            apply_el, refs, species_ids, compartment_ids,
            merged_params, compartment_sizes, species_rename, global_params, stag
        )

    if not reversible:
        if kl_el is not None:
            math_el = kl_el.find(mtag("math"))
            if math_el is not None:
                top = math_el.find(mtag("apply"))
                if top is not None:
                    rate_col, w = _try_user(_unwrap_compartment_times(top), r_refs)
                    all_warnings.extend(w)
                    if rate_col is not None:
                        return (rate_col, "USER"), None, all_warnings
            # Fall back to existing extraction
            rate_val, rate_name, kw = extract_irreversible_rate(
                kl_el, species_ids, compartment_ids, stag, global_params
            )
            all_warnings.extend(kw)
        else:
            rate_val, rate_name = "1.0", "?"
            all_warnings.append("No <kineticLaw> element; rate set to 1.0")
        return (rate_val, rate_name), None, all_warnings

    else:   # reversible
        fwd_apply, rev_apply = _split_reversible_math(kl_el, stag)
        fwd_col = rev_col = None

        if fwd_apply is not None:
            fwd_col, w = _try_user(fwd_apply, r_refs)
            all_warnings.extend(w)
        if rev_apply is not None:
            # For the reverse direction the products become the "reactants"
            rev_col, w = _try_user(rev_apply, p_refs)
            all_warnings.extend(w)

        if fwd_col is not None and rev_col is not None:
            return (fwd_col, "USER"), (rev_col, "USER"), all_warnings

        # Either USER failed or expression is simple — use existing extraction
        if kl_el is not None:
            fv, fn, rv, rn, kw = extract_reversible_rates(
                kl_el, species_ids, compartment_ids, stag, global_params
            )
            all_warnings.extend(kw)
        else:
            fv, fn, rv, rn = "1.0", "?", "1.0", "?"
            all_warnings.append("No <kineticLaw> element; both rates set to 1.0")
        return (fv, fn), (rv, rn), all_warnings


# ── Non-kinetic rules → comment lines ────────────────────────────────────────

def collect_rule_comments(
    model_el,
    stag,
    species_ids:      set,
    compartment_ids:  set,
    global_params:    dict,
    compartment_sizes: dict,
    species_rename:   dict,
) -> list[str]:
    """
    Parse <listOfRules> (assignmentRule, rateRule, algebraicRule) and return
    a list of '#'-prefixed comment lines for inclusion in MODEL.DAT.

    These rule types define algebraic or ODE constraints that Kintecus cannot
    represent directly.  Outputting them as comments preserves the information
    for the modeller.

    Returns an empty list if no rules are present.
    """
    lor = model_el.find(stag("listOfRules"))
    if lor is None:
        return []

    _RULE_LABEL = {
        "assignmentRule": "AssignmentRule",
        "rateRule":       "RateRule",
        "algebraicRule":  "AlgebraicRule",
    }

    rule_lines: list[str] = []

    for el in lor:
        local_tag = el.tag.split("}")[-1] if "}" in el.tag else el.tag
        label = _RULE_LABEL.get(local_tag)
        if label is None:
            continue
        symbol   = el.get("variable") or el.get("symbol") or "?"
        math_el  = el.find(mtag("math"))

        if math_el is None:
            rule_lines.append(f"# {label}: {symbol} = [no <math> element]")
            continue

        used_vars: set[str] = set()
        try:
            expr = _mathml_to_kintecus_expr(
                math_el, species_ids, compartment_ids, global_params,
                compartment_sizes, species_rename, used_vars,
                substitute_params=False   # keep names for readability
            )
            if label == "RateRule":
                rule_lines.append(f"# {label}: d({symbol})/dt = {expr}")
            elif label == "AlgebraicRule":
                rule_lines.append(f"# {label}: 0 = {expr}")
            else:
                rule_lines.append(f"# {label}: {symbol} = {expr}")
        except ValueError as e:
            rule_lines.append(
                f"# {label}: {symbol} = [expression not convertible: {e}]"
            )

    if not rule_lines:
        return []

    header = [
        "#",
        "# ── Non-kinetic rules (not representable in Kintecus — for reference only) ──",
    ]
    return header + rule_lines


# ── Notes extraction ──────────────────────────────────────────────────────────

def extract_notes(model_el, stag) -> list[str]:
    """
    Extract human-readable text from the SBML <notes> element and return it
    as a list of '#'-prefixed comment lines, each capped at 700 characters.

    The notes section typically contains XHTML (<body>, <p>, <table>, <h1> …).
    We walk every element with itertext() to collect raw text, then:
      - split on any embedded newlines
      - strip leading/trailing whitespace from each fragment
      - drop blank fragments
      - prefix each surviving line with '# '
      - hard-truncate at 700 characters (including the '# ' prefix)
    """
    MAX_LEN = 700

    notes_el = model_el.find(stag("notes"))
    if notes_el is None:
        return []

    # itertext() walks the whole subtree and yields text + tail strings
    raw_fragments: list[str] = []
    for text in notes_el.itertext():
        # Split on any embedded newlines so each physical line is separate
        for part in text.splitlines():
            part = part.strip()
            if part:
                raw_fragments.append(part)

    if not raw_fragments:
        return []

    comment_lines: list[str] = []
    comment_lines.append("# " + "-" * 68)
    comment_lines.append("# Notes:")
    for fragment in raw_fragments:
        line = f"# {fragment}"
        # Hard-truncate at MAX_LEN characters
        if len(line) > MAX_LEN:
            line = line[:MAX_LEN]
        comment_lines.append(line)
    comment_lines.append("# " + "-" * 68)

    return comment_lines


# ── Species.dat generation ────────────────────────────────────────────────────

def generate_species_dat(species_info: dict, species_rename: dict,
                         output_path: str | None = None) -> str:
    """
    Generate a Kintecus species.dat file from the parsed species data.

    Columns (tab-delimited):
        Species   ResidenceTime   InitialConc   DisplayOutput   ExternalConc
        SpecialDirectives   ConstantFile   Comments

    species_rename maps original species IDs to their Kintecus display names
    (populated only for multi-compartment models; identity otherwise).
    """
    out: list[str] = []

    # Header (mirrors the reference format exactly)
    out.append("#\tSpecies Description Spreadsheet")
    out.append("# Species\tResidence\tInitial\tDisplay Output\t  External\tSpecies Special\tConstant File?\t")
    out.append("#     \tTime in CSTR(s)\tConc.\t(Y/N) ?\tConc.\tDirectives ? (/N)\t(Filename/#/No)\tComments")

    # Determine column width for species name alignment (match reference right-alignment)
    display_names = [species_rename.get(sid, sid) for sid in species_info]
    name_width = max((len(n) for n in display_names), default=8)
    name_width = max(name_width, 8)  # minimum width of 8

    for sid, sdata in species_info.items():
        display = species_rename.get(sid, sid)
        ic      = sdata["ic"]
        # Right-align the species name to match the reference style
        name_col = display.rjust(name_width)
        row = f"{name_col}\t0\t{ic}\tY\t0\tNo\tNo\t"
        out.append(row)

    out.append("END")
    out.append("")  # trailing newline

    output_text = "\n".join(out)

    if output_path:
        write_kintecus_file(output_path, output_text)
        print(f"Written to: {output_path}")
    else:
        print(output_text)

    return output_text

def generate_parm_dat(global_params: dict, output_path: str | None = None) -> str:
    """
    Generate a Kintecus PARM.DAT file.

    Values extracted from SBML global parameters where recognisable, otherwise
    the following defaults apply:
        Temperature      : 298 K
        Ea units         : KCAL
        Simulation length: 10 minutes
        Accuracy         : 1e-8
        Conc. units      : MOLES/LITER

    Format is tab-delimited, matching the reference PARM.DAT layout exactly.
    """

    # ── Extract or default each field ────────────────────────────────────────

    # Temperature — look for a global param whose id/name suggests temperature
    temperature = "298"
    for pid, val in global_params.items():
        if pid.lower() in ("t", "temperature", "temp"):
            temperature = val
            break

    # Concentration units — look for a hint in global params, else default
    conc_units = "MOLES/LITER"
    for pid in global_params:
        if "molecule" in pid.lower() or "molec" in pid.lower():
            conc_units = "MOLECULES/CC"
            break

    # Ea units — always default unless a global param says otherwise
    ea_units = "KCAL"

    # Simulation length — default 10 minutes
    sim_days    = "0"
    sim_hours   = "0"
    sim_minutes = "10"
    sim_seconds = "0"
    sim_pico    = "0"

    # Accuracy
    accuracy = "0.00000001"

    # ── Build output ─────────────────────────────────────────────────────────
    out: list[str] = []

    out.append("#\tParameter Description SpreadSheet\t\t\t\t")
    out.append("# See .pdf file for an explanation of each field\t\t\t\t\t")
    out.append("#\t\t\t\t\t")
    out.append("# Starting Integration Time\tMaximum Integration Time\t"
               "Ea UNITS (Kcal, KJ,J,CAL,K)\tConc. Units(moles/L, moles/cc, molecules/cc)\tX0\t")
    out.append(f"0.000001\t1\t{ea_units}\t{conc_units}\t0\t")
    out.append("#\t\t\tExternal Heat Source/Sink # OR Profile (Filename) \t\t")
    out.append("#Temperature (K) or Filename\tPressure (Constant ? (Yes/No) )\t"
               "Volume Profile (Filename/No)\t"
               "OR Condutance:Extern. Temp OR Conductance:Ext.Temp. Profile (Filename)\tX0\t")
    out.append(f"{temperature}\tNOPE\tNO\tNo\t0\t")
    out.append("# Simulation Length:\t\t\t\t\t")
    out.append("#       DAYS\tHours\tMinutes\tSeconds\tPicoSeconds\t")
    out.append(f"{sim_days}\t{sim_hours}\t{sim_minutes}\t{sim_seconds}\t{sim_pico}\t")
    out.append("#\t                                                     \t\t\t\t")
    out.append("#hv(filename)\tSampling Interval (s)\tPercent(%)\tAccuracy\tX0\t")
    out.append(f"None1\t1\t0\t{accuracy}\t0\t")
    out.append("#  THESE FIELDS BELOW (X0) ARE CURRENTLY NOT USED, LEAVE THEM AT 0\t\t\t\t\t")
    out.append("# CSTR/PFR inlet Flow Temperature\tX0\tX0\tX0\tX0\t")
    out.append("0\t0\t0\t0\t0\t")
    out.append("#\t\t\t\t\t")
    out.append("END\t\t\t\t\t")
    out.append("")  # trailing newline

    output_text = "\n".join(out)

    if output_path:
        write_kintecus_file(output_path, output_text)
        print(f"Written to: {output_path}")
    else:
        print(output_text)

    return output_text




def convert(xml_path: str, output_path: str | None = None,
            species_path: str | None = None, parm_path: str | None = None,
            verbose: bool = False) -> str:
    """
    Parse an SBML XML file and return (or write) a Kintecus model.dat string.
    """
    tree = ET.parse(xml_path)
    root = tree.getroot()

    sbml_ns = detect_sbml_ns(root.tag)
    stag = make_tagger(sbml_ns)

    # ── SBML level / version ──
    sbml_level, sbml_version = get_model_level_version(root)

    # ── L3 package warnings (collected here, reported with other warnings) ──
    pkg_warnings = check_unsupported_packages(root, xml_path)

    model_el = root.find(stag("model"))
    if model_el is None:
        sys.exit("ERROR: No <model> element found in SBML file.")

    model_id   = model_el.get("id",   "UnknownModel")
    model_name = model_el.get("name", model_id)

    # ── Compartments ──
    compartment_ids:    set[str] = set()
    compartment_sizes:  dict[str, str] = {}
    for c in model_el.findall(f".//{stag('compartment')}"):
        cid  = c.get("id")
        if cid:
            compartment_ids.add(cid)
            # L2 uses 'volume', L3 uses 'size'; fall back to 1 if absent
            compartment_sizes[cid] = c.get("size") or c.get("volume") or "1"

    # ── Global parameters (listOfParameters at model level) ──
    global_params: dict[str, str] = {}
    glop = model_el.find(stag("listOfParameters"))
    if glop is not None:
        for p in glop.findall(stag("parameter")):
            pid = p.get("id")
            val = p.get("value", "1.0")
            global_params[pid] = val

    # ── Species (for header info) ──
    species_ids: set[str] = set()
    species_info: dict[str, dict] = {}
    species_using_amount: set[str] = set()   # track species whose IC came from initialAmount
    los = model_el.find(stag("listOfSpecies"))
    if los is not None:
        for sp in los.findall(stag("species")):
            sid = sp.get("id")
            species_ids.add(sid)
            # Prefer initialConcentration; fall back to initialAmount (L3 common)
            ic = sp.get("initialConcentration")
            if ic is None:
                ic = sp.get("initialAmount")
                if ic is not None:
                    species_using_amount.add(sid)
            if ic is None:
                ic = "0"
            species_info[sid] = {
                "name":       sp.get("name", sid),
                "ic":         ic,
                "compartment": sp.get("compartment", ""),
                "boundary":   sp.get("boundaryCondition", "false").lower() == "true",
                "constant":   sp.get("constant",          "false").lower() == "true",
            }

    # ── Multi-compartment species renaming ──
    # When the model has more than one compartment, every species is suffixed
    # with {compartmentID} so Kintecus can distinguish e.g. C{BZ1} from C{BZ2}.
    # With a single compartment the dict stays empty (no rename applied).
    species_rename: dict[str, str] = {}
    if len(compartment_ids) > 1:
        for sid, sdata in species_info.items():
            comp = sdata["compartment"]
            if comp:
                base = sdata["name"] if sdata["name"] else sid
                species_rename[sid] = f"{base}{{{comp}}}"

    # ── Initial assignments (SBML L2v2+ / L3) ──
    # These may override the IC values read from species attributes above.
    ia_values, ia_warnings = parse_initial_assignments(model_el, stag, global_params)
    # Apply only to species; ignore compartment-size or parameter assignments
    ia_applied: list[str] = []
    for sid in species_info:
        if sid in ia_values:
            old_ic = species_info[sid]["ic"]
            species_info[sid]["ic"] = ia_values[sid]
            ia_applied.append(f"{sid}: {old_ic} → {ia_values[sid]}")

    # ── Build output lines ──
    out: list[str] = []
    all_warnings: list[str] = []

    # Prepend any package / IA warnings so they appear first in the warning block
    all_warnings.extend(pkg_warnings)
    all_warnings.extend(ia_warnings)
    if species_using_amount:
        all_warnings.append(
            "The following species used 'initialAmount' (not 'initialConcentration'). "
            "Kintecus expects concentration values — verify or convert manually: "
            + ", ".join(sorted(species_using_amount))
        )

    out.append(f"# Kintecus model.dat — converted from SBML")
    out.append(f"# Model : {model_name}  (id={model_id})")
    out.append(f"# SBML  : Level {sbml_level} Version {sbml_version}")
    out.append(f"# Source: {Path(xml_path).name}")

    # ── Notes ──
    notes_lines = extract_notes(model_el, stag)
    if notes_lines:
        out.append("#")
        out.extend(notes_lines)
    if compartment_ids:
        if len(compartment_ids) == 1:
            out.append(f"# Compartment (ignored in Kintecus): {', '.join(sorted(compartment_ids))}")
        else:
            out.append(f"# Compartments ({len(compartment_ids)} found — species renamed speciesID{{compartmentID}}): "
                       f"{', '.join(sorted(compartment_ids))}")
    if global_params:
        gp_str = "  ".join(f"{k}={v}" for k, v in global_params.items())
        out.append(f"# Global params: {gp_str}")
    if species_info:
        out.append("#")
        out.append("# Species initial concentrations:")
        for sid, sdata in species_info.items():
            display = species_rename.get(sid, sid)  # show renamed form if applicable
            flags = []
            if sdata["boundary"]: flags.append("boundary")
            if sdata["constant"]: flags.append("constant")
            if sid in species_using_amount: flags.append("initialAmount — verify units")
            if sid in ia_values: flags.append("IC from initialAssignment")
            flag_str = f"  [{', '.join(flags)}]" if flags else ""
            out.append(f"#   {display:16s}  IC={sdata['ic']:>12}  ({sdata['name']}){flag_str}")
    out.append("#")
    out.append("# Rate constant  Reaction               Comment")
    out.append("#" + "-" * 70)

    lor = model_el.find(stag("listOfReactions"))
    if lor is None:
        all_warnings.append("No <listOfReactions> found in model.")
    else:
        for rxn in lor.findall(stag("reaction")):
            rxn_id   = rxn.get("id",   "?")
            rxn_name = rxn.get("name", rxn_id)
            reversible = rxn.get("reversible", "true").lower() == "true"

            rxn_warnings: list[str] = []

            # Reactants
            lor_el = rxn.find(stag("listOfReactants"))
            r_refs = lor_el.findall(stag("speciesReference")) if lor_el is not None else []
            lhs, lw = side_string(r_refs, global_params, stag, species_rename, is_reactant=True)
            rxn_warnings.extend(lw)

            # Products
            lop_el = rxn.find(stag("listOfProducts"))
            p_refs = lop_el.findall(stag("speciesReference")) if lop_el is not None else []
            rhs, pw = side_string(p_refs, global_params, stag, species_rename)
            rxn_warnings.extend(pw)

            kl_el = rxn.find(stag("kineticLaw"))

            # Merge global + local params so USER expressions can substitute both
            local_params_raw = _load_local_params(kl_el, stag) if kl_el is not None else {}
            merged_params = {
                **global_params,
                **{k: v["value"] for k, v in local_params_raw.items()},
            }

            if reversible:
                # ── Two lines: forward then reverse ──────────────────────────
                fwd_info, rev_info, kw = get_rate_column(
                    kl_el, r_refs, p_refs, True,
                    species_ids, compartment_ids, merged_params,
                    compartment_sizes, species_rename, global_params, stag,
                )
                rxn_warnings.extend(kw)
                fwd_val, fwd_name = fwd_info
                rev_val, rev_name = rev_info

                fwd_comment = f"# {rxn_name} ({rxn_id}) forward  k={fwd_name}"
                rev_comment = f"# {rxn_name} ({rxn_id}) reverse  k={rev_name}"
                out.append(f"{fwd_val}\t{lhs}==>{rhs}\t{fwd_comment}")
                out.append(f"{rev_val}\t{rhs}==>{lhs}\t{rev_comment}")

            else:
                # ── Single line: irreversible ─────────────────────────────────
                rate_info, _, kw = get_rate_column(
                    kl_el, r_refs, p_refs, False,
                    species_ids, compartment_ids, merged_params,
                    compartment_sizes, species_rename, global_params, stag,
                )
                rxn_warnings.extend(kw)
                rate_val, rate_name = rate_info

                comment = f"# {rxn_name} ({rxn_id})  k={rate_name}"
                out.append(f"{rate_val}\t{lhs}==>{rhs}\t{comment}")

            if rxn_warnings:
                all_warnings.append(f"Reaction {rxn_id} ({rxn_name}):")
                all_warnings.extend(rxn_warnings)

    # ── Non-kinetic rules (assignmentRule / rateRule / algebraicRule) ──
    rule_comments = collect_rule_comments(
        model_el, stag, species_ids, compartment_ids,
        global_params, compartment_sizes, species_rename
    )
    if rule_comments:
        out.extend(rule_comments)
        all_warnings.append(
            "Non-kinetic rules found (assignmentRule/rateRule/algebraicRule). "
            "Kintecus cannot represent these directly — they have been preserved "
            "as comments in MODEL.DAT. Simulation results may differ from the "
            "original SBML model without manual adjustment."
        )

    out.append("END")
    out.append("")  # trailing newline
    output_text = "\n".join(out)

    if output_path:
        write_kintecus_file(output_path, output_text)
        print(f"Written to: {output_path}")
    else:
        print(output_text)

    # ── Species.dat ──
    generate_species_dat(species_info, species_rename, species_path)

    # ── PARM.DAT ──
    generate_parm_dat(global_params, parm_path)

    if all_warnings:
        print("\n--- CONVERSION WARNINGS ---", file=sys.stderr)
        for w in all_warnings:
            print(w, file=sys.stderr)
        if not verbose:
            print(
                "(Use -v / --verbose to always show warnings, or check stderr above)",
                file=sys.stderr
            )

    return output_text


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Convert an SBML Level-2 XML file to a Kintecus model.dat file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python sbml2kintecus.py model.xml                                   # writes MODEL.DAT, SPECIES.DAT, PARM.DAT
  python sbml2kintecus.py model.xml -o model.dat -s species.dat -p PARM.DAT
  python sbml2kintecus.py model.xml -o model.dat -s species.dat -p PARM.DAT -v
""",
    )
    parser.add_argument("input",              help="Input SBML XML file")
    parser.add_argument("-o", "--output",     help="Output model.dat   (default: MODEL.DAT)",   default="MODEL.DAT")
    parser.add_argument("-s", "--species",    help="Output species.dat (default: SPECIES.DAT)", default="SPECIES.DAT")
    parser.add_argument("-p", "--parm",       help="Output PARM.DAT   (default: PARM.DAT)",    default="PARM.DAT")
    parser.add_argument("-v", "--verbose",    action="store_true",
                        help="Always print conversion warnings to stderr")
    args = parser.parse_args()

    if not Path(args.input).is_file():
        sys.exit(f"ERROR: File not found: {args.input}")

    convert(args.input, args.output, args.species, args.parm, args.verbose)


if __name__ == "__main__":
    main()
