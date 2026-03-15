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
                species_rename: dict) -> tuple[str, list[str]]:
    """
    Build a '+'-joined species string (with repetition for stoichiometry > 1)
    from a list of speciesReference elements.
    species_rename maps original species IDs to their Kintecus display name
    (e.g. "C" -> "C{BZ1}" in multi-compartment models; identity otherwise).
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
    return "+".join(parts) if parts else "NULL", warnings


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

    model_el = root.find(stag("model"))
    if model_el is None:
        sys.exit("ERROR: No <model> element found in SBML file.")

    model_id   = model_el.get("id",   "UnknownModel")
    model_name = model_el.get("name", model_id)

    # ── Compartments ──
    compartment_ids = set()
    for c in model_el.findall(f".//{stag('compartment')}"):
        compartment_ids.add(c.get("id"))

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
    los = model_el.find(stag("listOfSpecies"))
    if los is not None:
        for sp in los.findall(stag("species")):
            sid = sp.get("id")
            species_ids.add(sid)
            species_info[sid] = {
                "name":       sp.get("name", sid),
                "ic":         sp.get("initialConcentration", "0"),
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

    # ── Build output lines ──
    out: list[str] = []
    all_warnings: list[str] = []

    out.append(f"# Kintecus model.dat — converted from SBML")
    out.append(f"# Model : {model_name}  (id={model_id})")
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
            lhs, lw = side_string(r_refs, global_params, stag, species_rename)
            rxn_warnings.extend(lw)

            # Products
            lop_el = rxn.find(stag("listOfProducts"))
            p_refs = lop_el.findall(stag("speciesReference")) if lop_el is not None else []
            rhs, pw = side_string(p_refs, global_params, stag, species_rename)
            rxn_warnings.extend(pw)

            kl_el = rxn.find(stag("kineticLaw"))

            if reversible:
                # ── Two lines: forward then reverse ──────────────────────────
                if kl_el is not None:
                    fwd_val, fwd_name, rev_val, rev_name, kw = extract_reversible_rates(
                        kl_el, species_ids, compartment_ids, stag, global_params
                    )
                    rxn_warnings.extend(kw)
                else:
                    fwd_val, fwd_name = "1.0", "?"
                    rev_val, rev_name = "1.0", "?"
                    rxn_warnings.append("No <kineticLaw> element; both rates set to 1.0")

                fwd_comment = f"# {rxn_name} ({rxn_id}) forward  k={fwd_name}"
                rev_comment = f"# {rxn_name} ({rxn_id}) reverse  k={rev_name}"
                out.append(f"{fwd_val}\t{lhs}==>{rhs}\t{fwd_comment}")
                out.append(f"{rev_val}\t{rhs}==>{lhs}\t{rev_comment}")

            else:
                # ── Single line: irreversible ─────────────────────────────────
                if kl_el is not None:
                    rate_val, rate_name, kw = extract_irreversible_rate(
                        kl_el, species_ids, compartment_ids, stag, global_params
                    )
                    rxn_warnings.extend(kw)
                else:
                    rate_val, rate_name = "1.0", "?"
                    rxn_warnings.append("No <kineticLaw> element; rate set to 1.0")

                comment = f"# {rxn_name} ({rxn_id})  k={rate_name}"
                out.append(f"{rate_val}\t{lhs}==>{rhs}\t{comment}")

            if rxn_warnings:
                all_warnings.append(f"Reaction {rxn_id} ({rxn_name}):")
                all_warnings.extend(rxn_warnings)

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
