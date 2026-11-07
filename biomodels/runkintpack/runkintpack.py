#!/usr/bin/env python3
"""
runkintpack.py -- unpack a Kintecus .kintpack into its text files and run Kintecus.exe on it.
The .kintpack itself is never modified: everything is written into a separate run folder.

    runkintpack.py PACK.kintpack [options]

Options
  -k "SWITCHES"     Kintecus switches; override the switches stored in the pack
  -p PATH           path to Kintecus.exe; overrides the path stored in the pack (settings.kpath)
  -i                international inputs: numbers written with a decimal comma (4,34) are written with a
                    decimal point (4.34) in the text files Kintecus reads (also used if the pack's own
                    "nation" setting is on)
  --set FIELD=VALUE override one field in the unpacked files (repeatable; the pack is not changed):
        species.NAME.COLUMN=value        species.H2O.initialconcentration=4.34
        parm.FIELD=value                 parm.temperature=897      parm.samplinginterval=0.5
        model.LINE="A,m,Ea,reaction,comment"  or  model.LINE="A,reaction,comment"
        model.LINE.FIELD=value           model.12.A=3.4e10    model.12.reaction=A+B==>C
        SHEET.LINE="v1,v2,..." / "#comment"   any other sheet or file: therm, transport, fitdata,
                                         fitdata001, a perturbation/profile file name, ...
        SHEET.LINE.COLUMN=value          COLUMN as a 1-based number (or a header name for sheets with headers)
      LINE is the line number in the written text file (1 = first line), exactly as Kintecus sees it.
      Fields in a VALUE are separated by commas (by semicolons with -i, where the comma is a decimal mark);
      they are written tab-delimited. Lines past the end of a data file (fitdata, profiles) are appended.
  -o DIR            run folder (default: <pack name>_run next to the pack)
  -u, --unpack-only write the files but do not run Kintecus
  --list            show the settable fields (species, columns, PARM fields, model lines, files) and stop
  --check           validate the --set options and stop (nothing written)
  --batch           run unattended: Kintecus -i in its own hidden console, watched through KDONE.OUT /
                    KERR.OUT / CONC.TXT; stops a stuck run (see --stall, --cap). Default: Kintecus runs in this
                    console exactly as from a command prompt.
  --stall SEC       (--batch) stop when CONC.TXT gets no new row for SEC seconds [120]
  --cap SEC         (--batch) stop after SEC seconds [0 = no limit]

Note: with -obeymaxint Kintecus writes CONC.TXT rows every "Maximum Integration Time" when that is larger
than the "Sampling Interval" - set both (parm.maxint=..., parm.sampling=...) to change the output spacing.

Exit code: 0 = Kintecus finished, 1 = Kintecus stopped with a FATAL error / was stopped, 2 = usage or --set error.
"""
import difflib, json, os, re, shutil, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))

# sheet name -> (file name, number of columns or None = as in the sheet, trailing tab)
FILES = {"model": ("MODEL.DAT", 5, False), "species": ("SPECIES.DAT", 9, True), "parm": ("PARM.DAT", 6, True),
         "therm": ("THERM.DAT", None, True), "transport": ("TRANSPORT.DAT", None, True),
         "fitdata": ("FITDATA.TXT", None, False), "spname": ("SPNAME.DAT", None, True),
         "period": ("PERIOD.DAT", None, True)}

SPECIES_COLS = {0: ["species", "name"], 1: ["residence", "residencetime", "cstr"],
                2: ["initialconcentration", "initialconc", "initial", "init", "concentration", "conc"],
                3: ["display", "disp"], 4: ["externalconcentration", "externalconc", "external", "ext"],
                5: ["ssa", "speciesspecial", "special"],
                6: ["constantfile", "constfile", "constant", "const", "constantfileno", "pin"],
                7: ["comments", "comment"]}
MODEL_FIELDS = {"a": 0, "k": 0, "rateconstant": 0, "field1": 0, "m": 1, "tm": 1, "ea": 2, "reaction": 3,
                "equation": 3, "comment": 4, "comments": 4}
NUM_COMMA = re.compile(r"^\s*[+-]?\d*,\d+([eEdD][+-]?\d+)?\s*$|^\s*[+-]?\d+,\d*([eEdD][+-]?\d+)?\s*$")


class SetError(Exception):
    pass


def norm(s):
    s = s.strip().lstrip("#").strip().lower()
    return re.sub(r"[^a-z0-9]", "", s)


def label_key(s):
    """'Sampling Interval (s)' -> 'samplinginterval'; '#Temperature (K) or Filename' -> 'temperature'."""
    s = s.strip().lstrip("#").strip()
    s = re.split(r"\(|\bor\b", s)[0]
    return norm(s)


def intl(cell):
    return cell.replace(",", ".") if NUM_COMMA.match(cell) else cell


def is_num(s):
    try:
        float(s.replace("D", "E").replace("d", "e"))
        return True
    except ValueError:
        return False


# ------------------------------------------------------------------ the unpacked "files" (lists of rows)
class Deck:
    def __init__(self, pack):
        self.files = {}          # key -> dict(name, rows, ncols, trailing, header)
        g = pack.get("grids", {})
        for key, grid in g.items():
            rows = [list(map(lambda c: "" if c is None else str(c), r)) for r in grid.get("rows", [])]
            if not rows:
                continue
            fn, n, tr = FILES.get(key.lower(), (key.upper() + ".DAT", None, True))
            cols = grid.get("cols") or []
            if n is None:
                n = max(len(cols), max(len(r) for r in rows))
            self.files[key.lower()] = {"name": fn, "rows": rows, "ncols": n, "trailing": tr, "cols": cols,
                                       "grid": True}
        for of in pack.get("ofiles", []):
            stem = os.path.splitext(of["filename"])[0].lower()
            rows = [list(map(str, r)) for r in of.get("rows", [])]
            self.files[stem] = {"name": of["filename"], "rows": rows, "ncols": None, "trailing": False,
                                "cols": [], "grid": False}

    def write(self, d, international):
        written = []
        for key, f in self.files.items():
            out = []
            for r in f["rows"]:
                cells = list(r)
                if f["grid"]:
                    n = f["ncols"] if key in ("model", "species", "parm") else max(f["ncols"], len(cells))
                    cells = (cells + [""] * n)[:n]
                    if cells and cells[0] == "":
                        cells[0] = "# "
                if international:
                    cells = [intl(c) for c in cells]
                line = "\t".join(cells)
                out.append(line + ("\t" if f["trailing"] else ""))
            with open(os.path.join(d, f["name"]), "w", newline="", encoding="latin-1", errors="replace") as fh:
                fh.write("\r\n".join(out) + "\r\n")
            written.append(f["name"])
        return written

    # -------------------------------------------------------------- lookups
    def sheet(self, key):
        k = key.lower()
        if k in self.files:
            return k, self.files[k]
        for kk, f in self.files.items():                      # by file name, e.g. FITDATA001.TXT
            if f["name"].lower() == k or os.path.splitext(f["name"])[0].lower() == k:
                return kk, f
        raise SetError(f"no sheet or file '{key}' in this pack; available: {', '.join(sorted(self.files))}")

    def species_row(self, name):
        f = self.files.get("species")
        if not f:
            raise SetError("this pack has no species sheet")
        names = [(i, r[0].strip()) for i, r in enumerate(f["rows"])
                 if r and r[0].strip() and not r[0].lstrip().startswith("#") and r[0].strip().upper() != "END"]
        hit = [i for i, n in names if n == name] or [i for i, n in names if n.upper() == name.upper()]
        if len(hit) == 1:
            return hit[0]
        if len(hit) > 1:
            raise SetError(f"species '{name}' matches more than one row (case-insensitive); use the exact case")
        close = difflib.get_close_matches(name, [n for _, n in names], n=5, cutoff=0.6)
        raise SetError(f"species '{name}' is not in SPECIES.DAT" + (f"; did you mean {', '.join(close)}?" if close else
                                                                   f" ({len(names)} species; see --list)"))

    def parm_fields(self):
        """normalized label -> list of (row, col, label): each '#' description row labels the value row below it."""
        f = self.files.get("parm")
        out = {}
        if not f:
            return out
        rows = f["rows"]
        for i, r in enumerate(rows[:-1]):
            if not (r and r[0].lstrip().startswith("#")):
                continue
            nxt = rows[i + 1]
            if not nxt or nxt[0].lstrip().startswith("#") or nxt[0].strip().upper() == "END":
                continue
            for c, lab in enumerate(r):
                k = label_key(lab)
                if k and c < max(len(nxt), f["ncols"]):
                    out.setdefault(k, []).append((i + 1, c, lab.strip().lstrip("#").strip()))
        return out


ALIASES = {"temp": "temperature", "t": "temperature", "maxtime": "maximumintegrationtime",
           "maxint": "maximumintegrationtime", "maxstep": "maximumintegrationtime",
           "starttime": "startingintegrationtime", "start": "startingintegrationtime",
           "sampling": "samplinginterval", "samp": "samplinginterval", "acc": "accuracy",
           "concunits": "concunits", "units": "concunits", "eaunits": "eaunits", "pressure": "pressure",
           "sec": "seconds", "secs": "seconds", "day": "days", "hour": "hours", "minute": "minutes"}


SEP = [","]          # field separator for --set values: ',' normally, ';' with -i (the comma is then a decimal mark)


def split_fields(value):
    if "\t" in value:
        return value.split("\t")
    return [p.strip() for p in value.split(SEP[0])]


def apply_set(deck, spec, international, log):
    if "=" not in spec:
        raise SetError(f"--set {spec}: expected FIELD=VALUE")
    lhs, value = spec.split("=", 1)
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in "\"'":
        value = value[1:-1]
    num_ok = (lambda v: is_num(intl(v))) if international else is_num
    parts = lhs.strip().split(".")
    sheet = parts[0].lower()
    if not parts[0]:
        raise SetError(f"--set {spec}: missing sheet name")

    if sheet == "species":
        if len(parts) < 3:
            raise SetError(f"--set {spec}: use species.NAME.COLUMN=value (e.g. species.H2O.initialconcentration=4.34)")
        name, col = ".".join(parts[1:-1]), parts[-1]
        try:
            i = deck.species_row(name)
        except SetError as ex:
            raise SetError(f"--set {spec}: {ex}")
        c = None
        if col.isdigit():
            c = int(col) - 1
        else:
            nk = norm(col)
            for idx, names in SPECIES_COLS.items():
                if nk in names:
                    c = idx
            if c is None:
                close = difflib.get_close_matches(nk, [n for v in SPECIES_COLS.values() for n in v], n=1, cutoff=0.7)
                if close:
                    raise SetError(f"--set {spec}: unknown species column '{col}'; did you mean '{close[0]}'?")
                raise SetError(f"--set {spec}: unknown species column '{col}'; use one of "
                               f"{', '.join(n[0] for n in SPECIES_COLS.values())} or a column number 1-8")
        if not 0 <= c < 8:
            raise SetError(f"--set {spec}: species column must be 1-8")
        if c in (1, 2, 4) and not num_ok(value):
            raise SetError(f"--set {spec}: '{value}' is not a number")
        if c in (3, 5) and value.upper() not in ("Y", "N", "YES", "NO"):
            raise SetError(f"--set {spec}: expected Y or N")
        row = deck.files["species"]["rows"][i]
        row += [""] * (c + 1 - len(row))
        old, row[c] = row[c], value
        log.append(f"SPECIES.DAT line {i + 1} ({name}), column {c + 1}: '{old}' -> '{value}'")
        return

    if sheet == "parm":
        fields = deck.parm_fields()
        if len(parts) == 3 and parts[1].isdigit() and parts[2].isdigit():
            i, c = int(parts[1]) - 1, int(parts[2]) - 1
            rows = deck.files["parm"]["rows"]
            if not (0 <= i < len(rows)) or rows[i][0].lstrip().startswith("#"):
                raise SetError(f"--set {spec}: PARM.DAT line {i + 1} is not a value line")
        else:
            if len(parts) != 2:
                raise SetError(f"--set {spec}: use parm.FIELD=value (e.g. parm.temperature=897) or parm.LINE.COLUMN=value")
            k = ALIASES.get(norm(parts[1]), norm(parts[1]))
            hits = fields.get(k) or [h for kk, v in fields.items() if kk.startswith(k) for h in v]
            if not hits:
                close = difflib.get_close_matches(k, list(fields), n=5, cutoff=0.5)
                raise SetError(f"--set {spec}: no PARM field '{parts[1]}'" +
                               (f"; did you mean {', '.join(close)}?" if close else "; see --list"))
            if len({(r, c) for r, c, _ in hits}) > 1:
                raise SetError(f"--set {spec}: '{parts[1]}' matches several PARM fields "
                               f"({'; '.join(f'{lab} (line {r + 1}, column {c + 1})' for r, c, lab in hits)}); "
                               f"use parm.LINE.COLUMN=value")
            i, c, lab = hits[0]
        row = deck.files["parm"]["rows"][i]
        row += [""] * (c + 1 - len(row))
        old, row[c] = row[c], value
        log.append(f"PARM.DAT line {i + 1}, column {c + 1}: '{old}' -> '{value}'")
        return

    key, f = deck.sheet(parts[0])
    if len(parts) < 2 or not parts[1].isdigit():
        raise SetError(f"--set {spec}: use {parts[0]}.LINE=... with LINE the line number in {f['name']}")
    i = int(parts[1]) - 1
    rows = f["rows"]
    if i < 0:
        raise SetError(f"--set {spec}: line numbers start at 1")
    if i >= len(rows):
        if key == "model":
            raise SetError(f"--set {spec}: MODEL.DAT has only {len(rows)} lines")
        rows.extend([[""] for _ in range(i + 1 - len(rows))])        # data files may be extended
    row = rows[i]

    if key == "model":
        if row and row[0].strip().upper() == "END":
            raise SetError(f"--set {spec}: MODEL.DAT line {i + 1} is the END line")
        if len(parts) == 3:
            fk = norm(parts[2])
            if fk not in MODEL_FIELDS:
                raise SetError(f"--set {spec}: model field must be one of A, m, Ea, reaction, comment")
            c = MODEL_FIELDS[fk]
            if c in (0, 1, 2) and not num_ok(value) and not (c == 0 and value.startswith("#")):
                raise SetError(f"--set {spec}: '{value}' is not a number")
            new = (row + [""] * 5)[:5]
            new[c] = value
        elif len(parts) == 2:
            ps = split_fields(value)
            new = (row + [""] * 5)[:5]
            if "\t" in value:
                vals = (ps + [""] * 5)[:5] if len(ps) >= 5 else [ps[0], new[1], new[2]] + (ps[1:] + ["", ""])[:2]
                new = vals
            elif len(ps) >= 5 and all(num_ok(x) or not x.strip() for x in ps[1:3]):
                new = [ps[0], ps[1], ps[2], SEP[0].join(ps[3:-1]), ps[-1]]
            elif len(ps) >= 3:
                new = [ps[0], new[1], new[2], SEP[0].join(ps[1:-1]), ps[-1]]
            elif len(ps) == 2:
                new = [ps[0], new[1], new[2], ps[1], new[4]]
            else:
                raise SetError(f'--set {spec}: use model.LINE="A,m,Ea,reaction,comment" or "A,reaction,comment"')
        else:
            raise SetError(f"--set {spec}: use model.LINE=... or model.LINE.FIELD=...")
        a0 = new[0].strip()
        if not (a0.startswith("#") or num_ok(a0)):
            raise SetError(f"--set {spec}: first field '{a0}' is not a number (or # for a comment line)")
        if not a0.startswith("#") and "=" not in new[3]:
            raise SetError(f"--set {spec}: '{new[3]}' does not look like a reaction (no '=' or '==>')")
        if not a0.startswith("#"):
            for c in (1, 2):
                if new[c].strip() and not num_ok(new[c]):
                    raise SetError(f"--set {spec}: field {c + 1} '{new[c]}' is not a number")
            known = None
            if "species" in deck.files:
                known = {r[0].strip().upper() for r in deck.files["species"]["rows"] if r and r[0].strip()}
            if known:
                rx = re.sub(r"\+?\s*[A-Z]\[[^\]]*\]", "", new[3], flags=re.I)
                for s in re.split(r"==>|=|\+", rx):
                    s = re.sub(r"^[0-9.]+(?=[A-Za-z])", "", s.strip())
                    if s and s.upper() not in known and s.upper() not in ("NULL", "0 NULL", "M", "HV"):
                        log.append(f"warning: MODEL.DAT line {i + 1}: '{s}' is not in SPECIES.DAT "
                                   f"(Kintecus can create it with -c)")
        old = "\t".join(row)
        rows[i] = new
        if row and row[0].lstrip().startswith("#") and not new[0].startswith("#"):
            log.append(f"note: MODEL.DAT line {i + 1} was a comment line and is now a reaction")
        log.append(f"MODEL.DAT line {i + 1}: '{old}' -> '{chr(9).join(new)}'")
        return

    # any other sheet / file
    if len(parts) == 3:
        col = parts[2]
        if col.isdigit():
            c = int(col) - 1
        else:
            heads = [norm(h) for h in f["cols"]]
            if norm(col) not in heads:
                raise SetError(f"--set {spec}: column '{col}' not found in {f['name']}; use a column number")
            c = heads.index(norm(col))
        row += [""] * (c + 1 - len(row))
        old = row[c]; row[c] = value
        log.append(f"{f['name']} line {i + 1}, column {c + 1}: '{old}' -> '{value}'")
    elif len(parts) == 2:
        old = "\t".join(row)
        rows[i] = [value] if value.startswith("#") else split_fields(value)
        log.append(f"{f['name']} line {i + 1}: '{old}' -> '{chr(9).join(rows[i])}'")
    else:
        raise SetError(f"--set {spec}: use {parts[0]}.LINE=... or {parts[0]}.LINE.COLUMN=...")


def show_fields(deck):
    print("Settable fields (use with --set):")
    if "species" in deck.files:
        rows = deck.files["species"]["rows"]
        names = [r[0].strip() for r in rows if r and r[0].strip() and not r[0].startswith("#") and r[0].strip() != "END"]
        print(f"\n species.NAME.COLUMN   ({len(names)} species)")
        print("   columns: " + ", ".join(f"{n[0]} ({i + 1})" for i, n in SPECIES_COLS.items()))
        print("   species: " + ", ".join(names[:60]) + (" ..." if len(names) > 60 else ""))
    pf = deck.parm_fields()
    if pf:
        print("\n parm.FIELD   (or parm.LINE.COLUMN)")
        for k, hits in pf.items():
            if k == "x0":
                continue
            for r, c, lab in hits:
                v = deck.files["parm"]["rows"][r]
                print(f"   {k:32s} line {r + 1:2d} col {c + 1}  {lab[:40]:40s} = {v[c] if c < len(v) else ''}")
    for key, f in deck.files.items():
        if key in ("species", "parm"):
            continue
        print(f"\n {key}.LINE  ->  {f['name']} ({len(f['rows'])} lines)")
        if key == "model":
            for i, r in enumerate(f["rows"]):
                if r and not r[0].lstrip().startswith("#") and r[0].strip() != "END":
                    print(f"   line {i + 1:4d}: {r[0]}\t{r[1] if len(r) > 1 else ''}\t{r[2] if len(r) > 2 else ''}\t"
                          f"{(r[3] if len(r) > 3 else '')[:90]}")


# ------------------------------------------------------------------ running
def find_exe(opt, pack):
    cands = [opt] if opt else []
    kp = (pack.get("settings") or {}).get("kpath") or ""
    if kp:
        cands += [kp, os.path.join(kp, "Kintecus.exe"), os.path.join(kp, "kintecus.exe")]
    cands += [os.environ.get("KINTECUS_EXE", ""), os.path.join(HERE, "Kintecus.exe"), os.path.join(os.getcwd(), "Kintecus.exe")]
    for c in cands:
        if c and os.path.isfile(c):
            return os.path.abspath(c)
    w = shutil.which("Kintecus.exe") or shutil.which("kintecus")
    if w:
        return w
    return None


def conc_rows(d):
    try:
        with open(os.path.join(d, "CONC.TXT"), "rb") as fh:
            lines = [l for l in fh.read().split(b"\n") if l.strip()]
    except OSError:
        return 0, 0.0
    for l in reversed(lines[1:]):
        try:
            return len(lines) - 1, float(l.split(b"\t")[0])
        except ValueError:
            pass
    return max(len(lines) - 1, 0), 0.0


def run(exe, switches, d, batch, stall, cap, runner):
    for f in ("KDONE.OUT", "KERR.OUT", "KRUN.TMP"):
        try:
            os.remove(os.path.join(d, f))
        except OSError:
            pass
    args = [exe] + switches.split()
    if not batch:
        cmd = (runner.replace("{cmd}", " ".join(args)) if runner else args)
        rc = subprocess.call(cmd, cwd=d, shell=bool(runner))
        return ("fatal" if os.path.exists(os.path.join(d, "KERR.OUT")) else "done"), f"exit code {rc}"
    if "-i" not in [s.lower() for s in switches.split()]:
        args.append("-i")
    kw = dict(cwd=d, stdin=subprocess.DEVNULL)
    if runner:
        target = runner.replace("{cmd}", " ".join(args)); kw.update(shell=True, stdout=subprocess.DEVNULL,
                                                                    stderr=subprocess.DEVNULL)
        if os.name != "nt":
            kw["start_new_session"] = True
    else:
        target = args
        if os.name == "nt":
            kw["creationflags"] = subprocess.CREATE_NEW_CONSOLE
            si = subprocess.STARTUPINFO(); si.dwFlags |= subprocess.STARTF_USESHOWWINDOW; si.wShowWindow = 0
            kw["startupinfo"] = si
        else:
            kw.update(stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    p = subprocess.Popen(target, **kw)
    t0 = time.time(); last, last_t, state, msg = 0, t0, None, ""
    while True:
        time.sleep(0.5)
        now = time.time()
        if os.path.exists(os.path.join(d, "KERR.OUT")):
            time.sleep(1); state = "fatal"
            msg = " ".join(open(os.path.join(d, "KERR.OUT"), encoding="latin-1", errors="replace").read().split())
            break
        if os.path.exists(os.path.join(d, "KDONE.OUT")):
            time.sleep(0.5)
            txt = " ".join(open(os.path.join(d, "KDONE.OUT"), encoding="latin-1", errors="replace").read().split())
            state, msg = ("fatal", txt) if "FATAL" in txt.upper() else ("done", txt)
            break
        if p.poll() is not None:
            time.sleep(0.5)
            if not any(os.path.exists(os.path.join(d, f)) for f in ("KDONE.OUT", "KERR.OUT")):
                state, msg = "done", f"Kintecus exited (code {p.returncode})"
                break
            continue
        n, ts = conc_rows(d)
        if n > last:
            last, last_t = n, now
        if stall and now - last_t > stall and now - t0 > stall:
            state, msg = "stopped", f"no new CONC.TXT row for {stall:g} s (simulated time {ts:g})"
        elif cap and now - t0 > cap:
            state, msg = "stopped", f"over the {cap:g} s limit (simulated time {ts:g})"
        if state:
            break
    if p.poll() is None:
        try:
            if runner and os.name != "nt":
                import signal
                os.killpg(p.pid, signal.SIGKILL)
            else:
                p.kill()
        except Exception:
            pass
    if runner and os.name != "nt" and os.path.isdir("/proc"):   # Wine: Kintecus.exe outlives its launcher
        real = os.path.realpath(d)
        for pid in os.listdir("/proc"):
            try:
                if pid.isdigit() and open(f"/proc/{pid}/comm").read().strip() == "Kintecus.exe" and \
                        os.path.realpath(os.readlink(f"/proc/{pid}/cwd")) == real:
                    os.kill(int(pid), 9)
            except OSError:
                pass
    return state, msg


def parse_args(argv):
    """Hand-rolled so that -k "-INT:5 -obeymaxint" (a value starting with '-') works."""
    a = {"pack": None, "k": None, "p": None, "i": False, "sets": [], "o": None, "unpack": False, "list": False,
         "check": False, "batch": False, "stall": 120.0, "cap": 0.0, "runner": None}
    i = 0
    while i < len(argv):
        x = argv[i]
        def val():
            nonlocal i
            if "=" in x and x.startswith("--") and x.split("=", 1)[0] in ("--switches", "--set", "--kintecus", "--out",
                                                                           "--stall", "--cap", "--runner"):
                return x.split("=", 1)[1]
            if i + 1 >= len(argv):
                raise SetError(f"{x} needs a value")
            i += 1
            return argv[i]
        if x in ("-h", "--help", "/?"):
            print(__doc__); sys.exit(0)
        elif x in ("-k", "--switches") or x.startswith("--switches="):
            a["k"] = val()
        elif x in ("-p", "--kintecus") or x.startswith("--kintecus="):
            a["p"] = val()
        elif x in ("-i", "--international"):
            a["i"] = True
        elif x == "--set" or x.startswith("--set="):
            a["sets"].append(val())
        elif x in ("-o", "--out") or x.startswith("--out="):
            a["o"] = val()
        elif x in ("-u", "--unpack-only"):
            a["unpack"] = True
        elif x == "--list":
            a["list"] = True
        elif x == "--check":
            a["check"] = True
        elif x == "--batch":
            a["batch"] = True
        elif x == "--stall" or x.startswith("--stall="):
            a["stall"] = float(val())
        elif x == "--cap" or x.startswith("--cap="):
            a["cap"] = float(val())
        elif x == "--runner" or x.startswith("--runner="):
            a["runner"] = val()
        elif x.startswith("-") and a["pack"] is not None:
            raise SetError(f"unknown option {x} (see --help)")
        elif a["pack"] is None:
            a["pack"] = x
        else:
            raise SetError(f"unexpected argument {x}")
        i += 1
    if not a["pack"]:
        raise SetError("no .kintpack given (see --help)")
    return a


def main(argv=None):
    try:
        a = parse_args(sys.argv[1:] if argv is None else argv)
    except SetError as ex:
        print(f"runkintpack: {ex}", file=sys.stderr); return 2
    pk = os.path.abspath(a["pack"])
    try:
        with open(pk, encoding="utf-8") as fh:
            pack = json.load(fh)
    except (OSError, ValueError) as ex:
        print(f"runkintpack: cannot read {pk}: {ex}", file=sys.stderr); return 2
    if not isinstance(pack, dict) or "grids" not in pack:
        print(f"runkintpack: {pk} is not a .kintpack (no grids)", file=sys.stderr); return 2
    deck = Deck(pack)
    if a["list"]:
        show_fields(deck); return 0
    settings = pack.get("settings") or {}
    international = a["i"] or bool(settings.get("nation"))
    if international:
        SEP[0] = ";"
    log, errors = [], []
    for s in a["sets"]:
        try:
            apply_set(deck, s, international, log)
        except SetError as ex:
            errors.append(str(ex))
    for l in log:
        print("  " + l)
    if errors:
        for e in errors:
            print(f"runkintpack: {e}", file=sys.stderr)
        return 2
    if a["check"]:
        print(f"{len(a['sets'])} --set option(s) OK"); return 0
    switches = a["k"] if a["k"] is not None else (settings.get("switches") or "")
    out = os.path.abspath(a["o"] or os.path.splitext(pk)[0] + "_run")
    if os.path.abspath(out) == os.path.dirname(pk) and os.path.exists(os.path.join(out, "MODEL.DAT")):
        print("note: writing into the pack's own folder; existing DAT files there are overwritten")
    os.makedirs(out, exist_ok=True)
    written = deck.write(out, international)
    with open(os.path.join(out, "runkintpack_log.txt"), "w") as fh:
        fh.write(f"pack: {pk}\nswitches: {switches}\ninternational: {international}\n")
        fh.write("files: " + ", ".join(written) + "\n")
        fh.write("--set:\n" + "".join(f"  {l}\n" for l in log))
    print(f"unpacked {os.path.basename(pk)} -> {out}\n  {', '.join(written)}")
    if a["unpack"]:
        return 0
    exe = find_exe(a["p"], pack)
    if not exe:
        print("runkintpack: Kintecus.exe not found - use -p PATH (or set KINTECUS_EXE)", file=sys.stderr); return 2
    print(f"running {exe} {switches}" + ("   (batch mode)" if a["batch"] else "") + f"\n  in {out}")
    t0 = time.time()
    state, msg = run(exe, switches, out, a["batch"], a["stall"], a["cap"], a["runner"])
    n, ts = conc_rows(out)
    print(f"Kintecus {'finished' if state == 'done' else 'stopped' if state == 'stopped' else 'FATAL'}: {msg}  "
          f"({time.time() - t0:.1f} s; CONC.TXT {n} rows, last time {ts:g})")
    outs = sorted(f for f in os.listdir(out) if f.upper().endswith((".TXT", ".OUT")) and f not in written
                  and f != "runkintpack_log.txt")
    if outs:
        print("  output files: " + ", ".join(outs))
    return 0 if state == "done" else 1


if __name__ == "__main__":
    sys.exit(main())
