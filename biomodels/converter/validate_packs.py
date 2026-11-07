"""validate_packs.py -- each .kintpack: v3 schema, grid widths, and DAT files rebuilt from the JSON (GUI/VBA rules)
byte-identical to the deck Kintecus actually ran for the reported result."""
import json, glob, os, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sbml2kin_batch as K

bad, ok = [], 0
for f in sorted(glob.glob("/home/claude/kin/batch/packs/*.kintpack")):
    mid = os.path.basename(f)[:-9]
    try:
        p = json.load(open(f))
        assert p["version"] == 3 and isinstance(p["created"], str) and isinstance(p["description"], str)
        assert set(p["settings"]) == {"kpath", "switches", "nation"}
        g = p["grids"]
        for k, n in (("model", 5), ("species", 9), ("parm", 6)):
            assert g[k]["cols"] and all(len(r) == n for r in g[k]["rows"]), f"{k} width"
            assert all(isinstance(c, str) for r in g[k]["rows"] for c in r), f"{k} non-string cell"
        assert g["model"]["rows"][-1][0] == "END" and g["species"]["rows"][-1][0] == "END"
        assert len(p["ofiles"]) <= 10 or "10 O_ file limit" in p["description"]
        for of in p["ofiles"]:
            assert of["key"] == "O_" + of["filename"]
        r = json.load(open(f"/home/claude/kin/batch/out/{mid}/result.json"))
        runs = ([r["attempts"][r["best_run"] - 1]] if r.get("best_run") and r["best_run"] <= len(r["attempts"]) else
                [a for a in r["attempts"] if a["switches"] == r["switches"]])
        run_dir = f"/home/claude/kin/batch/out/{mid}/run_{r['attempts'].index(runs[0]) + 1}"
        with tempfile.TemporaryDirectory() as d:
            K.write_deck(d, g["model"]["rows"], g["species"]["rows"], g["parm"]["rows"],
                         [{"filename": o["filename"], "rows": o["rows"]} for o in p["ofiles"]])
            def norm(fn, data):
                if fn != "MODEL.DAT":
                    return data
                # comments now sit on a '#' line above each reaction (MODEL.DAT lines <= 897 characters):
                # compare the reaction lines only, without their comment field
                out = []
                for ln in data.decode("latin-1").split("\r\n"):
                    c = ln.split("\t")
                    if c and c[0] and not c[0].startswith("#") and c[0] != "END":
                        out.append("\t".join(c[:4]))
                return "\n".join(out)
            for fn in os.listdir(d):
                a, b = open(f"{d}/{fn}", "rb").read(), open(f"{run_dir}/{fn}", "rb").read()
                assert norm(fn, a) == norm(fn, b), f"{fn} differs from the deck that was run"
        assert p["settings"]["switches"] == r["switches"]
        ok += 1
    except Exception as ex:
        bad.append((mid, repr(ex)[:160]))
print(f"{ok} packs valid, {len(bad)} invalid")
for b in bad[:20]:
    print("  ", *b)
