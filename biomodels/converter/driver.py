"""driver.py IDFILE [workers] -- run run_model.py over ids with a per-model wall limit; skips ids already done."""
import subprocess, sys, os, json, time
from concurrent.futures import ThreadPoolExecutor
ids = [l.strip() for l in open(sys.argv[1]) if l.strip()]
W = int(sys.argv[2]) if len(sys.argv) > 2 else 2
LIMIT = int(os.environ.get("MODEL_LIMIT", "700"))
def one(mid):
    OUT = os.environ.get("KIN_OUT", "/home/claude/kin/batch/out")
    rj = f"{OUT}/{mid}/result.json"
    if os.environ.get("SKIP_DONE") and os.path.exists(rj):
        return
    t = time.time()
    try:
        subprocess.run(["python3", "/home/claude/kin/batch/run_model.py", mid], timeout=LIMIT,
                       capture_output=True, start_new_session=True)
    except subprocess.TimeoutExpired:
        os.makedirs(os.path.dirname(rj), exist_ok=True)
        json.dump({"id": mid, "status": "TIMEOUT_TOTAL", "seconds": LIMIT, "flags": [], "attempts": []}, open(rj, "w"))
        base = os.path.realpath(f"{OUT}/{mid}")      # kill only this model's Kintecus
        for pid in os.listdir("/proc"):
            try:
                if pid.isdigit() and open(f"/proc/{pid}/comm").read().strip() == "Kintecus.exe" and \
                        os.path.realpath(os.readlink(f"/proc/{pid}/cwd")).startswith(base + "/"):
                    os.kill(int(pid), 9)
            except OSError:
                pass
    try:
        r = json.load(open(rj)); print(f"{mid} {r['status']} err={r.get('error')} {round(time.time()-t)}s", flush=True)
    except Exception as ex:
        print(f"{mid} NO_RESULT {ex}", flush=True)
with ThreadPoolExecutor(W) as ex:
    list(ex.map(one, ids))
