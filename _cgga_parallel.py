"""Run preprocess.py in parallel across N worker processes.

Splits the input manifest into N shards (stride-sharded), spawns one
preprocess.py subprocess per shard. brainles_preprocessing uses per-call
tempfile dirs so concurrent runs don't collide.

Usage:  python _cgga_parallel.py <manifest> <output> <n_workers>
"""
from __future__ import annotations
import json, subprocess, sys, time
from pathlib import Path


def main():
    manifest = Path(sys.argv[1])
    output   = Path(sys.argv[2])
    workers  = int(sys.argv[3])

    subjects = json.loads(manifest.read_text())
    shard_dir = manifest.parent / "_shards"
    shard_dir.mkdir(exist_ok=True)
    shard_paths = []
    for r in range(workers):
        shard = subjects[r::workers]
        p = shard_dir / f"shard_{r:02d}.json"
        p.write_text(json.dumps(shard, indent=2))
        shard_paths.append(p)
        print(f"shard {r}: {len(shard)} subjects")

    log_dir = output / "_parallel_logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    procs = []
    for r, sp in enumerate(shard_paths):
        log_path = log_dir / f"shard_{r:02d}.log"
        log_file = open(log_path, "w")
        proc = subprocess.Popen(
            [sys.executable, "preprocess.py", "--manifest", str(sp),
             "--output", str(output), "--device", "cpu"],
            stdout=log_file, stderr=subprocess.STDOUT,
        )
        procs.append((proc, log_file, log_path))
        print(f"launched shard {r}: pid={proc.pid}  log={log_path}")
    print(f"\nrunning {workers} workers in parallel — see {log_dir} for progress\n")

    t0 = time.time()
    while any(p.poll() is None for p, _, _ in procs):
        alive = sum(1 for p, _, _ in procs if p.poll() is None)
        # count finished subjects across all shard outputs
        done = len([f for f in output.glob("*_t1_preprocessed.nii.gz")])
        print(f"[{int(time.time()-t0):>5d}s]  alive workers: {alive}/{workers}  "
              f"finished subjects (any T1): {done}/{len(subjects)}", flush=True)
        time.sleep(60)
    elapsed = time.time() - t0

    for proc, log_file, _ in procs:
        log_file.close()
    rcs = [p.returncode for p, _, _ in procs]
    print(f"\nall workers done in {elapsed:.0f}s. exit codes: {rcs}")


if __name__ == "__main__":
    main()
