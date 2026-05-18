"""Parallel HTTP-Range downloader. 8 chunks, infinite retries."""
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
import urllib.request
import urllib.error

URL = sys.argv[1]
OUT = sys.argv[2]
EXPECTED = int(sys.argv[3])
N_PARTS = 8


def head_size(url):
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req, timeout=30) as r:
        return int(r.headers["Content-Length"])


def download_range(url, start, end, fout, retry_pause=10):
    while True:
        try:
            req = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}"})
            with urllib.request.urlopen(req, timeout=120) as r:
                pos = start
                while True:
                    chunk = r.read(1 << 16)
                    if not chunk:
                        break
                    with open(fout, "r+b") as fh:
                        fh.seek(pos)
                        fh.write(chunk)
                    pos += len(chunk)
                if pos - 1 == end:
                    return
                # short read -> retry remainder
                start = pos
        except (urllib.error.URLError, ConnectionError, TimeoutError, OSError) as e:
            print(f"part {start}-{end} retry after {e}", flush=True)
            time.sleep(retry_pause)


def main():
    size = head_size(URL)
    if size != EXPECTED:
        print(f"WARNING: expected {EXPECTED}, server says {size}; using server", flush=True)
        EXPECTED_LOCAL = size
    else:
        EXPECTED_LOCAL = EXPECTED
    if not os.path.exists(OUT) or os.path.getsize(OUT) != EXPECTED_LOCAL:
        with open(OUT, "wb") as fh:
            fh.truncate(EXPECTED_LOCAL)
    chunk = EXPECTED_LOCAL // N_PARTS
    parts = []
    for i in range(N_PARTS):
        s = i * chunk
        e = (i + 1) * chunk - 1 if i < N_PARTS - 1 else EXPECTED_LOCAL - 1
        parts.append((s, e))
    print(f"downloading {EXPECTED_LOCAL} bytes in {N_PARTS} parts", flush=True)
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=N_PARTS) as ex:
        futs = [ex.submit(download_range, URL, s, e, OUT) for s, e in parts]
        for f in as_completed(futs):
            f.result()
    print(f"done in {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
