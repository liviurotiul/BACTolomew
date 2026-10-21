#!/usr/bin/env python3
"""Upload a database bundle to Zenodo.

Safety design:
  * defaults to the SANDBOX instance (sandbox.zenodo.org), which is throwaway
  * NEVER publishes automatically. Publishing mints a permanent DOI and a
    published record cannot be deleted, so it is a separate explicit command.
  * computes and records SHA-256 for every file
  * supports creating a new VERSION of an existing record

Token: create at
  https://sandbox.zenodo.org/account/settings/applications/tokens/new/   (sandbox)
  https://zenodo.org/account/settings/applications/tokens/new/           (real)
with scopes: deposit:write, deposit:actions
Then:  export ZENODO_TOKEN=...

Usage
-----
  # 1. inspect what would happen
  python zenodo_upload.py create --metadata meta.json --file db.tar.xz --dry-run

  # 2. create a draft on the sandbox and upload
  python zenodo_upload.py create --metadata meta.json --file db.tar.xz

  # 3. when you are certain (add --production for the real Zenodo)
  python zenodo_upload.py publish --id <deposition_id>

  # later: a new version of an existing record
  python zenodo_upload.py newversion --id <old_id> --file db_v2.tar.xz --metadata meta.json
"""
from __future__ import annotations
import argparse, hashlib, json, os, sys
from pathlib import Path

try:
    import requests
except ImportError:
    sys.exit("pip install requests")

SANDBOX = "https://sandbox.zenodo.org/api"
PRODUCTION = "https://zenodo.org/api"


def api_base(args) -> str:
    return PRODUCTION if args.production else SANDBOX


def token() -> str:
    t = os.environ.get("ZENODO_TOKEN")
    if not t:
        sys.exit("Set ZENODO_TOKEN (see the docstring for how to create one).")
    return t


def sha256(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def check(r: requests.Response, what: str):
    if r.status_code >= 300:
        sys.exit(f"{what} failed [{r.status_code}]: {r.text[:800]}")
    return r.json() if r.content else {}


def upload_files(base: str, bucket: str, files: list[Path], tok: str):
    """Upload via the bucket (files) API, which streams and handles large files."""
    for f in files:
        size_mb = f.stat().st_size / 1024 / 1024
        print(f"  uploading {f.name} ({size_mb:,.1f} MB) ...", flush=True)
        with open(f, "rb") as fh:
            r = requests.put(f"{bucket}/{f.name}", data=fh, params={"access_token": tok})
        check(r, f"upload of {f.name}")
        print(f"    sha256 {sha256(f)}")


def cmd_create(args):
    base, tok = api_base(args), token()
    meta = json.loads(Path(args.metadata).read_text())
    files = [Path(f) for f in args.file]
    for f in files:
        if not f.is_file():
            sys.exit(f"missing file: {f}")

    print(f"Target: {base}")
    print(f"Title : {meta.get('metadata', {}).get('title', '?')}")
    print("Files :")
    for f in files:
        print(f"   {f}  ({f.stat().st_size/1024/1024:,.1f} MB)")
    if args.dry_run:
        print("\n[dry run] nothing sent. Metadata that would be posted:")
        print(json.dumps(meta, indent=2)[:2000])
        return

    dep = check(requests.post(f"{base}/deposit/depositions",
                              params={"access_token": tok}, json={}), "create draft")
    dep_id, bucket = dep["id"], dep["links"]["bucket"]
    print(f"\ndraft id {dep_id}")

    upload_files(base, bucket, files, tok)

    check(requests.put(f"{base}/deposit/depositions/{dep_id}",
                       params={"access_token": tok}, json=meta), "set metadata")
    print(f"\nDraft ready (NOT published): {dep['links']['html']}")
    print(f"Review it in the browser, then:  python {Path(__file__).name} publish --id {dep_id}"
          + (" --production" if args.production else ""))


def cmd_newversion(args):
    base, tok = api_base(args), token()
    r = check(requests.post(f"{base}/deposit/depositions/{args.id}/actions/newversion",
                            params={"access_token": tok}), "new version")
    new_id = r["links"]["latest_draft"].rstrip("/").split("/")[-1]
    dep = check(requests.get(f"{base}/deposit/depositions/{new_id}",
                             params={"access_token": tok}), "fetch draft")
    # remove files inherited from the previous version
    for f in dep.get("files", []):
        requests.delete(f"{base}/deposit/depositions/{new_id}/files/{f['id']}",
                        params={"access_token": tok})
    upload_files(base, dep["links"]["bucket"], [Path(f) for f in args.file], tok)
    if args.metadata:
        check(requests.put(f"{base}/deposit/depositions/{new_id}",
                           params={"access_token": tok},
                           json=json.loads(Path(args.metadata).read_text())), "set metadata")
    print(f"New version draft {new_id} ready (NOT published): {dep['links']['html']}")


def cmd_publish(args):
    base, tok = api_base(args), token()
    where = "PRODUCTION Zenodo" if args.production else "sandbox"
    print(f"About to PUBLISH deposition {args.id} on {where}.")
    print("Publishing mints a permanent DOI. A published record CANNOT be deleted.")
    if input("Type the deposition id again to confirm: ").strip() != str(args.id):
        sys.exit("aborted")
    r = check(requests.post(f"{base}/deposit/depositions/{args.id}/actions/publish",
                            params={"access_token": tok}), "publish")
    print(f"Published. DOI: {r.get('doi')}   {r.get('links', {}).get('record_html')}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--production", action="store_true",
                   help="use the real Zenodo instead of the sandbox")
    sub = p.add_subparsers(dest="cmd", required=True)

    c = sub.add_parser("create"); c.set_defaults(fn=cmd_create)
    c.add_argument("--metadata", required=True); c.add_argument("--file", nargs="+", required=True)
    c.add_argument("--dry-run", action="store_true")

    n = sub.add_parser("newversion"); n.set_defaults(fn=cmd_newversion)
    n.add_argument("--id", required=True); n.add_argument("--file", nargs="+", required=True)
    n.add_argument("--metadata")

    b = sub.add_parser("publish"); b.set_defaults(fn=cmd_publish)
    b.add_argument("--id", required=True)

    args = p.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
