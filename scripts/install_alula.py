"""Explicit, pinned optional installation; never called by a solver execution."""

import argparse
import hashlib
import json
from pathlib import Path
from urllib.request import Request, urlopen
import zipfile

ROOT = Path(__file__).resolve().parents[1]
LOCK = json.loads((ROOT/"src/aeroopt/solvers/_alula/upstream-lock.json").read_text(encoding="utf-8"))


def digest(data):
    return hashlib.sha256(data).hexdigest()


def obtain(url, target, expected):
    if target.exists():
        data = target.read_bytes()
    else:
        with urlopen(Request(url,headers={"User-Agent":"AeroOpt-pinned-install"}),timeout=60) as response:
            data = response.read(100_000_001)
        if len(data) > 100_000_000 or digest(data) != expected:
            raise ValueError(f"Download size/hash mismatch: {url}")
        target.parent.mkdir(parents=True,exist_ok=True)
        with target.open("xb") as stream:
            stream.write(data)
    if digest(data) != expected:
        raise ValueError(f"Source hash mismatch: {target}")
    return data


def install(source=None, *, download=False, portable_node=False):
    if (source is None) == (not download):
        raise ValueError("Choose exactly one of --source-dir or --download")
    target = ROOT/".deps/alula"/LOCK["commit"]
    originals = ROOT/"outputs/dependency_sources_v040"/f"alula-{LOCK['commit']}"
    for name,expected in LOCK["files"].items():
        if source is None:
            data = obtain(f"https://raw.githubusercontent.com/kennyjensen/alula/{LOCK['commit']}/{name}",
                          originals/name,expected)
        else:
            data = (source/name).read_bytes()
            if digest(data) != expected:
                raise ValueError(f"Pinned source mismatch: {name}")
        dest = target/name
        if dest.exists():
            if dest.read_bytes() != data:
                raise ValueError(f"Refusing to overwrite changed source: {dest}")
        else:
            dest.parent.mkdir(parents=True,exist_ok=True)
            with dest.open("xb") as stream:
                stream.write(data)
    (target/"AEROOPT_UPSTREAM.json").write_text(json.dumps(LOCK,indent=2)+"\n",encoding="utf-8")
    print(f"Alula: {target}")
    if portable_node:
        if not download:
            raise ValueError("--portable-node requires explicit --download")
        node = LOCK["tested_node"]
        archive = ROOT/"outputs/dependency_sources_v040"/node["url"].rsplit("/",1)[1]
        obtain(node["url"],archive,node["archive_sha256"])
        prefix = archive.stem
        dest = ROOT/".deps"/prefix
        with zipfile.ZipFile(archive) as bundle:
            for name in ("node.exe","LICENSE","README.md","CHANGELOG.md"):
                data = bundle.read(f"{prefix}/{name}")
                path = dest/name
                if path.exists() and path.read_bytes() != data:
                    raise ValueError(f"Refusing to overwrite changed Node file: {path}")
                if not path.exists():
                    path.parent.mkdir(parents=True,exist_ok=True)
                    with path.open("xb") as stream:
                        stream.write(data)
        if digest((dest/"node.exe").read_bytes()) != node["node_exe_sha256"]:
            raise ValueError("Node executable hash mismatch")
        print(f"Node: {dest/'node.exe'}")
    return target


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir",type=Path,help="Offline upstream checkout/snapshot; every selected file is pinned")
    parser.add_argument("--download",action="store_true",help="Explicit network installation of fixed sources")
    parser.add_argument("--portable-node",action="store_true",help="Also install the pinned Windows x64 Node ZIP locally")
    args = parser.parse_args()
    install(args.source_dir,download=args.download,portable_node=args.portable_node)
