#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0
"""Verify published wheel hashes and install from PyPI in a fresh environment."""
import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import urlopen
import venv

from build_python_wheel import isolated_environment


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--fixture", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output = args.output.resolve()
    if args.output.exists():
        raise RuntimeError("Refusing to overwrite an existing index-install test")
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    version = manifest["version"]
    if manifest.get("schema") != "paacs-python-release-v1" or len(manifest.get("wheels", {})) != 3:
        raise RuntimeError("Expected the qualified three-platform publication manifest")
    if not re.fullmatch(r"[A-Za-z0-9.!+_-]+", version):
        raise RuntimeError("Invalid release version")
    for attempt in range(6):
        try:
            with urlopen("https://pypi.org/pypi/paacs-mesh/" + quote(version, safe="") + "/json", timeout=15) as response:
                index = json.load(response)
            actual = {record["filename"]: record["digests"]["sha256"] for record in index["urls"]
                      if record["packagetype"] == "bdist_wheel"}
            if actual == manifest["wheels"]:
                break
        except (HTTPError, URLError):
            if attempt == 5:
                raise
        if attempt == 5:
            raise RuntimeError("PyPI does not expose the exact three qualified wheels")
        time.sleep(5)
    args.output.mkdir(parents=True)
    environment = args.output / "environment"
    venv.EnvBuilder(with_pip=True).create(environment)
    python = environment / ("Scripts/python.exe" if sys.platform == "win32" else "bin/python")
    clean = isolated_environment(environment)
    requirement = environment / "release.txt"
    requirement.write_text("paacs-mesh==" + version + " " + " ".join(
        "--hash=sha256:" + digest for digest in manifest["wheels"].values()) + "\n", encoding="utf-8")
    subprocess.run([str(python), "-m", "pip", "--isolated", "install", "--index-url", "https://pypi.org/simple",
        "--only-binary=:all:", "--no-cache-dir", "--no-deps", "--require-hashes", "-r", str(requirement)], check=True,
        cwd=environment, env=clean)
    fixture = environment / "cylinder.step"
    shutil.copyfile(args.fixture, fixture)
    code = """import json, paacs_mesh as pm
model = pm.read_step('cylinder.step')
mesh = model.mesh(.3)
r = mesh.report
assert r['audit']['complete'] and r['audit']['valid'] and r['quality']['valid']
assert r['quality']['inverted_cells'] == 0 and mesh.body_tags == [1]
assert mesh.face_tags == [1, 2, 3]
assert r['bodies'][0]['initial_front']['nodes'] > 0
mesh.write('cylinder.paacsmesh')
mesh.write('cylinder.vtk')
print(json.dumps({'version': pm.__version__, 'nodes': mesh.node_count, 'tet4': mesh.tetrahedron_count, 'audit': r['audit']}))
"""
    subprocess.run([str(python), "-I", "-c", code], check=True, cwd=environment, env=clean)
    console = environment / ("Scripts/paacs-mesh.exe" if sys.platform == "win32" else "bin/paacs-mesh")
    subprocess.run([str(console), "--help"], check=True, cwd=environment, env=clean)
    receipt = {"version": version, "source_commit": manifest["source_commit"],
               "published_wheels": actual, "fresh_index_install": "passed"}
    (args.output / "index-validation.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt))


if __name__ == "__main__":
    main()
