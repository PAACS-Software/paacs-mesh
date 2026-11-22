# Python package releases

The first package-index distribution targets CPython 3.12 on Windows x64,
macOS 15 ARM64 and Linux x64 with glibc 2.28 or newer. Metadata restricts Python
to `>=3.12,<3.13` until other interpreter versions have qualified wheels.
Runtime CAD libraries are bundled. Only PAACS wheels enter PyPI; the matching
OCCT source archive and build recipes remain available with the GitHub release.

## One-time account setup

The configured PyPI Trusted Publisher uses the following identity. A new
maintainer can reproduce this mapping at the account's **Publishing** page:

| Setting | Value |
|---|---|
| PyPI project name | `paacs-mesh` |
| GitHub repository owner | `PAACS-Software` |
| GitHub repository | `paacs-mesh` |
| Workflow filename | `publish-python.yml` |
| Environment name | `pypi` |

The corresponding GitHub environment is also named `pypi`, with deployment
restricted to `main`. Publication is initiated by running the release workflow
with an already qualified version tag. The workflow uses a short-lived OIDC
credential; it does not need a stored PyPI
API token. A pending publisher creates the project on first upload and does not
reserve its name beforehand. See [PyPI's setup instructions](https://docs.pypi.org/trusted-publishers/creating-a-project-through-oidc/).

## Prepare a release

1. Merge and qualify the release source. Update the version consistently using
   `scripts/verify_version.py` as the check. Keep existing release tags immutable.
2. Create the new version tag and publish its GitHub release. For example,
   `VERSION.txt` containing `1.2.0-rc.2` requires tag `v1.2.0-rc.2`.
3. Run **Publish Python package** from `main`, providing that release tag. The
   tag must identify the same main revision selected for the workflow. A version
   bump or tag move after qualification requires a new build.

The workflow resolves one immutable commit and uses it for every producer and
consumer checkout. It builds three repaired wheels and tests them on separate
clean runners, without a configured OCCT SDK. Installation tests exercise import,
meshing, controls, tags, exports, optional NumPy versions and the installed command.

Assembly checks each wheel's metadata and hashes, successful builder/consumer
receipts, matching native-CLI parity and the pinned OCCT source. Wheels are staged
separately from source assets. Matching source and qualification receipts are
attached to the public GitHub release before PyPI publication. Existing assets
are reused only when their hashes match; they are never silently replaced.

The `pypi` job rechecks the selected wheel hashes and uploads through Trusted
Publishing. It does not skip existing PyPI filenames or replace published files.
After upload, a separate job verifies all three index hashes and performs a fresh,
hash-pinned install from PyPI, then imports, meshes, exports and runs the command.
That job can be rerun independently if index propagation delays the first check.

## Installation

For the release candidate:

```sh
python -m pip install --pre paacs-mesh
```

For an exact qualified version:

```sh
python -m pip install paacs-mesh==1.2.0rc2
```

Stable releases use `python -m pip install paacs-mesh`. Wheel users do not need a
compiler or a separate CAD installation. Source builds and additional platforms
follow the [build guide](BUILDING.md) and require their build dependencies.
