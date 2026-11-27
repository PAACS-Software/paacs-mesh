# Python package releases

Version `1.2.1` uses the same source and PyPI version and the release
qualification gates below. Its matching GitHub release carries the publication
receipts. It targets regular,
64-bit CPython 3.12–3.15 on Windows x64/ARM64, Linux x64/ARM64 with glibc 2.28 or
newer, and macOS ARM64 with a minimum of macOS 15. Its metadata restricts Python
to `>=3.12,<3.16`.
Runtime CAD libraries are bundled. Only PAACS wheels enter PyPI; the matching
OCCT source archive and build recipes remain available with the GitHub release.

The previous published `1.2.0rc2` release supports CPython 3.12 on Windows x64,
macOS ARM64 and Linux x64 only, with metadata `>=3.12,<3.13`. Its tag and files
remain immutable, as do the later published RC3, RC4 and stable 1.2.0 tags and
files. Version 1.2.1 requires its own versioned build and fresh tagged-source
qualification; development artifacts are never relabelled as a different release.

Version 1.2.1 adds optional [pipeline access](PIPELINE.md) for qualified
boundaries, CAD-edge exports, pre-recovery classical point clouds and prescribed
point replay. It contains no learned model or training dependency and does not
change ordinary meshing algorithms or defaults. Its release checks must exercise
these staged interfaces as well as the existing import, sizing, audit and export
paths. The first staged scope remains one original solid with uniform sizing.

## Expanded qualification matrix

The release pipeline uses [python/compatibility.json](../python/compatibility.json)
for **20 repaired wheels and 28 independent installation results**: CPython
3.12–3.15 on Windows x64/ARM64, Linux x64/ARM64 and macOS ARM64. macOS wheels retain
deployment target 15.0 and are each installed on macOS 15, 26 and 27. The macOS 27
runner is a public preview, with its actual OS version checked and its result
required. Intel macOS and free-threaded Python are excluded.

The implementation passed all 55 compatibility jobs and seven native jobs in
[PR #11](https://github.com/PAACS-Software/paacs-mesh/pull/11), followed by a local
assembly of all 48 builder/consumer artifacts and strict Twine checks on all 20
wheels. The [compatibility evidence](PYTHON_COMPATIBILITY.md#qualification-evidence)
records those historical source identities and workflow runs. They do not
qualify 1.2.1; it must repeat qualification from its own versioned source before
publication. All required combinations and native regression checks remain
mandatory; a passing development PR does not
replace release qualification. See [Python compatibility](PYTHON_COMPATIBILITY.md)
for exact tags, NumPy versions, CAD SDK reuse and upstream runner references.

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
   `scripts/verify_version.py --final` for a stable release. Keep existing release
   tags immutable.
   Before tagging, record successful runs for all seven native qualification jobs
   and all 55 Python qualification jobs. Check that the tested checkout tree
   matches the tag target and merged main tree; a pull-request run can check out
   a virtual merge rather than its reported head commit. The publication workflow
   repeats the wheel matrix; it does not replace this native regression gate.
2. Create the new version tag and publish its GitHub release. For example,
   `VERSION.txt` containing `1.2.1` requires tag `v1.2.1` and publishes PyPI
   version `1.2.1`. Publish its GitHub release with `prerelease=false`.
3. Run **Publish Python package** from `main`, providing that release tag. The
   annotated tag must identify a commit on the history of the selected main
   revision. The tagged `VERSION.txt` determines the package version. Keep
   `qualified-run-id` empty for a new build.

For the separately checksummed PAACS source ZIP, run the following from the clean
qualified tag checkout. The flag preserves the tracked source manifest:

```sh
python scripts/package_source.py --preserve-source-manifest --output out/paacs-mesh-1.2.1.zip
```

Attach that ZIP and its `.sha256` file to the matching public GitHub release.
The publication workflow separately attaches the verified OCCT source, release
manifest and builder/consumer receipts; it does not generate this PAACS ZIP.

The workflow resolves one immutable release commit and uses it for every producer
and consumer checkout. Release operations run from the selected main revision,
which can contain a later packaging fix without changing the tagged engine.
It builds the complete manifest-defined wheel set and tests it on separate
clean runners, without a configured OCCT SDK. Installation tests exercise import,
meshing, controls, tags, exports, optional NumPy versions and the installed command.

Assembly checks each wheel's metadata and hashes, successful builder/consumer
receipts, actual compiler/configuration evidence, matching native-CLI parity and
the pinned OCCT source. Wheels are staged
separately from source assets. Matching source and qualification receipts are
attached to the public GitHub release before PyPI publication. Existing assets
are reused only when their hashes match; they are never silently replaced.

The `pypi` job rechecks the selected wheel hashes and uploads through Trusted
Publishing. It does not skip existing PyPI filenames or replace published files.
After upload, a separate job verifies every qualified wheel's index hash and performs a fresh,
hash-pinned install from PyPI, then imports, meshes, exports and runs the command.
That job can be rerun independently if index propagation delays the first check.

## Resume after an assembly failure

If all required qualification jobs passed but assembly failed before source assets
or wheels were published, a packaging fix can reuse those exact artifacts.
Merge the fix, then run **Publish Python package** from `main` with the unchanged
release tag and the previous publication run's numeric `qualified-run-id`.
Do not move the release tag or rebuild a different engine under that version.

Preparation verifies that the previous run is a completed `workflow_dispatch`
of this repository's publication workflow on `main`, at the tagged source
commit. Its preparation and the complete qualification job set must have passed;
assembly must have failed and source publication, PyPI upload and index checks
must have been skipped. A version already present on PyPI cannot be resumed.
Every required builder and consumer artifact must be unexpired and have GitHub SHA-256
digests. Downloads bind their immutable artifact IDs as well as the run ID.

Only new qualification is skipped. Assembly still checks the original wheel
and source hashes, exact clean source receipts, independent consumer tests,
metadata and native parity. Corrected release scripts run from the new workflow
revision while reading the original tag's version and OCCT pin. Source assets,
Trusted Publishing and the fresh PyPI installation then use the usual gates.
An attempted upload or incomplete qualification requires a separate diagnosis;
this resume option does not retry or overwrite published package files.

## Installation

Install the pinned stable release:

```sh
python -m pip install paacs-mesh==1.2.1
```

To allow prereleases when selecting the newest available version:

```sh
python -m pip install --pre paacs-mesh
```

Stable releases use `python -m pip install paacs-mesh`. Wheel users do not need a
compiler or a separate CAD installation. Source builds and additional platforms
follow the [build guide](BUILDING.md) and require their build dependencies.
