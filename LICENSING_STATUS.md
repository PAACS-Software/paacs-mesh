# Licensing

PAACS Software Ltd licenses the PAACS-authored software source under the
[Mozilla Public License 2.0](LICENSE). The [copyright notice](NOTICE) applies
to that source. MPL applies its source-sharing terms to covered files when
copies are distributed. A larger application can keep its own separate
source files proprietary. See the [Mozilla MPL FAQ](https://www.mozilla.org/en-US/MPL/2.0/FAQ/).

The source headers credit Paul LeDewhurst as project author and identify
PAACS Software Ltd and MPL 2.0.

The optional CAD adapter uses Open CASCADE Technology, which has its own
[LGPL-2.1 license](https://github.com/Open-Cascade-SAS/OCCT/blob/V8_0_1/LICENSE_LGPL_21.txt)
and [OCCT exception](https://github.com/Open-Cascade-SAS/OCCT/blob/V8_0_1/OCCT_LGPL_EXCEPTION.txt).
The source package contains a build recipe that fetches OCCT directly from
upstream and installs it as separate shared libraries in the local build tree.
The automatic build installs the two OCCT license files under
`build/step/occt/share/doc/opencascade`.
For a future binary distribution with OCCT, include the applicable notices
and satisfy the LGPL source and library replacement requirements.

Boost headers are used by the tests. NumPy, Matplotlib and SciPy are dependencies
of selected Python examples and visualization tools.

The STEP coverage benchmarks use the [ABC Dataset](https://deep-geometry.github.io/abc-dataset/).
The [NYU archive record](https://archive.nyu.edu/handle/2451/44309?mode=full)
describes its dataset license. The project page identifies the CAD model
creators as copyright holders and refers readers to the Onshape terms for
the models. The [asset provenance](docs/ASSET_PROVENANCE.md) identifies the
input cases, original figures and geometry fixtures included in the source package.
