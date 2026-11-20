# Toolchain sources checked for this upgrade

- GCC C++ language support: https://gcc.gnu.org/projects/cxx-status.html
- GCC 16 changes: https://gcc.gnu.org/gcc-16/changes.html
- CMake download/release listing: https://cmake.org/download/
- Required language standard behavior: https://cmake.org/cmake/help/latest/prop_tgt/CXX_STANDARD.html
- Checkout action: https://github.com/actions/checkout
- Micromamba action: https://github.com/mamba-org/setup-micromamba

GCC labels its C++26 implementation experimental. This build uses the C++26
mode and tested library facilities without requiring experimental modules,
reflection, or contracts. CMake's compiler feature negotiation plus
CXX_STANDARD_REQUIRED prevents silent language downgrade.

Actual installed versions and local execution outcomes are in TOOLCHAIN_LOCK.json
and the root VALIDATION.md; web support tables are not substitutes for those
local tests.
