# Build and test V0.10

## Release build

```bash
CC=clang CXX=clang++ cmake -S source/standalone -B build-v010 \
  -DCMAKE_BUILD_TYPE=Release \
  -DPAACS_V010_BUILD_CLI=ON \
  -DPAACS_V010_BUILD_TEST_PROBES=ON \
  -DPAACS_V010_BUILD_BENCHMARKS=ON \
  -DPAACS_V010_WARNINGS_AS_ERRORS=ON
cmake --build build-v010 -j
ctest --test-dir build-v010 --output-on-failure
```

## Sanitizer build

```bash
CC=clang CXX=clang++ cmake -G Ninja -S source/standalone -B build-v010-san \
  -DCMAKE_BUILD_TYPE=Debug \
  -DPAACS_V010_BUILD_CLI=OFF \
  -DPAACS_V010_BUILD_TEST_PROBES=ON \
  -DPAACS_V010_BUILD_BENCHMARKS=OFF \
  -DPAACS_V010_WARNINGS_AS_ERRORS=ON \
  -DCMAKE_CXX_FLAGS='-O0 -g -fno-omit-frame-pointer -fsanitize=address,undefined' \
  -DCMAKE_EXE_LINKER_FLAGS='-fsanitize=address,undefined'
ninja -C build-v010-san
ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 \
UBSAN_OPTIONS=halt_on_error=1 \
ctest --test-dir build-v010-san \
  -R 'pre_occt_core|delaunay_goldens' --output-on-failure
```

## Install and consume

```bash
cmake --install build-v010 --prefix /path/to/install
```

Then:

```cmake
find_package(PAACSMeshV010 CONFIG REQUIRED)
target_link_libraries(my_target PRIVATE PAACS::mesh)
```

## CLI examples

```bash
paacs-mesh-v010 tet-box 0.45 box.paacsmesh --threads 1
paacs-mesh-v010 tet-box 0.45 box10.paacsmesh --tet10
paacs-mesh-v010 tet-cylinder 1.0 1.5 0.45 cylinder.paacsmesh
paacs-mesh-v010 audit box10.paacsmesh
paacs-mesh-v010 hex 16 12 8 hex.paacsmesh --threads 4
```

## Patch application

From the package directory:

```bash
./verify_patch.sh /path/to/clean/v093/source
./apply.sh /path/to/clean/v093/source
```
