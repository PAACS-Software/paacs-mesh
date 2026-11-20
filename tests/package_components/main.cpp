// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/geometry/Kernel.hpp"

int main() {
    const auto kernel = paacs::geo::makeBuiltinKernel();
    const auto box = kernel->makeBox({0, 0, 0}, {1, 2, 3});
    if (!box || box->solids().size() != 1 || box->faces().size() != 6)
        return 1;
    return box->solids().front()->isInside({0.5, 1, 1.5}) ? 0 : 2;
}
