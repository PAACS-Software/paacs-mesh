// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/geometry/OcctImport.hpp"

int main() {
    const auto result = paacs::geo::importPlanarCad("missing-package-check.step");
    return result.reason == paacs::geo::CadImportReason::FileUnavailable ? 0 : 1;
}
