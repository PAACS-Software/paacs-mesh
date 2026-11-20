// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Private optional route for topology-gated paired curved sheets.
#include "OcctNativeSweep.hpp"
namespace paacs::geo::detail {
struct OcctPairedShellOptions : OcctNativeSweepOptions {
    unsigned maximumNativeValidationThreads = 1; // Explicit cap1..8.
};
struct OcctPairedShellOutcome : OcctNativeSweepOutcome {
    unsigned nativeValidationThreads = 0;
};
using OcctPairedShellReason = OcctNativeSweepReason;
// Exactly two curved disk faces with native collapsed poles and closed rims,
// plus their connecting face. Actual owned topology establishes every tag;
// nearest trimmed projection proposes correspondence, never entity identity.
// The pulling Tet4 candidate must pass exact exterior embedding, full oriented
// volume incidence, every original vertex/face, owning-entity point tolerances,
// caller facet error, stable original GK volume, and every native node/Tet
// centroid before a nonthrowing swap. No uniform-extrusion claim is made.
// Source topology and tolerances remain immutable. Counterpart parameter seams
// may lie inside a face triangulation; the immutable source still retains their
// native edge/use identities. Work and working bytes cover logical PAACS visits
// and payloads, not native solver internals or process memory/time.
OcctPairedShellOutcome buildOcctPairedShellVolume(const OcctModel &, int solidTag,
                                                  const OcctPairedShellOptions &,
                                                  mesh::Mesh &destination) noexcept;
} // namespace paacs::geo::detail
