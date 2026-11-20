// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Private optional-adapter support for separate, coarse native box meshes.
// This does not weld body nodes or approve a general bonded interface.
#include "paacs/geometry/OcctModel.hpp"
#include "paacs/algo/MixedMeshQuality.hpp"
#include "paacs/mesh/LinearVolumeAudit.hpp"
#include <atomic>
#include <memory>

namespace paacs::geo::detail {
enum class OcctBoxContactReason : std::uint8_t {
    None, InvalidRequest, UnsupportedBody, PartialFaceContact, MaterialOverlap,
    NativeBooleanFailed, GeometryError, UnsupportedResolution, InvalidMesh,
    BudgetExceeded, Cancelled, AllocationFailed, SourceChanged, InternalError
};
// Optional per-meshing-call cache. Use on the calling thread only. Cached
// certificates apply only to the same imported model and exact source bytes.
class OcctBoxCertificateSession {
public:
    OcctBoxCertificateSession() noexcept;
    ~OcctBoxCertificateSession();
    OcctBoxCertificateSession(const OcctBoxCertificateSession&) = delete;
    OcctBoxCertificateSession& operator=(const OcctBoxCertificateSession&) = delete;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
    friend struct OcctBoxCertificateAccess;
};
struct OcctBoxContactOptions {
    std::uint64_t maximumWork = 2'000'000;
    std::uint64_t maximumSnapshotBytes = 64ULL * 1024 * 1024;
    std::uint64_t maximumWorkingBytes = 512ULL * 1024 * 1024;
    const std::atomic<bool>* cancellation = nullptr;
    OcctBoxCertificateSession* certificateSession = nullptr;
};
struct OcctBoxFaceContactOutcome {
    bool admitted = false, sourceImmutabilityVerified = false;
    OcctBoxContactReason reason = OcctBoxContactReason::InternalError;
    int firstBodyTag = -1, secondBodyTag = -1;
    std::uint64_t workUsed = 0;
    double nativeSkinArea = 0, expectedFullFaceArea = 0;
};
// Certifies both actual native solids as axis-aligned boxes with exact bbox
// corners and complete native expected-box Common volume. Material overlap,
// partial faces and zero-area contacts are refused. Native operations use deep
// copies. The caller must still pair every contact triangle exactly, uniquely
// and with opposite winding, and compare its summed area with nativeSkinArea.
OcctBoxFaceContactOutcome admitOcctBoxFaceContact(const OcctModel&, int firstBodyTag,
    int secondBodyTag, const OcctBoxContactOptions& = {}) noexcept;

struct OcctBoxCornerEvaluationOutcome {
    bool published = false, destinationUnchanged = true;
    bool sourceImmutabilityVerified = false;
    OcctBoxContactReason reason = OcctBoxContactReason::InternalError;
    std::uint64_t workUsed = 0, evaluatedVertices = 0, nativeInclusionSamples = 0;
    double maximumEndpointCorrection = 0;
    double maximumBoundaryNodeError = 0, maximumBoundarySampleError = 0;
    double cadVolume = 0, meshVolume = 0;
    mesh::LinearVolumeAudit audit;
    algo::MixedMeshQualityReport quality;
};
// Only an eight-node, six-Tet4/twelve-Tri3 body is supported. Each distinct node
// must already carry a unique owned original CAD vertex tag. Coordinates are
// evaluated from that identity, never inferred by proximity. Counts, cell
// connectivity and all node/cell classifications are preserved. Exact positive
// Tet signs, native node/centroid inclusion, original face coverage/projection,
// incidence, volume and quality are checked before one transactional swap.
// Call only for bodies belonging to an admitted full-face contact; the caller
// retains whole-output/contact/audit/quality/cancellation publication gates.
// Work counts explicit PAACS visits/native calls; cancellation is cooperative.
// Logical snapshot/working bounds exclude OCCT allocation and process RSS.
OcctBoxCornerEvaluationOutcome evaluateOcctBoxContactCorners(const OcctModel&,
    int bodyTag, const OcctBoxContactOptions&, mesh::Mesh& body) noexcept;
const char* occtBoxContactReasonName(OcctBoxContactReason) noexcept;
} // namespace paacs::geo::detail
