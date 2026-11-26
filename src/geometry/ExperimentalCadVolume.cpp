// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0
#include "paacs/geometry/ExperimentalCadVolume.hpp"
#include "OcctFrontMesher.hpp"

namespace paacs::geo {
ExperimentalCadBoundaryOutcome prepareExperimentalCadBoundary(const OcctModel& model,
    const NativeCadMeshOptions& options,ExperimentalCadBoundary& destination) noexcept {
    return detail::prepareExperimentalCadBoundaryImpl(model,nullptr,options,destination);
}
ExperimentalCadBoundaryOutcome rebindExperimentalCadBoundary(const OcctModel& model,
    const mesh::Mesh& boundary,const NativeCadMeshOptions& options,
    ExperimentalCadBoundary& destination) noexcept {
    return detail::prepareExperimentalCadBoundaryImpl(model,&boundary,options,destination);
}
ExperimentalCadVolumeOutcome replayExperimentalCadVolume(const OcctModel& model,
    const ExperimentalCadBoundary& boundary,const NativeCadMeshOptions& options,
    ExperimentalInteriorMode mode,std::span<const Vec3> points,mesh::Mesh& destination) noexcept {
    return detail::replayExperimentalCadVolumeImpl(model,boundary,options,mode,points,destination);
}
} // namespace paacs::geo
