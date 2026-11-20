// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/mesh/LinearVolumeAudit.hpp"
namespace paacs::mesh {
LinearVolumeAudit referenceLinearVolumeBoundary(const Mesh&, const LinearVolumeAuditLimits& = {}) noexcept;
}
