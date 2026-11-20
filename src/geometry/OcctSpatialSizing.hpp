// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "OcctModelInternal.hpp"
#include "paacs/geometry/OcctBoundary.hpp"
#include "paacs/geometry/OcctMesher.hpp"
#include "paacs/sizing/RegionSizeField.hpp"
#include <BRepBndLib.hxx>
#include <Bnd_Box.hxx>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace paacs::geo::detail {

class OcctSpatialSizing {
public:
    OcctSpatialSizing(const OcctModel& model, const NativeCadMeshOptions& options,
                      sizing::RegionSizeField boxes)
        : model_(model), boxes_(std::move(boxes)) {
        faces_.reserve(options.faceSizes.size());
        for (const auto& control : options.faceSizes) {
            Bnd_Box bounds;
            BRepBndLib::AddOptimal(OcctModelAccess::face(model, control.faceTag), bounds, false, true);
            Vec3 lower{}, upper{};
            bounds.Get(lower.x, lower.y, lower.z, upper.x, upper.y, upper.z);
            faces_.push_back({control.faceTag, lower, upper,
                std::max(options.minimumSize, control.size),
                2 * (1 - 2 / (control.growth + 1))});
        }
    }

    std::uint64_t queryWork() const { return boxes_.boxQueryWorkBound() + faces_.size(); }
    const sizing::RegionSizeField& boxes() const { return boxes_; }

    double eval(Vec3 point, int onFace = -1) const {
        auto h = boxes_.eval(point);
        for (const auto& face : faces_) {
            if (face.size >= h) continue;
            if (face.tag == onFace || face.slope == 0.) {
                h = face.size;
                continue;
            }
            const auto separation = [](double p, double lo, double hi) {
                return std::max({0., lo - p, p - hi});
            };
            const auto lowerBound = std::hypot(separation(point.x, face.lower.x, face.upper.x),
                separation(point.y, face.lower.y, face.upper.y),
                separation(point.z, face.lower.z, face.upper.z));
            if (face.size + face.slope * lowerBound >= h) continue;
            Vec3 closest;
            if (!closestPointOnOcctFace(model_, face.tag, point, closest))
                throw std::runtime_error("native sizing could not project onto the selected trimmed face");
            h = std::min(h, face.size + face.slope * (point - closest).norm());
        }
        return h;
    }

private:
    struct Face {
        int tag;
        Vec3 lower, upper;
        double size, slope;
    };
    const OcctModel& model_;
    sizing::RegionSizeField boxes_;
    std::vector<Face> faces_;
};
}
