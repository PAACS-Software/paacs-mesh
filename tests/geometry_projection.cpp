// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/geometry/PolyhedralModel.hpp"
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
void require(bool good, const char *message) {
    if (!good) throw std::runtime_error(message);
}
class BrokenFace final : public paacs::geo::GFace {
public:
    enum class Failure { NonfiniteUv, KernelException, Allocation };
    explicit BrokenFace(Failure failure) : GFace(1), failure_(failure) {}
    mutable unsigned evaluations = 0;
    paacs::Vec2 paramAt(const paacs::Vec3&) const override {
        if (failure_ == Failure::KernelException) throw std::runtime_error("kernel failure");
        if (failure_ == Failure::Allocation) throw std::bad_alloc();
        return {std::numeric_limits<double>::quiet_NaN(), 0};
    }
    paacs::Vec3 pointAt(const paacs::Vec2&) const override { ++evaluations; return {}; }
    paacs::Vec3 normalAt(const paacs::Vec2&) const override { return {0, 0, 1}; }
    void paramBounds(paacs::Vec2& lo, paacs::Vec2& hi) const override { lo = {}; hi = {1, 1}; }
    std::vector<int> boundaryEdgeTags() const override { return {}; }
    double curvatureAt(const paacs::Vec2&) const override { return 0; }
private:
    Failure failure_;
};
}
int main() {
    using namespace paacs;
    try {
        const std::vector<std::vector<Vec2>> loops{
            {{0, 0}, {3, 0}, {3, 3}, {0, 3}},
            {{1, 1}, {1, 2}, {2, 2}, {2, 1}}};
        const auto model = geo::makePolyhedralModel(geo::extrudedPolygon(loops, {0, 0, 1}));
        require(model.built(), "projection fixture");
        bool checkedCap = false;
        for (const auto *face : model.model->faces()) {
            const auto normal = face->normalAt(face->paramAt({0, 0, 1}));
            if (normal.z < .99) continue;
            checkedCap = true;
            Vec2 uv;
            Vec3 projected;
            require(face->projectPoint({.5, .5, 2}, uv, projected), "inside trimmed cap projects");
            require((projected - Vec3{.5, .5, 1}).norm2() < 1e-20, "projection point");
            require(!face->projectPoint({1.5, 1.5, 2}, uv, projected), "hole projection refused");
            require(!face->projectPoint({4, .5, 2}, uv, projected), "outside trim projection refused");
            require(!face->projectPoint({std::numeric_limits<double>::quiet_NaN(), .5, 2}, uv, projected),
                    "nonfinite projection refused");
            require(!face->projectPoint({std::numeric_limits<double>::infinity(), .5, 2}, uv, projected),
                    "infinite projection refused");
        }
        require(checkedCap, "upper cap exercised");
        for (const auto failure : {BrokenFace::Failure::NonfiniteUv, BrokenFace::Failure::KernelException}) {
            BrokenFace face(failure);
            Vec2 uv{7, 8};
            Vec3 projected{9, 10, 11};
            require(!face.projectPoint({}, uv, projected) && face.evaluations == 0,
                    "malformed backend projection refused before evaluation");
            require(uv.u == 7 && uv.v == 8 && projected.x == 9 && projected.y == 10 && projected.z == 11,
                    "failed projection outputs changed");
        }
        BrokenFace allocation(BrokenFace::Failure::Allocation);
        bool propagated = false;
        try {
            Vec2 uv;
            Vec3 projected;
            allocation.projectPoint({}, uv, projected);
        } catch (const std::bad_alloc&) {
            propagated = true;
        }
        require(propagated, "allocation failure lost its diagnostic");
        std::cout << "Trimmed face projection and nonfinite input regressions passed\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
