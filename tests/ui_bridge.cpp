// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/mesh/MeshAudit.hpp"
#include "paacs/sdk/JobService.hpp"
#include "paacs/sdk/MeshArchive.hpp"
#include "paacs/sdk/UIBridge.h"
#include "paacs/sizing/ControlField.hpp"
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
using namespace paacs;
namespace {
int checks = 0;
void require(bool b, const char *text) {
    ++checks;
    if (!b)
        throw std::runtime_error(text);
}
using Request = std::unique_ptr<paacs_ui_request, decltype(&paacs_ui_request_destroy)>;
using Service = std::unique_ptr<paacs_ui_service, decltype(&paacs_ui_service_destroy)>;
using Job = std::unique_ptr<paacs_ui_job, decltype(&paacs_ui_job_release)>;
using Result = std::unique_ptr<paacs_ui_result, decltype(&paacs_ui_result_release)>;
Request request(const char *route) {
    return {paacs_ui_request_create(1, route), paacs_ui_request_destroy};
}
void set(paacs_ui_request *r, const char *key, double value) {
    require(paacs_ui_request_set(r, key, value) == PAACS_UI_OK, "set option");
}
Job submit(paacs_ui_service *s, const paacs_ui_request *r) {
    paacs_ui_job *job = nullptr;
    char error[192];
    auto status = paacs_ui_submit(s, r, &job, error, sizeof(error));
    if (status)
        std::cerr << error << '\n';
    require(status == PAACS_UI_OK && job, "submitted through C ABI");
    return {job, paacs_ui_job_release};
}
Result complete(const paacs_ui_job *j) {
    uint32_t done = 0;
    require(paacs_ui_job_wait(j, 20000, &done) == PAACS_UI_OK && done, "job completed");
    paacs_ui_snapshot snapshot{};
    require(paacs_ui_job_snapshot(j, &snapshot) == PAACS_UI_OK, "snapshot");
    if (snapshot.state != 3)
        std::cerr << snapshot.diagnostic << '\n';
    require(snapshot.state == 3, "job succeeded");
    paacs_ui_result *r = nullptr;
    require(paacs_ui_job_result(j, &r) == PAACS_UI_OK && r, "result handle");
    return {r, paacs_ui_result_release};
}
void compare(const paacs_ui_result *r, const mesh::Mesh &direct) {
    paacs_ui_summary summary{};
    require(paacs_ui_result_summary(r, &summary) == PAACS_UI_OK, "summary");
    require(summary.hash_available && summary.hash == mesh::meshHash(direct),
            "opt-in strict hash agrees");
    require(summary.nodes == static_cast<uint64_t>(direct.nodeCount()) &&
                summary.cells == static_cast<uint64_t>(direct.cellCount()),
            "mesh counts identical");
    // Every C page entry checked against the independent native invocation.
    std::array<paacs_ui_node, 37> nodes;
    for (uint64_t offset = 0; offset < summary.nodes; offset += nodes.size()) {
        size_t n = 0;
        require(paacs_ui_result_nodes(r, offset, nodes.data(), nodes.size(), &n) == PAACS_UI_OK,
                "node page");
        for (size_t i = 0; i < n; ++i) {
            const auto id = static_cast<Index>(offset + i);
            const auto p = direct.node(id);
            const auto c = direct.nodeClass(id);
            require(nodes[i].position[0] == p.x && nodes[i].position[1] == p.y &&
                        nodes[i].position[2] == p.z,
                    "same node coordinates");
            require(nodes[i].dimension == c.dim && nodes[i].tag == c.tag,
                    "node classification preserved");
        }
    }
    std::array<paacs_ui_cell, 29> cells;
    for (uint64_t offset = 0; offset < summary.cells; offset += cells.size()) {
        size_t n = 0;
        require(paacs_ui_result_cells(r, offset, cells.data(), cells.size(), &n) == PAACS_UI_OK,
                "cell page");
        for (size_t i = 0; i < n; ++i) {
            const auto id = static_cast<Index>(offset + i);
            const auto &c = direct.cell(id);
            const auto &v = cells[i];
            require(v.type == static_cast<uint32_t>(c.type) &&
                        v.node_count == static_cast<uint32_t>(mesh::nodesPerCell(c.type)),
                    "same cell type");
            require(v.dimension == c.onEntity.dim && v.tag == c.onEntity.tag,
                    "cell classification preserved");
            for (size_t k = 0; k < v.node_count; ++k)
                require(v.nodes[k] == static_cast<uint64_t>(direct.cellNodes(id)[k]),
                        "same connectivity");
        }
    }
    size_t written = 123;
    require(paacs_ui_result_nodes(r, summary.nodes, nullptr, 0, &written) == PAACS_UI_OK &&
                written == 0,
            "empty end page");
    require(paacs_ui_result_cells(r, UINT64_MAX, cells.data(), 1, &written) ==
                    PAACS_UI_OUT_OF_RANGE &&
                written == 0,
            "overflow-safe page offset");
    require(paacs_ui_result_nodes(r, 0, nodes.data(), PAACS_UI_MAX_PAGE + 1, &written) ==
                PAACS_UI_LIMIT,
            "page cap enforced before write");
    require(paacs_ui_result_nodes(r, 0, nullptr, 1, &written) == PAACS_UI_INVALID,
            "null page rejected");
}
void descriptors() {
    for (const char *route : {"primitive", "adaptive_box"}) {
        auto r = request(route);
        require(bool(r), "create request");
        const auto n = paacs_ui_option_count(r.get());
        for (size_t i = 0; i < n; ++i) {
            paacs_ui_option o{};
            double current = 0;
            require(paacs_ui_option_at(r.get(), i, &o) == PAACS_UI_OK, "option descriptor");
            require(paacs_ui_request_get(r.get(), o.key, &current) == PAACS_UI_OK &&
                        current == o.default_value,
                    "descriptor default matches implementation");
            set(r.get(), o.key, current);
            require(
                paacs_ui_request_set(r.get(), o.key, std::numeric_limits<double>::quiet_NaN()) ==
                    PAACS_UI_INVALID,
                "NaN rejected");
            require(paacs_ui_request_set(r.get(), o.key, std::nextafter(o.maximum, INFINITY)) ==
                        PAACS_UI_INVALID,
                    "upper bound rejected");
            if (o.kind)
                require(paacs_ui_request_set(r.get(), o.key, .5) == PAACS_UI_INVALID,
                        "fractional integer rejected");
            require(paacs_ui_request_get(r.get(), o.key, &current) == PAACS_UI_OK &&
                        current == o.default_value,
                    "failed setter preserves request");
        }
        require(paacs_ui_request_set(r.get(), "typo", 1) == PAACS_UI_INVALID,
                "unknown option rejected");
    }
    require(!paacs_ui_request_create(2, "primitive") && !paacs_ui_request_create(1, "missing"),
            "version/route mismatch rejected");
    require(!paacs_ui_service_create(1, 0, 1, 1), "invalid service limits");
}
void roundTrips() {
    Service s(paacs_ui_service_create_with_flags(1, 8, 64 * 1024 * 1024, 128 * 1024 * 1024,
                                                 PAACS_UI_COMPUTE_HASH),
              paacs_ui_service_destroy);
    require(bool(s), "service created");
    for (int elements = 0; elements < 3; ++elements) {
        auto r = request("primitive");
        set(r.get(), "elements", elements);
        set(r.get(), "target_size", .4);
        set(r.get(), "maximum_chord_error", .03);
        set(r.get(), "warning_quality", .99);
        auto j = submit(s.get(), r.get());
        set(r.get(), "radius", -0.0);
        r.reset(); // Job keeps original request.
        auto result = complete(j.get());
        algo::PrimitiveMeshingOptions o;
        o.elements = static_cast<algo::PrimitiveElements>(elements);
        o.targetSize = .4;
        o.maximumChordError = .03;
        o.warningQuality = .99;
        mesh::Mesh direct;
        require(algo::meshPrimitive({}, o, direct).published, "direct primitive");
        compare(result.get(), direct);
        paacs_ui_summary sum{};
        require(paacs_ui_result_summary(result.get(), &sum) == PAACS_UI_OK && sum.quality_warning,
                "warning preserved");
    }
    {
        auto r = request("primitive");
        set(r.get(), "height", 1);
        set(r.get(), "elements", 1);
        set(r.get(), "target_size", .4);
        set(r.get(), "maximum_chord_error", .04);
        set(r.get(), "profile_global_size", .4);
        set(r.get(), "profile_enabled", 1);
        paacs_ui_profile_control control{0, 0, .2, 1.3};
        require(paacs_ui_request_profile_controls(r.get(), &control, 1) == PAACS_UI_OK,
                "profile controls");
        auto j = submit(s.get(), r.get());
        r.reset();
        auto result = complete(j.get());
        sizing::CappedSizeSpecification ps;
        ps.globalSize = .4;
        ps.controls = {{0, 0, .2, 1.3}};
        auto p = sizing::compileCappedSizeProfile(ps);
        algo::PrimitiveShape shape;
        shape.height = 1;
        algo::PrimitiveMeshingOptions o;
        o.elements = algo::PrimitiveElements::Tet4;
        o.targetSize = .4;
        o.maximumChordError = .04;
        o.sizeProfile = &p.profile;
        mesh::Mesh direct;
        require(algo::meshPrimitive(shape, o, direct).published, "direct profiled primitive");
        compare(result.get(), direct);
    }
    {
        auto r = request("adaptive_box");
        set(r.get(), "target_size", .6);
        paacs_ui_field_node nodes[2]{};
        nodes[0].operation = 5;
        nodes[0].first = nodes[0].second = UINT32_MAX;
        nodes[0].a[0] = nodes[0].a[1] = nodes[0].a[2] = .5;
        nodes[0].parameters[0] = .1;
        nodes[1].operation = 6;
        nodes[1].first = 0;
        nodes[1].second = UINT32_MAX;
        nodes[1].parameters[0] = 0;
        nodes[1].parameters[1] = .35;
        nodes[1].parameters[2] = .2;
        nodes[1].parameters[3] = .6;
        require(paacs_ui_request_field(r.get(), nodes, 2, 1, .001, 1, 1) == PAACS_UI_OK,
                "set field");
        auto j = submit(s.get(), r.get());
        r.reset();
        auto result = complete(j.get());
        sizing::ControlFieldSpecification fs;
        fs.nodes = {sizing::field::distanceBall({.5, .5, .5}, .1),
                    sizing::field::threshold(0, 0, .35, .2, .6)};
        fs.output = 1;
        auto f = sizing::compileControlField(fs);
        algo::AdaptiveBoxRequest o;
        o.targetSize = .6;
        o.field = &f.field;
        mesh::Mesh direct;
        require(algo::meshAdaptiveBox(o, direct).published, "direct adaptive field");
        compare(result.get(), direct);
        // C result is independently owned after service and job destruction.
        j.reset();
        s.reset();
        compare(result.get(), direct);
    }
}
} // namespace
int main() {
    try {
        descriptors();
        roundTrips();
        std::cout << "{\"checks\":" << checks << ",\"direct_mesh_pairs\":5,\"passed\":true}\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
