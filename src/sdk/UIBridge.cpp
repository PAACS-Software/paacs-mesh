// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/sdk/UIBridge.h"
#include "paacs/sdk/JobService.hpp"
#include "paacs/sdk/MeshArchive.hpp"
#include "paacs/sdk/Version.hpp"
#include "paacs/sizing/ControlField.hpp"
#include <cmath>
#include <cstring>
#include <optional>
#include <type_traits>

using namespace paacs;
struct paacs_ui_service {
    sdk::JobService service;
    explicit paacs_ui_service(sdk::JobServiceLimits l) : service(l) {}
};
struct paacs_ui_request {
    bool primitive = true, profileEnabled = false;
    algo::PrimitiveShape shape;
    algo::PrimitiveMeshingOptions p;
    algo::AdaptiveBoxRequest a;
    sizing::CappedSizeSpecification profile;
    std::vector<double> positions;
    std::optional<sizing::ControlFieldSpecification> field;
};
struct paacs_ui_job {
    sdk::JobHandle job;
};
struct paacs_ui_result {
    std::shared_ptr<const sdk::JobResult> result;
};
namespace {
using Request = paacs_ui_request;
struct Entry {
    const char *name;
    double minimum, maximum;
    const char *choices;
    std::uint32_t kind;
    double (*get)(const Request &);
    void (*set)(Request &, double);
};
template <class T> constexpr std::uint32_t kind() {
    if constexpr (std::is_same_v<T, bool>)
        return 2;
    if constexpr (std::is_integral_v<T> || std::is_enum_v<T>)
        return 1;
    return 0;
}
template <class T> void assign(T &out, double v) {
    if constexpr (std::is_enum_v<T>)
        out = static_cast<T>(static_cast<int>(v));
    else
        out = static_cast<T>(v);
}
constexpr double maxInteger = 9007199254740991.0;
// The same table owns defaults, bounds, discovery and typed setters.
#define OPT(key, member, lo, hi, choices)                                                          \
    Entry {                                                                                        \
        key, lo, hi, choices, kind<decltype(Request{}.member)>(),                                  \
            +[](const Request &r) { return static_cast<double>(r.member); },                       \
            +[](Request &r, double v) { assign(r.member, v); }                                     \
    }
#define REAL(key, member) OPT(key, member, 0, 1e100, "")
#define LIMIT(key, member) OPT(key, member, 0, maxInteger, "")
#define VEC(key, member)                                                                           \
    OPT(key "_x", member.x, -1e100, 1e100, ""), OPT(key "_y", member.y, -1e100, 1e100, ""),        \
        OPT(key "_z", member.z, -1e100, 1e100, "")
constexpr Entry primitiveOptions[]{
    OPT("kind", shape.kind, 0, 5, "cylinder,cone,frustum,tube,torus,elbow"),
    VEC("origin", shape.origin),
    VEC("axis", shape.axis),
    VEC("x_direction", shape.xDirection),
    REAL("radius", shape.radius),
    REAL("height", shape.height),
    REAL("end_radius", shape.endRadius),
    REAL("inner_radius", shape.innerRadius),
    REAL("major_radius", shape.majorRadius),
    REAL("sweep_angle", shape.sweepAngle),
    REAL("target_size", p.targetSize),
    REAL("maximum_chord_error", p.maximumChordError),
    REAL("maximum_relative_volume_error", p.maximumRelativeVolumeError),
    OPT("minimum_quality", p.minimumQuality, 0, 1, ""),
    OPT("warning_quality", p.warningQuality, 0, 1, ""),
    OPT("elements", p.elements, 0, 2, "hex_dominant,tet4,tet10"),
    OPT("cone_strategy", p.coneStrategy, 0, 2, "automatic,swept,graded"),
    OPT("cylinder_strategy", p.cylinderStrategy, 0, 2, "automatic,swept,uniform"),
    OPT("circular_strategy", p.circularStrategy, 0, 2, "automatic,swept,uniform"),
    OPT("quality_policy", p.qualityPolicy, 0, 1, "warn,reject"),
    OPT("threads", p.threads, 0, 2147483647, ""),
    LIMIT("maximum_nodes", p.maximumNodes),
    LIMIT("maximum_cells", p.maximumCells),
    LIMIT("maximum_output_bytes", p.maximumOutputBytes),
    LIMIT("maximum_working_bytes", p.maximumWorkingBytes),
    LIMIT("audit_maximum_face_uses", p.auditLimits.maximumFaceUses),
    LIMIT("audit_maximum_temporary_bytes", p.auditLimits.maximumTemporaryBytes),
    OPT("profile_strategy", p.profileStrategy, 0, 1, "sweep,volume"),
    REAL("maximum_axial_growth", p.maximumAxialGrowth),
    LIMIT("maximum_profile_section_nodes", p.maximumProfileSectionNodes),
    LIMIT("maximum_profile_slab_work", p.maximumProfileSlabWork),
    OPT("reuse_profile_slabs", p.reuseProfileSlabs, 0, 1, ""),
    LIMIT("profile_maximum_cells", p.profileLimits.maximumCells),
    LIMIT("profile_distribution_maximum_bytes", p.profileLimits.maximumBytes),
    REAL("profile_maximum_adjacent_growth", p.profileLimits.maximumAdjacentGrowth),
    OPT("profile_enabled", profileEnabled, 0, 1, ""),
    REAL("profile_length", profile.length),
    REAL("profile_global_size", profile.globalSize),
    REAL("profile_minimum_size", profile.minimumSize),
    OPT("profile_maximum_controls", profile.maximumControls, 0, 1024, ""),
    LIMIT("profile_maximum_knots", profile.maximumKnots),
    LIMIT("profile_maximum_work", profile.maximumWork),
    LIMIT("profile_maximum_bytes", profile.maximumBytes)};
constexpr Entry adaptiveOptions[]{
    VEC("lower", a.lower),
    VEC("upper", a.upper),
    REAL("target_size", a.targetSize),
    OPT("maximum_depth", a.maximumDepth, 0, UINT32_MAX, ""),
    LIMIT("maximum_leaves", a.maximumLeaves),
    LIMIT("maximum_nodes", a.maximumNodes),
    LIMIT("maximum_cells", a.maximumCells),
    LIMIT("maximum_output_bytes", a.maximumOutputBytes),
    LIMIT("maximum_working_bytes", a.maximumWorkingBytes),
    OPT("minimum_quality", a.minimumQuality, 0, 1, ""),
    OPT("threads", a.threads, 0, 2147483647, ""),
    LIMIT("audit_maximum_face_uses", a.auditLimits.maximumFaceUses),
    LIMIT("audit_maximum_temporary_bytes", a.auditLimits.maximumTemporaryBytes),
    REAL("maximum_adjacent_size_ratio", a.maximumAdjacentSizeRatio),
    LIMIT("maximum_sizing_work", a.maximumSizingWork)};
#undef VEC
#undef LIMIT
#undef REAL
#undef OPT
std::span<const Entry> entries(const Request &r) {
    if (r.primitive)
        return primitiveOptions;
    return adaptiveOptions;
}
const Entry *find(const Request &r, const char *key) {
    if (key)
        for (const auto &e : entries(r))
            if (std::strcmp(e.name, key) == 0)
                return &e;
    return nullptr;
}
void diagnostic(char *buffer, std::size_t capacity, const char *text) {
    if (!buffer || !capacity)
        return;
    const auto n = std::min(capacity - 1, std::strlen(text));
    std::memcpy(buffer, text, n);
    buffer[n] = 0;
}
Vec3 vec(const double *a) {
    return {a[0], a[1], a[2]};
}
void vec(double *out, Vec3 p) {
    out[0] = p.x;
    out[1] = p.y;
    out[2] = p.z;
}
} // namespace
extern "C" {
uint32_t paacs_ui_abi_version(void) {
    return PAACS_UI_ABI_VERSION;
}
const char *paacs_ui_release_version(void) {
    return sdk::kApiVersion;
}
size_t paacs_ui_abi_size(uint32_t t) {
    constexpr std::size_t sizes[]{sizeof(paacs_ui_option),    sizeof(paacs_ui_snapshot),
                                  sizeof(paacs_ui_summary),   sizeof(paacs_ui_node),
                                  sizeof(paacs_ui_cell),      sizeof(paacs_ui_profile_control),
                                  sizeof(paacs_ui_field_node)};
    return t < std::size(sizes) ? sizes[t] : 0;
}
paacs_ui_service *paacs_ui_service_create(uint32_t abi, uint32_t jobs, uint64_t bytes,
                                          uint64_t retained) {
    return paacs_ui_service_create_with_flags(abi, jobs, bytes, retained, 0);
}
paacs_ui_service *paacs_ui_service_create_with_flags(uint32_t abi, uint32_t jobs, uint64_t bytes,
                                                     uint64_t retained, uint32_t flags) {
    if (abi != PAACS_UI_ABI_VERSION || (flags & ~PAACS_UI_COMPUTE_HASH))
        return nullptr;
    try {
        return new paacs_ui_service({jobs, bytes, retained, bool(flags & PAACS_UI_COMPUTE_HASH)});
    } catch (...) {
        return nullptr;
    }
}
void paacs_ui_service_destroy(paacs_ui_service *s) {
    delete s;
}
uint64_t paacs_ui_service_retained_bytes(const paacs_ui_service *s) {
    return s ? s->service.retainedResultBytes() : 0;
}
paacs_ui_request *paacs_ui_request_create(uint32_t abi, const char *route) {
    if (abi != PAACS_UI_ABI_VERSION || !route ||
        (std::strcmp(route, "primitive") && std::strcmp(route, "adaptive_box")))
        return nullptr;
    try {
        auto *r = new Request;
        r->primitive = std::strcmp(route, "primitive") == 0;
        return r;
    } catch (...) {
        return nullptr;
    }
}
void paacs_ui_request_destroy(paacs_ui_request *r) {
    delete r;
}
size_t paacs_ui_option_count(const paacs_ui_request *r) {
    return r ? entries(*r).size() : 0;
}
paacs_ui_status paacs_ui_option_at(const paacs_ui_request *r, size_t i, paacs_ui_option *out) {
    if (!r || !out)
        return PAACS_UI_INVALID;
    auto options = entries(*r);
    if (i >= options.size())
        return PAACS_UI_OUT_OF_RANGE;
    Request defaults;
    const auto &e = options[i];
    *out = {e.name, e.kind, e.get(defaults), e.minimum, e.maximum, e.choices};
    return PAACS_UI_OK;
}
paacs_ui_status paacs_ui_request_get(const paacs_ui_request *r, const char *key, double *out) {
    if (!r || !out)
        return PAACS_UI_INVALID;
    auto *e = find(*r, key);
    if (!e)
        return PAACS_UI_INVALID;
    *out = e->get(*r);
    return PAACS_UI_OK;
}
paacs_ui_status paacs_ui_request_set(paacs_ui_request *r, const char *key, double v) {
    if (!r)
        return PAACS_UI_INVALID;
    const auto *e = find(*r, key);
    if (!e || !std::isfinite(v) || v < e->minimum || v > e->maximum ||
        (e->kind && std::trunc(v) != v))
        return PAACS_UI_INVALID;
    e->set(*r, v);
    return PAACS_UI_OK;
}
paacs_ui_status paacs_ui_request_axial_positions(paacs_ui_request *r, const double *values,
                                                 size_t n) {
    if (!r || !r->primitive || (n && !values))
        return PAACS_UI_INVALID;
    if (n > 100001)
        return PAACS_UI_LIMIT;
    try {
        std::vector<double> copy;
        if (n)
            copy.assign(values, values + n);
        for (auto v : copy)
            if (!std::isfinite(v))
                return PAACS_UI_INVALID;
        r->positions.swap(copy);
        return PAACS_UI_OK;
    } catch (...) {
        return PAACS_UI_ALLOCATION_FAILED;
    }
}
paacs_ui_status paacs_ui_request_profile_controls(paacs_ui_request *r,
                                                  const paacs_ui_profile_control *controls,
                                                  size_t n) {
    if (!r || !r->primitive || (n && !controls))
        return PAACS_UI_INVALID;
    if (n > 1024)
        return PAACS_UI_LIMIT;
    try {
        std::vector<sizing::ProfileSizeControl> copy;
        copy.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            const auto &v = controls[i];
            if (!std::isfinite(v.begin) || !std::isfinite(v.end) || !std::isfinite(v.size) ||
                !std::isfinite(v.growth))
                return PAACS_UI_INVALID;
            copy.push_back({v.begin, v.end, v.size, v.growth});
        }
        r->profile.controls.swap(copy);
        return PAACS_UI_OK;
    } catch (...) {
        return PAACS_UI_ALLOCATION_FAILED;
    }
}
paacs_ui_status paacs_ui_request_field(paacs_ui_request *r, const paacs_ui_field_node *nodes,
                                       size_t n, uint32_t output, double minimum, double maximum,
                                       double factor) {
    if (!r || r->primitive || (n && !nodes))
        return PAACS_UI_INVALID;
    if (n > sizing::maximumControlFieldNodes)
        return PAACS_UI_LIMIT;
    if (!n) {
        r->field.reset();
        return PAACS_UI_OK;
    }
    try {
        sizing::ControlFieldSpecification spec;
        spec.output = output;
        spec.minimumSize = minimum;
        spec.maximumSize = maximum;
        spec.sizeFactor = factor;
        spec.nodes.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            if (nodes[i].operation > static_cast<unsigned>(sizing::FieldOperation::Clamp))
                return PAACS_UI_INVALID;
            sizing::FieldNode v;
            v.operation = static_cast<sizing::FieldOperation>(nodes[i].operation);
            v.first = nodes[i].first;
            v.second = nodes[i].second;
            v.a = vec(nodes[i].a);
            v.b = vec(nodes[i].b);
            v.c = vec(nodes[i].c);
            std::copy_n(nodes[i].parameters, 4, v.parameters.begin());
            spec.nodes.push_back(v);
        }
        // All semantic field checks run at submit, returning a named diagnostic.
        r->field = std::move(spec);
        return PAACS_UI_OK;
    } catch (...) {
        return PAACS_UI_ALLOCATION_FAILED;
    }
}
paacs_ui_status paacs_ui_submit(paacs_ui_service *s, const paacs_ui_request *r, paacs_ui_job **out,
                                char *error, size_t capacity) {
    if (out)
        *out = nullptr;
    diagnostic(error, capacity, "");
    if (!s || !r || !out)
        return PAACS_UI_INVALID;
    try {
        sdk::MeshTask task;
        if (r->primitive) {
            auto options = r->p;
            options.axialPositions = r->positions;
            sizing::CappedSizeOutcome profile;
            if (r->profileEnabled) {
                profile = sizing::compileCappedSizeProfile(r->profile);
                if (!profile.compiled) {
                    diagnostic(error, capacity, sizing::profileReasonName(profile.reason));
                    return PAACS_UI_INVALID;
                }
                options.sizeProfile = &profile.profile;
            }
            task = sdk::primitiveTask(r->shape, options);
        } else {
            std::shared_ptr<const sizing::BoundedSizeField> field;
            if (r->field) {
                auto built = sizing::compileControlField(*r->field);
                if (!built.compiled) {
                    diagnostic(error, capacity, sizing::controlFieldReasonName(built.reason));
                    return PAACS_UI_INVALID;
                }
                field = std::make_shared<sizing::ControlField>(std::move(built.field));
            }
            task = sdk::adaptiveBoxTask(r->a, std::move(field));
        }
        auto handle = std::make_unique<paacs_ui_job>();
        auto submitted = s->service.submit(std::move(task));
        if (!submitted.job) {
            diagnostic(error, capacity,
                       submitted.reason == sdk::SubmitReason::JobLimit
                           ? "job handle limit reached; release completed jobs"
                           : "service refused job");
            return submitted.reason == sdk::SubmitReason::AllocationFailed
                       ? PAACS_UI_ALLOCATION_FAILED
                       : PAACS_UI_LIMIT;
        }
        handle->job = std::move(submitted.job);
        *out = handle.release();
        return PAACS_UI_OK;
    } catch (const std::bad_alloc &) {
        diagnostic(error, capacity, "allocation failed");
        return PAACS_UI_ALLOCATION_FAILED;
    } catch (...) {
        diagnostic(error, capacity, "request snapshot failed");
        return PAACS_UI_INVALID;
    }
}
void paacs_ui_job_release(paacs_ui_job *j) {
    delete j;
}
paacs_ui_status paacs_ui_job_snapshot(const paacs_ui_job *j, paacs_ui_snapshot *out) {
    if (!j || !out)
        return PAACS_UI_INVALID;
    try {
        auto s = j->job.snapshot();
        paacs_ui_snapshot value{};
        value.id = s.id;
        value.completed = s.completed;
        value.total = s.total;
        value.elapsed_seconds = s.elapsedSeconds;
        value.state = static_cast<uint32_t>(s.state);
        value.stage = static_cast<uint32_t>(s.stage);
        value.cancellation_requested = s.cancellationRequested;
        value.error = static_cast<uint32_t>(s.error);
        value.domain = static_cast<uint32_t>(s.domain);
        value.engine_reason = s.engineReason;
        std::copy(s.diagnostic.begin(), s.diagnostic.end(), value.diagnostic);
        *out = value;
        return PAACS_UI_OK;
    } catch (...) {
        return PAACS_UI_FAILED;
    }
}
paacs_ui_status paacs_ui_job_cancel(paacs_ui_job *j, uint32_t *out) {
    if (!j || !out)
        return PAACS_UI_INVALID;
    try {
        *out = j->job.cancel();
        return PAACS_UI_OK;
    } catch (...) {
        return PAACS_UI_FAILED;
    }
}
paacs_ui_status paacs_ui_job_wait(const paacs_ui_job *j, uint32_t ms, uint32_t *out) {
    if (!j || !out)
        return PAACS_UI_INVALID;
    try {
        *out = j->job.waitFor(std::chrono::milliseconds(ms));
        return PAACS_UI_OK;
    } catch (...) {
        return PAACS_UI_FAILED;
    }
}
paacs_ui_status paacs_ui_job_result(const paacs_ui_job *j, paacs_ui_result **out) {
    if (out)
        *out = nullptr;
    if (!j || !out)
        return PAACS_UI_INVALID;
    try {
        auto r = j->job.result();
        if (!r)
            return PAACS_UI_NOT_READY;
        *out = new paacs_ui_result{std::move(r)};
        return PAACS_UI_OK;
    } catch (...) {
        return PAACS_UI_ALLOCATION_FAILED;
    }
}
void paacs_ui_result_release(paacs_ui_result *r) {
    delete r;
}
paacs_ui_status paacs_ui_result_summary(const paacs_ui_result *r, paacs_ui_summary *out) {
    if (!r || !out)
        return PAACS_UI_INVALID;
    const auto &s = r->result->summary;
    paacs_ui_summary value{};
    value.nodes = s.nodes;
    value.cells = s.cells;
    value.connectivity = s.connectivity;
    value.logical_bytes = s.logicalBytes;
    value.hash = s.hash;
    value.hash_available = s.hashAvailable;
    std::copy(s.cellTypes.begin(), s.cellTypes.end(), value.cell_types);
    vec(value.lower, s.lower);
    vec(value.upper, s.upper);
    if (auto *p = std::get_if<algo::PrimitiveMeshingResult>(&r->result->report)) {
        value.minimum_quality = p->minimumQuality;
        value.mesh_volume = p->meshVolume;
        value.quality_warning = p->qualityWarning;
    } else if (auto *a = std::get_if<algo::AdaptiveBoxResult>(&r->result->report)) {
        value.minimum_quality = a->minimumQuality;
        value.mesh_volume = a->meshVolume;
    }
    *out = value;
    return PAACS_UI_OK;
}
paacs_ui_status paacs_ui_result_nodes(const paacs_ui_result *r, uint64_t offset, paacs_ui_node *out,
                                      size_t capacity, size_t *written) {
    if (written)
        *written = 0;
    if (!r || !written || (capacity && !out))
        return PAACS_UI_INVALID;
    if (capacity > PAACS_UI_MAX_PAGE)
        return PAACS_UI_LIMIT;
    const auto &m = r->result->mesh;
    if (offset > static_cast<uint64_t>(m.nodeCount()))
        return PAACS_UI_OUT_OF_RANGE;
    auto n = static_cast<size_t>(
        std::min<uint64_t>(capacity, static_cast<uint64_t>(m.nodeCount()) - offset));
    for (size_t i = 0; i < n; ++i) {
        auto id = static_cast<Index>(offset + i);
        vec(out[i].position, m.node(id));
        out[i].dimension = m.nodeClass(id).dim;
        out[i].tag = m.nodeClass(id).tag;
    }
    *written = n;
    return PAACS_UI_OK;
}
paacs_ui_status paacs_ui_result_cells(const paacs_ui_result *r, uint64_t offset, paacs_ui_cell *out,
                                      size_t capacity, size_t *written) {
    if (written)
        *written = 0;
    if (!r || !written || (capacity && !out))
        return PAACS_UI_INVALID;
    if (capacity > PAACS_UI_MAX_PAGE)
        return PAACS_UI_LIMIT;
    const auto &m = r->result->mesh;
    if (offset > static_cast<uint64_t>(m.cellCount()))
        return PAACS_UI_OUT_OF_RANGE;
    auto n = static_cast<size_t>(
        std::min<uint64_t>(capacity, static_cast<uint64_t>(m.cellCount()) - offset));
    for (size_t i = 0; i < n; ++i) {
        auto id = static_cast<Index>(offset + i);
        const auto &cell = m.cell(id);
        paacs_ui_cell value{};
        value.type = static_cast<uint32_t>(cell.type);
        value.node_count = static_cast<uint32_t>(mesh::nodesPerCell(cell.type));
        value.dimension = cell.onEntity.dim;
        value.tag = cell.onEntity.tag;
        for (uint32_t k = 0; k < value.node_count; ++k)
            value.nodes[k] = static_cast<uint64_t>(m.cellNodes(id)[k]);
        out[i] = value;
    }
    *written = n;
    return PAACS_UI_OK;
}
paacs_ui_status paacs_ui_result_write_archive(const paacs_ui_result *r, const char *path) {
    if (!r || !path || !*path)
        return PAACS_UI_INVALID;
    try {
        return sdk::writeMeshArchive(r->result->mesh, path).completed ? PAACS_UI_OK
                                                                      : PAACS_UI_FAILED;
    } catch (...) {
        return PAACS_UI_FAILED;
    }
}
}
