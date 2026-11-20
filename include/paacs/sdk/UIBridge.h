// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
/* Optional shared-library adapter. C11 ABI v1, independent of the C++ ABI.
 * Requests are thread-confined. Job polling/cancellation and immutable results
 * are safe for concurrent reads. Destroy/release must not race any use of that
 * same C handle. The native C++ JobService exposes additional routes.
 */
#include <stddef.h>
#include <stdint.h>
#if defined(_WIN32)
#if defined(PAACS_UI_BUILD)
#define PAACS_UI_API __declspec(dllexport)
#else
#define PAACS_UI_API __declspec(dllimport)
#endif
#else
#define PAACS_UI_API __attribute__((visibility("default")))
#endif
#ifdef __cplusplus
extern "C" {
#endif
#define PAACS_UI_ABI_VERSION 1u
#define PAACS_UI_MAX_PAGE 65536u
#define PAACS_UI_COMPUTE_HASH 1u
typedef struct paacs_ui_service paacs_ui_service;
typedef struct paacs_ui_request paacs_ui_request;
typedef struct paacs_ui_job paacs_ui_job;
typedef struct paacs_ui_result paacs_ui_result;
typedef enum paacs_ui_status {
    PAACS_UI_OK = 0,
    PAACS_UI_INVALID = 1,
    PAACS_UI_FAILED = 2,
    PAACS_UI_ALLOCATION_FAILED = 3,
    PAACS_UI_OUT_OF_RANGE = 4,
    PAACS_UI_NOT_READY = 5,
    PAACS_UI_LIMIT = 6
} paacs_ui_status;
typedef enum paacs_ui_job_state {
    PAACS_UI_JOB_QUEUED = 0,
    PAACS_UI_JOB_WAITING_FOR_ENGINE = 1,
    PAACS_UI_JOB_RUNNING = 2,
    PAACS_UI_JOB_SUCCEEDED = 3,
    PAACS_UI_JOB_FAILED = 4,
    PAACS_UI_JOB_CANCELLED = 5
} paacs_ui_job_state;
typedef enum paacs_ui_job_stage {
    PAACS_UI_STAGE_QUEUED = 0,
    PAACS_UI_STAGE_PREPARING = 1,
    PAACS_UI_STAGE_GENERATING = 2,
    PAACS_UI_STAGE_FINALIZING = 3,
    PAACS_UI_STAGE_COMPLETE = 4
} paacs_ui_job_stage;
typedef enum paacs_ui_job_error {
    PAACS_UI_ERROR_NONE = 0,
    PAACS_UI_ERROR_ENGINE_REJECTED = 1,
    PAACS_UI_ERROR_CANCELLED = 2,
    PAACS_UI_ERROR_RESULT_BUDGET = 3,
    PAACS_UI_ERROR_ALLOCATION = 4,
    PAACS_UI_ERROR_EXCEPTION = 5
} paacs_ui_job_error;
typedef enum paacs_ui_job_domain {
    PAACS_UI_DOMAIN_SERVICE = 0,
    PAACS_UI_DOMAIN_PRIMITIVE = 1,
    PAACS_UI_DOMAIN_ADAPTIVE_BOX = 2,
    PAACS_UI_DOMAIN_GEOMETRY = 3,
    PAACS_UI_DOMAIN_LAYERED_PATCH = 4,
    PAACS_UI_DOMAIN_MAPPED_HEX = 5,
    PAACS_UI_DOMAIN_CUSTOM = 6
} paacs_ui_job_domain;
typedef struct paacs_ui_option {
    const char *key; /* Strings have library lifetime. */
    uint32_t kind;   /* 0 real, 1 integer, 2 boolean. */
    double default_value, minimum, maximum;
    const char *choices; /* Comma-separated enum names, ordinal from zero. */
} paacs_ui_option;
typedef struct paacs_ui_snapshot {
    uint64_t id, completed, total;
    double elapsed_seconds;
    uint32_t state, stage, cancellation_requested, error, domain, engine_reason;
    char diagnostic[192];
} paacs_ui_snapshot;
typedef struct paacs_ui_summary {
    uint64_t nodes, cells, connectivity, logical_bytes, hash, cell_types[10];
    double lower[3], upper[3], minimum_quality, mesh_volume;
    uint32_t quality_warning, hash_available;
} paacs_ui_summary;
typedef struct paacs_ui_node {
    double position[3];
    int32_t dimension, tag;
} paacs_ui_node;
typedef struct paacs_ui_cell {
    uint32_t type, node_count;
    int32_t dimension, tag;
    uint64_t nodes[10];
} paacs_ui_cell;
typedef struct paacs_ui_profile_control {
    double begin, end, size, growth;
} paacs_ui_profile_control;
/* Operations are FieldOperation ordinals: constant, distance_point,
 * distance_segment, distance_triangle, distance_box, distance_ball,
 * threshold_linear, threshold_smooth, minimum, maximum, scale, clamp.
 * Dependencies use zero-based indices; UINT32_MAX means unused.
 */
typedef struct paacs_ui_field_node {
    uint32_t operation, first, second;
    double a[3], b[3], c[3], parameters[4];
} paacs_ui_field_node;
PAACS_UI_API uint32_t paacs_ui_abi_version(void);
PAACS_UI_API const char *paacs_ui_release_version(void);
/* sizeof values for bindings: option=0, snapshot=1, summary=2, node=3,
 * cell=4, profile_control=5, field_node=6. Unknown returns zero. */
PAACS_UI_API size_t paacs_ui_abi_size(uint32_t type);
PAACS_UI_API paacs_ui_service *paacs_ui_service_create(uint32_t abi, uint32_t maximum_jobs,
                                                       uint64_t maximum_result_bytes,
                                                       uint64_t maximum_retained_bytes);
/* Default service skips strict hashing. Unknown flag bits reject. */
PAACS_UI_API paacs_ui_service *paacs_ui_service_create_with_flags(uint32_t abi,
                                                                  uint32_t maximum_jobs,
                                                                  uint64_t maximum_result_bytes,
                                                                  uint64_t maximum_retained_bytes,
                                                                  uint32_t flags);
/* Cancels and joins work. Cooperative shutdown may block during an engine phase. */
PAACS_UI_API void paacs_ui_service_destroy(paacs_ui_service *);
PAACS_UI_API uint64_t paacs_ui_service_retained_bytes(const paacs_ui_service *);
/* Routes in ABI v1: "primitive", "adaptive_box". Unknown/version mismatch -> NULL. */
PAACS_UI_API paacs_ui_request *paacs_ui_request_create(uint32_t abi, const char *route);
PAACS_UI_API void paacs_ui_request_destroy(paacs_ui_request *);
PAACS_UI_API size_t paacs_ui_option_count(const paacs_ui_request *);
PAACS_UI_API paacs_ui_status paacs_ui_option_at(const paacs_ui_request *, size_t index,
                                                paacs_ui_option *);
PAACS_UI_API paacs_ui_status paacs_ui_request_get(const paacs_ui_request *, const char *key,
                                                  double *);
/* Unknown keys, nonfinite/fractional integers and out-of-range values reject.
 * Setters are transactional. Bounds are representation limits; cross-option
 * and geometric validation is performed by the unchanged engine on execution. */
PAACS_UI_API paacs_ui_status paacs_ui_request_set(paacs_ui_request *, const char *key, double);
PAACS_UI_API paacs_ui_status paacs_ui_request_axial_positions(paacs_ui_request *, const double *,
                                                              size_t);
/* Stores controls; enable them with request_set("profile_enabled", 1). */
PAACS_UI_API paacs_ui_status paacs_ui_request_profile_controls(paacs_ui_request *,
                                                               const paacs_ui_profile_control *,
                                                               size_t);
PAACS_UI_API paacs_ui_status paacs_ui_request_field(paacs_ui_request *, const paacs_ui_field_node *,
                                                    size_t, uint32_t output, double minimum_size,
                                                    double maximum_size, double size_factor);
/* Deep snapshot before returning. Caller may edit/destroy request immediately.
 * *job becomes NULL on failure; diagnostic is bounded, NUL terminated. */
PAACS_UI_API paacs_ui_status paacs_ui_submit(paacs_ui_service *, const paacs_ui_request *,
                                             paacs_ui_job **job, char *diagnostic,
                                             size_t diagnostic_capacity);
PAACS_UI_API void paacs_ui_job_release(paacs_ui_job *);
PAACS_UI_API paacs_ui_status paacs_ui_job_snapshot(const paacs_ui_job *, paacs_ui_snapshot *);
/* accepted=1 guarantees the job cannot subsequently publish; late calls -> 0. */
PAACS_UI_API paacs_ui_status paacs_ui_job_cancel(paacs_ui_job *, uint32_t *accepted);
PAACS_UI_API paacs_ui_status paacs_ui_job_wait(const paacs_ui_job *, uint32_t milliseconds,
                                               uint32_t *terminal);
PAACS_UI_API paacs_ui_status paacs_ui_job_result(const paacs_ui_job *, paacs_ui_result **);
/* Result retains the mesh, including after job/service destruction. No mesh copy. */
PAACS_UI_API void paacs_ui_result_release(paacs_ui_result *);
PAACS_UI_API paacs_ui_status paacs_ui_result_summary(const paacs_ui_result *, paacs_ui_summary *);
/* At most PAACS_UI_MAX_PAGE per call. offset==count returns an empty page;
 * offset>count rejects. Cell ordinals are offset+i. Connectivity uses global
 * node indices; type ordinals match mesh::CellType, tags remain unchanged. */
PAACS_UI_API paacs_ui_status paacs_ui_result_nodes(const paacs_ui_result *, uint64_t offset,
                                                   paacs_ui_node *, size_t capacity,
                                                   size_t *written);
PAACS_UI_API paacs_ui_status paacs_ui_result_cells(const paacs_ui_result *, uint64_t offset,
                                                   paacs_ui_cell *, size_t capacity,
                                                   size_t *written);
PAACS_UI_API paacs_ui_status paacs_ui_result_write_archive(const paacs_ui_result *,
                                                           const char *path);
#ifdef __cplusplus
}
#endif
