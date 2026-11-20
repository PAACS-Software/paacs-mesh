// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// PAACS 1.0 retains the existing V8/V0.10 C symbol names for source compatibility.
// All ownership stays inside an opaque session. Returned strings remain valid
// until the next mutating call on that session or session destruction.
typedef struct paacs_v8_session paacs_v8_session;

typedef struct paacs_v8_vec3 {
    double x;
    double y;
    double z;
} paacs_v8_vec3;

typedef struct paacs_v8_classification {
    int32_t dimension;
    int32_t tag;
} paacs_v8_classification;

typedef enum paacs_v8_cell_type {
    PAACS_V8_CELL_LINE2 = 0,
    PAACS_V8_CELL_LINE3 = 1,
    PAACS_V8_CELL_TRI3 = 2,
    PAACS_V8_CELL_TRI6 = 3,
    PAACS_V8_CELL_QUAD4 = 4,
    PAACS_V8_CELL_TET4 = 5,
    PAACS_V8_CELL_TET10 = 6,
    PAACS_V8_CELL_PRISM6 = 7,
    PAACS_V8_CELL_HEX8 = 8,
    PAACS_V8_CELL_PYRAMID5 = 9
} paacs_v8_cell_type;

// Fixed-width read-only cell view.  node_count states how many entries of
// nodes[] are meaningful; ten covers every standard V8 element family.
typedef struct paacs_v8_cell_view {
    paacs_v8_cell_type type;
    uint32_t node_count;
    paacs_v8_classification classification;
    uint64_t nodes[10];
} paacs_v8_cell_view;

typedef struct paacs_v8_hex_request {
    paacs_v8_vec3 origin;
    paacs_v8_vec3 axis_u;
    paacs_v8_vec3 axis_v;
    paacs_v8_vec3 axis_w;
    uint32_t cells_u;
    uint32_t cells_v;
    uint32_t cells_w;
    int32_t volume_tag;
    int32_t boundary_tag_base;
    int32_t threads;
} paacs_v8_hex_request;

typedef enum paacs_v8_status {
    PAACS_V8_OK = 0,
    PAACS_V8_INVALID_ARGUMENT = 1,
    PAACS_V8_FAILED = 2,
    PAACS_V8_ALLOCATION_FAILED = 3,
    PAACS_V8_OUT_OF_RANGE = 4
} paacs_v8_status;

const char* paacs_v8_api_version(void);
paacs_v8_session* paacs_v8_session_create(void);
void paacs_v8_session_destroy(paacs_v8_session* session);

paacs_v8_status paacs_v8_generate_mapped_hex(
    paacs_v8_session* session,
    const paacs_v8_hex_request* request);

// V0.10 prepared-production Tet request.  The primitive helpers exercise the
// same geometry-independent pipeline that an OCCT-backed Model will use.
typedef enum paacs_v010_tet_output {
    PAACS_V010_TET4_CONFORMING = 0,
    PAACS_V010_TET4_REFINED = 1,
    PAACS_V010_TET4_QUALITY = 2,
    PAACS_V010_TET10_CERTIFIED = 3
} paacs_v010_tet_output;

typedef struct paacs_v010_tet_options {
    double target_size;
    paacs_v010_tet_output output;
    uint64_t maximum_insertions;
    uint64_t maximum_nodes;
    uint64_t maximum_tets;
    int32_t threads;
} paacs_v010_tet_options;

typedef struct paacs_v010_box_request {
    paacs_v8_vec3 lower;
    paacs_v8_vec3 upper;
    paacs_v010_tet_options options;
} paacs_v010_box_request;

typedef struct paacs_v010_cylinder_request {
    paacs_v8_vec3 base_center;
    double radius;
    double height;
    paacs_v010_tet_options options;
} paacs_v010_cylinder_request;

paacs_v8_status paacs_v010_generate_box_tet(
    paacs_v8_session* session,
    const paacs_v010_box_request* request);
paacs_v8_status paacs_v010_generate_cylinder_tet(
    paacs_v8_session* session,
    const paacs_v010_cylinder_request* request);

paacs_v8_status paacs_v8_session_read_archive(
    paacs_v8_session* session,
    const char* path);
paacs_v8_status paacs_v8_session_write_archive(
    const paacs_v8_session* session,
    const char* path);

uint64_t paacs_v8_session_node_count(const paacs_v8_session* session);
uint64_t paacs_v8_session_cell_count(const paacs_v8_session* session);
uint64_t paacs_v8_session_mesh_hash(const paacs_v8_session* session);

paacs_v8_status paacs_v8_session_node_at(
    const paacs_v8_session* session,
    uint64_t index,
    paacs_v8_vec3* position,
    paacs_v8_classification* classification);
paacs_v8_status paacs_v8_session_cell_at(
    const paacs_v8_session* session,
    uint64_t index,
    paacs_v8_cell_view* cell);

const char* paacs_v8_session_last_error(const paacs_v8_session* session);

#ifdef __cplusplus
}
#endif
