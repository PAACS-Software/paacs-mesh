// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include <paacs/sdk/CAPI.h>
#include <stdio.h>
int main(void) {
    paacs_v8_session* s=paacs_v8_session_create();
    if(!s)return 1;
    paacs_v8_hex_request r={0};
    r.axis_u.x=1;r.axis_v.y=1;r.axis_w.z=1;
    r.cells_u=2;r.cells_v=2;r.cells_w=2;
    r.volume_tag=1;r.boundary_tag_base=101;r.threads=1;
    const int ok=paacs_v8_generate_mapped_hex(s,&r)==PAACS_V8_OK &&
                 paacs_v8_session_node_count(s)==27 && paacs_v8_session_cell_count(s)>0;
    printf("%s C consumer %s\n",paacs_v8_api_version(),ok?"passed":"failed");
    paacs_v8_session_destroy(s);
    return ok?0:1;
}
