// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/sdk/UIBridge.h"
#include <stdio.h>
int main(void) {
    paacs_ui_service *s =
        paacs_ui_service_create(PAACS_UI_ABI_VERSION, 4, 16 * 1024 * 1024, 32 * 1024 * 1024);
    paacs_ui_request *r = paacs_ui_request_create(PAACS_UI_ABI_VERSION, "adaptive_box");
    paacs_ui_job *j = NULL;
    paacs_ui_result *result = NULL;
    char error[192];
    uint32_t done = 0;
    paacs_ui_summary summary;
    int failed = 1;
    if (!s || !r)
        goto cleanup;
    if (paacs_ui_request_set(r, "target_size", .5) != PAACS_UI_OK)
        goto cleanup;
    if (paacs_ui_submit(s, r, &j, error, sizeof(error)) != PAACS_UI_OK)
        goto cleanup;
    if (paacs_ui_job_wait(j, 10000, &done) != PAACS_UI_OK || !done)
        goto cleanup;
    if (paacs_ui_job_result(j, &result) != PAACS_UI_OK)
        goto cleanup;
    if (paacs_ui_result_summary(result, &summary) != PAACS_UI_OK || !summary.nodes ||
        !summary.cells)
        goto cleanup;
    if (summary.hash_available)
        goto cleanup;
    printf("C UI ABI %u: %llu nodes\n", paacs_ui_abi_version(), (unsigned long long)summary.nodes);
    failed = 0;
cleanup:
    paacs_ui_result_release(result);
    paacs_ui_job_release(j);
    paacs_ui_request_destroy(r);
    paacs_ui_service_destroy(s);
    return failed;
}
