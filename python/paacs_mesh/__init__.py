"""PAACS Mesh: native STEP import, audited Tet4 meshing and tagged export.

STEP coordinates and mesh sizes use millimetres. ``Mesh.report`` contains the
native audit, quality, timings and assembly interference results. The default
assembly policy preserves independently validated parts; overlapping parts do
not constitute a globally non-overlapping material domain.

Native import and mesh calls are serialized across models in this initial API.
The active mesh job may use the requested internal workers. Its cooperative
``max_seconds`` deadline starts after it acquires the native operation gate.

``FaceSize`` and ``BoxSize`` are read-only sizing controls. A control's unset
growth inherits ``mesh(growth=...)``. Face and box controls select the native
spatial sizing route by default, which currently requires one solid and an
explicit ``assembly="conforming"`` request. ``native_sizing=False`` retains the
legacy face-cap route; volume boxes require native spatial sizing.

``Mesh.to_numpy()`` makes independent, initially read-only coordinate,
connectivity and CAD-classification arrays. NumPy is optional and loaded only
for extraction; import, meshing and export work without it.
"""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import os
from pathlib import Path


def _configure_bundled_resources():
    resources = Path(__file__).resolve().parent / "_occt" / "resources"
    for variable, directory in (
        ("CSF_XSMessage", "XSMessage"),
        ("CSF_SHMessage", "SHMessage"),
        ("CSF_STEPDefaults", "XSTEPResource"),
        ("CSF_IGESDefaults", "XSTEPResource"),
    ):
        target = resources / directory
        if target.is_dir():
            os.environ.setdefault(variable, str(target))


_configure_bundled_resources()

from ._native import (  # noqa: E402: configure OCCT resources before loading it
    BoxSize, CadImportError, CadModel, ExportError, FaceSize, Mesh, MeshingError,
    engine_version as native_version, read_step,
)

__all__ = [
    "BoxSize", "CadImportError", "CadModel", "ExportError", "FaceSize", "Mesh", "MeshingError", "read_step",
    "native_version",
]

try:
    from importlib.metadata import PackageNotFoundError, version

    __version__ = version("paacs-mesh")
except PackageNotFoundError:
    __version__ = native_version.replace("-rc.", "rc")
