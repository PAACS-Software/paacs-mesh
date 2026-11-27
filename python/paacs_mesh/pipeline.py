# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0
"""Explicit CAD boundary → classical initial sites → Tet4 construction stages.

The supported 1.2.1 scope is uniform sizing on one original solid. Preparation
qualifies a frozen boundary; point generation stops before constraint recovery
or volume construction. Only the final mesh result can pass volume publication
checks. No learned model or training dependency is part of this module.

``paacs_mesh.experimental`` retains the same class/function identities for
compatibility with earlier development scripts.
"""
from .experimental import (
    InitialPointCloud, PreparedBoundary, initial_points, load_boundary,
    load_initial_points, prepare_boundary,
)

__all__ = [
    "InitialPointCloud", "PreparedBoundary", "initial_points", "load_boundary",
    "load_initial_points", "prepare_boundary",
]
