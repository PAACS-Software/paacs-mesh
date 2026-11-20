#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Reject mismatched source/build/SDK version metadata before packaging."""
import argparse
import json
from pathlib import Path
import re


def single(pattern, text, description):
    matches = re.findall(pattern, text, flags=re.MULTILINE)
    if len(matches) != 1:
        raise ValueError(f'Expected one {description}, found {len(matches)}')
    return matches[0]


def verify_version(root):
    version = (root / 'VERSION.txt').read_text().strip()
    match = re.fullmatch(r'(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)(?:-[0-9A-Za-z.-]+)?', version)
    if match is None:
        raise ValueError('VERSION.txt must contain a numeric version with an optional prerelease suffix')
    parts = match.groups()
    base = '.'.join(parts)
    cmake = (root / 'CMakeLists.txt').read_text()
    sdk = (root / 'include/paacs/sdk/Version.hpp').read_text()
    values = {
        'CMake project': (single(r'project\(PAACSMesh\s+VERSION\s+(\S+)\s+LANGUAGES', cmake, 'project version'), base),
        'CMake release': (single(r'^set\(PAACS_RELEASE_VERSION\s+"([^"]+)"\)', cmake, 'release version'), version),
        'SDK version': (single(r'kApiVersion\s*=\s*"([^"]+)"', sdk, 'SDK string'), version),
    }
    citation = root / 'CITATION.cff'
    if citation.is_file():
        citation_version = single(
            r'''^version:\s*['"]?([0-9A-Za-z.-]+)['"]?\s*(?:#.*)?$''',
            citation.read_text(encoding='utf-8-sig'), 'citation version')
        values['Citation version'] = (citation_version, version)
    for label, expected in zip(['Major', 'Minor', 'Patch'], parts):
        values['SDK ' + label] = (single(r'kApi' + label + r'\s*=\s*([0-9]+)', sdk, 'SDK ' + label), expected)
    for label, (actual, expected) in values.items():
        if actual != expected:
            raise ValueError(f'{label}: {actual!r} does not match {expected!r}')
    return version


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--final', action='store_true', help='also reject a prerelease suffix')
    args = parser.parse_args()
    try:
        version = verify_version(args.root.resolve())
        if args.final and '-' in version:
            raise ValueError('A final release cannot have a prerelease suffix')
    except (OSError, ValueError) as error:
        parser.exit(1, str(error) + '\n')
    print(json.dumps(dict(version=version, consistent=True, final='-' not in version)))
