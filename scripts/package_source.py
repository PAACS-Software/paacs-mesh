#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Create a clean deterministic source ZIP and SHA256 manifest (no build output)."""
import argparse
import hashlib
import os
from pathlib import Path
import zipfile
from verify_version import verify_version

parser = argparse.ArgumentParser()
parser.add_argument('--output', type=Path)
parser.add_argument('--include-local-evidence', action='store_true',
                    help='Include local evidence and reproducers for an internal archive')
parser.add_argument('--preserve-source-manifest', action='store_true',
                    help='Generate the archive manifest without replacing the source SHA256SUMS.txt')
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
version = verify_version(root)
name = f'paacs-mesh-{version}'
output = (args.output or root.parent / f'{name}.zip').resolve()
output.parent.mkdir(parents=True, exist_ok=True)
checksum_path = output.with_name(output.name + '.sha256')
source_roots = {
    '.github', 'apps', 'benchmarks', 'cmake', 'docs', 'examples', 'include',
    'regression_corpus', 'scripts', 'src', 'standalone', 'tests',
}
root_files = {
    '.clang-format', '.gitattributes', '.gitignore', 'CMakeLists.txt',
    'CMakePresets.json', 'CHANGELOG.md', 'CITATION.cff', 'CONTRIBUTING.md', 'LICENSING_STATUS.md',
    'README.md', 'RELEASE_NOTES.md', 'VALIDATION.md', 'VERSION.txt',
}
excluded = {'.git', '.vs', '__pycache__'}
source_suffixes = {
    '.brep', '.c', '.cc', '.cmake', '.cpp', '.csv', '.cxx', '.h', '.hpp', '.hxx',
    '.html', '.in', '.json', '.md', '.png', '.ps1', '.py', '.sh', '.step',
    '.txt', '.yaml', '.yml',
}
excluded_suffixes = {
    '.zip', '.paacsmesh', '.pyc', '.exe', '.obj', '.dll', '.lib', '.pdb',
    '.ilk', '.exp', '.pch', '.idb', '.ipch', '.o', '.a', '.so', '.dylib',
}
local_evidence_suffixes = {
    '.md', '.txt', '.json', '.jsonl', '.csv', '.tsv', '.log', '.xml',
    '.cpp', '.hpp', '.h', '.cc', '.cxx', '.c', '.py', '.ps1', '.cmd',
    '.bat', '.cmake', '.sh', '.toml', '.yaml', '.yml',
}


def include_local_evidence(path):
    parts = path.relative_to(root).parts
    if parts[0] == 'evidence':
        return args.include_local_evidence and path.suffix.lower() in local_evidence_suffixes
    return True


if args.include_local_evidence:
    source_roots.add('evidence')


def source_files():
    for path in root.iterdir():
        if (path.is_file() and not path.is_symlink() and
                (path.name in root_files or path.stem in {'LICENSE', 'NOTICE'})):
            yield path
    for directory, subdirectories, filenames in os.walk(root):
        if Path(directory) == root:
            subdirectories[:] = [name for name in subdirectories if name in source_roots]
            continue
        subdirectories[:] = [name for name in subdirectories
                             if name not in excluded and not name.startswith('.')]
        for name in filenames:
            path = Path(directory) / name
            suffix = path.suffix.lower()
            is_evidence = path.relative_to(root).parts[0] == 'evidence'
            if (path.is_file() and not path.is_symlink() and not name.startswith('.')
                and name not in {'SHA256SUMS.txt', 'CMakeUserPresets.json'}
                and suffix not in excluded_suffixes
                and (is_evidence or suffix in source_suffixes)
                and '.so.' not in name.lower()
                and include_local_evidence(path)
                and path.resolve() not in {output, checksum_path}):
                yield path


def source_bytes(path):
    data = path.read_bytes()
    # Match Git's canonical text form on every build host. PNGs stay binary.
    return data if path.suffix.lower() == '.png' else data.replace(b'\r\n', b'\n')


files = sorted(source_files())
# Read each input once so a local log changing during packaging cannot corrupt
# the relationship between the manifest and the archived content.
contents = {p.relative_to(root).as_posix(): source_bytes(p) for p in files}
manifest = ''.join(f'{hashlib.sha256(data).hexdigest()}  {filename}\n'
                   for filename, data in contents.items())
if not args.preserve_source_manifest:
    (root / 'SHA256SUMS.txt').write_text(manifest)
contents['SHA256SUMS.txt'] = manifest.encode('utf-8')
with zipfile.ZipFile(output, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=9) as z:
    for filename, data in sorted(contents.items()):
        info = zipfile.ZipInfo(f'{name}/{filename}', (2026, 10, 2, 0, 0, 0))
        info.compress_type = zipfile.ZIP_DEFLATED
        info.external_attr = (0o100755 if Path(filename).suffix in {'.sh', '.py'} else 0o100644) << 16
        z.writestr(info, data)
with zipfile.ZipFile(output) as z:
    if z.testzip() is not None:
        raise RuntimeError('ZIP CRC validation failed')
    for line in manifest.splitlines():
        expected, filename = line.split('  ', 1)
        if hashlib.sha256(z.read(f'{name}/{filename}')).hexdigest() != expected:
            raise RuntimeError(f'Checksum mismatch: {filename}')
print(f'{output}: {len(contents)} files, {output.stat().st_size} bytes')
digest = hashlib.sha256(output.read_bytes()).hexdigest()
checksum_path.write_text(f'{digest}  {output.name}\n')
print('SHA256:', digest)
print('Checksum:', checksum_path)
