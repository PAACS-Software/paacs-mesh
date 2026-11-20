# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

param(
    [Parameter(Mandatory = $true)][string]$Archive,
    [string]$OutputRoot,
    [switch]$ActivateRuntime
)
$ErrorActionPreference = 'Stop'
$workspaceRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
if (!$OutputRoot) { $OutputRoot = Join-Path $workspaceRoot 'build\deps\occt-8.0.1' }
$targetRoot = [IO.Path]::GetFullPath($OutputRoot)
if (!$targetRoot.StartsWith($workspaceRoot + '\', [StringComparison]::OrdinalIgnoreCase)) {
    throw 'The SDK extraction target must be inside this workspace.'
}
Add-Type -AssemblyName System.IO.Compression.FileSystem
$outer = [IO.Compression.ZipFile]::OpenRead([IO.Path]::GetFullPath($Archive))
$temporaryArchive = $null
$sdkArchive = $outer
try {
    # The downloaded combined package may be wrapped in another ZIP.
    if ($outer.Entries.Count -eq 1 -and $outer.Entries[0].FullName.EndsWith('.zip')) {
        $temporaryArchive = Join-Path ([IO.Path]::GetTempPath()) ('paacs-occt-' + [guid]::NewGuid() + '.zip')
        [IO.Compression.ZipFileExtensions]::ExtractToFile($outer.Entries[0], $temporaryArchive, $false)
        $sdkArchive = [IO.Compression.ZipFile]::OpenRead($temporaryArchive)
    }
    $sdkName = 'opencascade-8.0.1-vc14-64'
    if (!$sdkArchive.GetEntry($sdkName + '/cmake/OpenCASCADEConfig.cmake')) {
        throw 'Expected the OCCT 8.0.1 vc14 x64 combined binary SDK.'
    }
    $prefixes = @(
        "$sdkName/inc/", "$sdkName/cmake/", "$sdkName/src/",
        "$sdkName/win64/vc14/bin/", "$sdkName/win64/vc14/lib/",
        "$sdkName/win64/vc14/libd/",
        '3rdparty-vc14-64/tbb-2021.13.0-x64/',
        '3rdparty-vc14-64/jemalloc-vc14-64/',
        '3rdparty-vc14-64/freetype-2.13.3-x64/',
        '3rdparty-vc14-64/freeimage-3.18.0-x64/',
        '3rdparty-vc14-64/openvr-1.14.15-64/',
        '3rdparty-vc14-64/ffmpeg-3.3.4-64/',
        '3rdparty-vc14-64/msvc-vc14-64/'
    )
    $files = 0
    $bytes = 0L
    foreach ($entry in $sdkArchive.Entries) {
        $name = $entry.FullName
        $selected = $name.StartsWith("$sdkName/win64/vc14/bind/") -and $name.EndsWith('.dll')
        foreach ($prefix in $prefixes) { if ($name.StartsWith($prefix)) { $selected = $true; break } }
        if (!$selected -or $name.EndsWith('/')) { continue }
        $targetPath = [IO.Path]::GetFullPath((Join-Path $targetRoot ($name.Replace('/', '\'))))
        if (!$targetPath.StartsWith($targetRoot + '\', [StringComparison]::OrdinalIgnoreCase)) {
            throw 'Archive path traversal rejected.'
        }
        New-Item -ItemType Directory -Path ([IO.Path]::GetDirectoryName($targetPath)) -Force | Out-Null
        [IO.Compression.ZipFileExtensions]::ExtractToFile($entry, $targetPath, $true)
        ++$files
        $bytes += $entry.Length
    }
    $runtimeRelative = @(
        "$sdkName\win64\vc14\bin",
        '3rdparty-vc14-64\tbb-2021.13.0-x64\bin',
        '3rdparty-vc14-64\jemalloc-vc14-64\bin',
        '3rdparty-vc14-64\freetype-2.13.3-x64\bin',
        '3rdparty-vc14-64\freeimage-3.18.0-x64\bin',
        '3rdparty-vc14-64\openvr-1.14.15-64\bin\win64',
        '3rdparty-vc14-64\ffmpeg-3.3.4-64\bin',
        '3rdparty-vc14-64\msvc-vc14-64'
    )
    $runtimePaths = @($runtimeRelative | ForEach-Object { Join-Path $targetRoot $_ })
    if ($ActivateRuntime) { $env:PATH = ($runtimePaths -join ';') + ';' + $env:PATH }
    [pscustomobject]@{
        OpenCASCADE_DIR = Join-Path $targetRoot "$sdkName\cmake"
        TBB_DIR = Join-Path $targetRoot '3rdparty-vc14-64\tbb-2021.13.0-x64\lib\cmake\TBB'
        RuntimePaths = $runtimePaths
        Files = $files
        ExpandedBytes = $bytes
    }
} finally {
    if ($sdkArchive -ne $outer) { $sdkArchive.Dispose() }
    $outer.Dispose()
    if ($temporaryArchive) { Remove-Item -LiteralPath $temporaryArchive }
}
