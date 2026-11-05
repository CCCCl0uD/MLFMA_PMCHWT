param(
    [string]$DataRoot = ''
)
$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot
if (-not $DataRoot) {
    $DataRoot = Join-Path (Split-Path -Parent $root) 'DATA'
}
$vswhere = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"
$vs = & $vswhere -latest -products '*' -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
if (-not $vs) { throw 'MSVC C++ tools were not found.' }
$vcvars = Join-Path $vs 'VC\Auxiliary\Build\vcvars64.bat'
Push-Location $root
try {
    New-Item -ItemType Directory -Force 'output\nas_rwg_validation' | Out-Null
    & cmd /d /c "`"$vcvars`" >nul && cl /nologo /std:c++17 /EHsc /W4 /utf-8 /Od /Zi /fsanitize=address tests\nas_rwg_validation.cpp Code\NASReader.cpp Code\RWGGenerator.cpp Code\OCTree.cpp /Fe:output\nas_rwg_validation\nas_rwg_validation.exe /Fo:output\nas_rwg_validation\ /Fd:output\nas_rwg_validation\validation.pdb /link /INCREMENTAL:NO"
    if ($LASTEXITCODE -ne 0) { throw 'NAS/RWG validation compilation failed.' }
    & cmd /d /c "`"$vcvars`" >nul && output\nas_rwg_validation\nas_rwg_validation.exe output\nas_rwg_validation\fixtures `"$DataRoot`" > output\nas_rwg_validation\validation.log 2>&1"
    if ($LASTEXITCODE -ne 0) { throw 'Validation failed; see output/nas_rwg_validation/validation.log.' }
    Get-Content -LiteralPath 'output\nas_rwg_validation\validation.log' -Tail 1
}
finally {
    Pop-Location
}
