param(
    [string]$BaselineRef = '9a357cc7304f85439bf1d85926e3b9e9f866cb91',
    [switch]$Benchmark,
    [int]$Size = 100000,
    [int]$Threads = 4
)
$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot
Push-Location $root
try {
    $out = Join-Path $root 'output\solver_openmp'
    $baseline = Join-Path $root 'output\solver_openmp_baseline'
    New-Item -ItemType Directory -Force $out, $baseline | Out-Null
    foreach ($header in @('GMRES.h', 'CGS.h')) {
        # Use the byte stream; PowerShell text redirection can corrupt legacy
        # encoded comments in the original GMRES header.
        $process = New-Object System.Diagnostics.Process
        $process.StartInfo.FileName = 'git'
        $process.StartInfo.Arguments = "show ${BaselineRef}:Code/$header"
        $process.StartInfo.UseShellExecute = $false
        $process.StartInfo.RedirectStandardOutput = $true
        $process.StartInfo.RedirectStandardError = $true
        $null = $process.Start()
        $stream = [System.IO.File]::Create((Join-Path $baseline $header))
        try { $process.StandardOutput.BaseStream.CopyTo($stream) }
        finally { $stream.Dispose() }
        $errorText = $process.StandardError.ReadToEnd()
        $process.WaitForExit()
        if ($process.ExitCode -ne 0) { throw $errorText }
    }
    $vswhere = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"
    $vs = & $vswhere -latest -products '*' -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
    if (-not $vs) { throw 'MSVC C++ tools were not found.' }
    $vcvars = Join-Path $vs 'VC\Auxiliary\Build\vcvars64.bat'
    $flags = '/nologo /std:c++17 /EHsc /O2 /W4 /utf-8 /wd4819 /wd4828 /DMLFMA_TEST_BASELINE'
    foreach ($variant in @('openmp', 'disabled', 'serial')) {
        $extra = switch ($variant) {
            'openmp' { '/openmp' }
            'disabled' { '/openmp /DMLFMA_SOLVER_DISABLE_OPENMP' }
            'serial' { '' }
        }
        $exe = "output\solver_openmp\validation_$variant.exe"
        $obj = "output\solver_openmp\validation_$variant.obj"
        & cmd /d /c "`"$vcvars`" >nul && cl $flags $extra tests\solver_openmp_validation.cpp /Fe:$exe /Fo:$obj"
        if ($LASTEXITCODE -ne 0) { throw "Compilation failed: $variant" }
        & ".\$exe" > "output\solver_openmp\validation_$variant.txt" 2>&1
        if ($LASTEXITCODE -ne 0) { throw "Validation failed: $variant (see output/solver_openmp)." }
        Write-Output "$variant : PASS"
    }
    if ($Benchmark) {
        foreach ($variant in @('openmp', 'disabled')) {
            $extra = if ($variant -eq 'disabled') { '/DMLFMA_SOLVER_DISABLE_OPENMP' } else { '' }
            $exe = "output\solver_openmp\benchmark_$variant.exe"
            $obj = "output\solver_openmp\benchmark_$variant.obj"
            & cmd /d /c "`"$vcvars`" >nul && cl $flags /openmp $extra tests\solver_openmp_benchmark.cpp /Fe:$exe /Fo:$obj"
            if ($LASTEXITCODE -ne 0) { throw "Benchmark compilation failed: $variant" }
            & ".\$exe" $Size $Threads 5 50 > "output\solver_openmp\benchmark_${variant}_${Size}_${Threads}.csv" 2>&1
            if ($LASTEXITCODE -ne 0) { throw "Benchmark failed: $variant" }
        }
    }
}
finally { Pop-Location }
