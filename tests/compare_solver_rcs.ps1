param(
    [Parameter(Mandatory=$true)][string]$ParallelDirectory,
    [Parameter(Mandatory=$true)][string]$BaselineDirectory,
    [double]$RelativePowerTolerance = 1e-4
)
$ErrorActionPreference = 'Stop'
if ($RelativePowerTolerance -le 0) { throw 'RelativePowerTolerance must be positive.' }
$files = @(Get-ChildItem -LiteralPath $ParallelDirectory -Filter '*_RCS_*.txt' -Recurse)
if ($files.Count -eq 0) { throw 'No RCS output files found.' }
foreach ($file in $files) {
    $relativeFolder = $file.DirectoryName.Substring((Get-Item -LiteralPath $ParallelDirectory).FullName.Length).TrimStart('\','/')
    $referenceFolder = if ($relativeFolder) { Join-Path $BaselineDirectory $relativeFolder } else { $BaselineDirectory }
    $pattern = $file.Name -replace '_RCS_.*\.txt$', '_RCS_*.txt'
    $reference = @(Get-ChildItem -LiteralPath $referenceFolder -Filter $pattern)
    if ($reference.Count -ne 1) { throw "Missing or ambiguous reference: $pattern" }
    $currentRows = @(Get-Content -LiteralPath $file.FullName | Where-Object { $_ -match '^\s*[-+]?\d' })
    $referenceRows = @(Get-Content -LiteralPath $reference[0].FullName | Where-Object { $_ -match '^\s*[-+]?\d' })
    if ($currentRows.Count -eq 0 -or $currentRows.Count -ne $referenceRows.Count) { throw "RCS row count mismatch: $pattern" }
    $maxPowerError = 0.0
    $maxTotalDbError = 0.0
    $maxAnyDbError = 0.0
    for ($row=0; $row -lt $currentRows.Count; ++$row) {
        $a = $currentRows[$row].Trim() -split '\s+'
        $b = $referenceRows[$row].Trim() -split '\s+'
        if ($a.Count -ne 6 -or $b.Count -ne 6) { throw "Unexpected RCS format: $pattern" }
        for ($column=0; $column -lt 3; ++$column) {
            if ([double]$a[$column] -ne [double]$b[$column]) { throw "Angle/index mismatch: $pattern" }
        }
        # Normalize all polarization-channel power errors to total scattered
        # power. A relative dB test on an almost-zero cross channel is unstable.
        $totalPower = [Math]::Pow(10, [double]$b[3]/10)
        if (-not ($totalPower -gt 0) -or [double]::IsInfinity($totalPower)) { throw "Invalid reference total power: $pattern" }
        for ($column=3; $column -lt 6; ++$column) {
            $currentDb = [double]$a[$column]
            $referenceDb = [double]$b[$column]
            if ([double]::IsNaN($currentDb) -or [double]::IsNaN($referenceDb) -or
                [double]::IsInfinity($currentDb) -or [double]::IsInfinity($referenceDb)) { throw "Nonfinite RCS: $pattern" }
            $dbError = [Math]::Abs($currentDb - $referenceDb)
            $maxAnyDbError = [Math]::Max($maxAnyDbError, $dbError)
            if ($column -eq 3) { $maxTotalDbError = [Math]::Max($maxTotalDbError, $dbError) }
            $powerError = [Math]::Abs([Math]::Pow(10,$currentDb/10) - [Math]::Pow(10,$referenceDb/10)) / $totalPower
            $maxPowerError = [Math]::Max($maxPowerError, $powerError)
        }
    }
    [pscustomobject]@{
        Case = "$relativeFolder/$pattern"
        MaxRelativePowerError = $maxPowerError
        MaxTotalDbError = $maxTotalDbError
        MaxAnyDbError = $maxAnyDbError
    }
    if ($maxPowerError -gt $RelativePowerTolerance) { throw "RCS power mismatch: $pattern" }
}
