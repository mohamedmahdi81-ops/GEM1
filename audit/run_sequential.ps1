# Requires: conda envs "GEM1" (main pipeline) and "gem1-troppo" (Step 5/extraction)
# already created -- see repository README.md for setup. Resolves conda's own
# env directory via `conda info --base` so this script is portable across
# machines/checkouts rather than hardcoding an absolute install path.
$ErrorActionPreference = "Continue"
Set-Location $PSScriptRoot
$condaBase = (conda info --base).Trim()
$pyMain = Join-Path $condaBase "envs\GEM1\python.exe"
$pyTroppo = Join-Path $condaBase "envs\gem1-troppo\python.exe"

# Run sequentially, not concurrently -- this machine has 16GB RAM and each of
# these jobs (FastCC's dense matrices, cobra+Gurobi FVA, OptGP sampling) is
# memory-heavy enough that running them together caused an OS paging-file
# exhaustion / MemoryError (see audit/logs/item6_err.log from the first,
# concurrent attempt). Mirrors the project's own established
# scripts/_run_step7_then_step8.ps1 chaining pattern.

Write-Output "=== Item 6 (troppo re-extraction) started $(Get-Date -Format o) ==="
& $pyTroppo "item6_troppo_reextract.py" 1> "logs\item6.log" 2> "logs\item6_err.log"
$item6Exit = $LASTEXITCODE
Write-Output "=== Item 6 finished $(Get-Date -Format o), exit code $item6Exit ==="

Write-Output "=== Item 1 (iMAT+FVA baseline) started $(Get-Date -Format o) ==="
& $pyMain "item1_imat_fva_baseline.py" 1> "logs\item1.log" 2> "logs\item1_err.log"
$item1Exit = $LASTEXITCODE
Write-Output "=== Item 1 finished $(Get-Date -Format o), exit code $item1Exit ==="

Write-Output "=== Item 3 (flux sampling + perturbation) started $(Get-Date -Format o) ==="
& $pyMain "item3_flux_sampling_perturbation.py" 1> "logs\item3.log" 2> "logs\item3_err.log"
$item3Exit = $LASTEXITCODE
Write-Output "=== Item 3 finished $(Get-Date -Format o), exit code $item3Exit ==="

Write-Output "=== All sequential audit items done $(Get-Date -Format o) ==="
Write-Output "Exit codes: item6=$item6Exit item1=$item1Exit item3=$item3Exit"
