param(
    [ValidateSet('quick','full')][string]$Mode = 'full',
    [string]$PythonExe = ''
)
# 只验证模型，不同步论文、不编译、不提交。每次运行生成独立目录。
$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
$runsRoot = Join-Path $projectRoot 'data\验证运行'
$before = @()
if (Test-Path -LiteralPath $runsRoot) { $before = @(Get-ChildItem -LiteralPath $runsRoot -Directory | ForEach-Object Name) }
$matlabSource = $PSScriptRoot.Replace("'", "''")
& matlab -batch "addpath('$matlabSource'); report=Mai_p0_validate('$Mode');"
if ($LASTEXITCODE -ne 0) { throw 'MATLAB验证失败，不能宣称运行完成。' }
$newRuns = @(Get-ChildItem -LiteralPath $runsRoot -Directory | Where-Object { $_.Name -notin $before })
if ($newRuns.Count -ne 1) { throw '无法唯一确定新运行目录；请检查是否有并行任务。' }
$run = $newRuns[0].FullName
if ($PythonExe) {
    & $PythonExe (Join-Path $PSScriptRoot 'Mai独立复核.py') --matlab (Join-Path $run 'Mai罚球结果.json') --output (Join-Path $run '独立复核结果.json')
    if ($LASTEXITCODE -ne 0) { throw '独立数值复核失败。请检查Python路径和固定依赖清单。' }
} else {
    Write-Host 'MATLAB验证通过；本次未运行Python独立复核。可通过 -PythonExe 指定已安装固定依赖的解释器。'
}
Write-Host "本次证据目录：$run"
