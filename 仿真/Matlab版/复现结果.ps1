# 从项目中任意位置调用本脚本即可串行复现，无需同时启动多个程序。
$ErrorActionPreference = 'Stop'
# 本脚本位于 项目根/仿真/Matlab版/；项目根需向上三级（Matlab版→仿真→项目根）
$simRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$projectRoot = $simRoot
Push-Location -LiteralPath $projectRoot
try {
    & matlab -batch "addpath('$PSScriptRoot'); penalty('full'); scenario_checks; scenario_ui_checks; penalty_figures; penalty_view('export'); penalty_view('demo');"
    if ($LASTEXITCODE -ne 0) { throw 'MATLAB计算或界面检查失败。' }
    & python (Join-Path $PSScriptRoot '论文结果同步.py')
    if ($LASTEXITCODE -ne 0) { throw '中英文主结果同步失败。' }
    Push-Location -LiteralPath (Join-Path $projectRoot 'latex\AI')
    try {
        foreach ($paper in @('B题英文论文.tex','B题中文论文.tex')) {
            for ($pass=1; $pass -le 2; $pass++) {
                & xelatex -interaction=nonstopmode -halt-on-error $paper
                if ($LASTEXITCODE -ne 0) { throw "论文编译失败：$paper" }
            }
        }
    } finally { Pop-Location }
    Write-Host '全程模型、符合观测的轨迹、MATLAB检查、双语配图及论文已生成。'
} finally { Pop-Location }
