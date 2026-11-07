# 从项目中任意位置调用本脚本即可串行复现，无需同时启动多个程序。
$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
Push-Location -LiteralPath $projectRoot
try {
    & matlab -batch "addpath('src'); r=Mai_penalty('full'); assert(r.verification.rk4GravityMaxError<1e-10); assert(r.verification.dragAnalyticMaxError<1e-8); assert(r.verification.magnusPowerPass && r.verification.dragPowerPass);"
    if ($LASTEXITCODE -ne 0) { throw 'MATLAB计算或数值检验失败。' }
    & python (Join-Path $PSScriptRoot 'Mai独立复核.py')
    if ($LASTEXITCODE -ne 0) { throw 'Python独立数值复核失败。' }
    & python (Join-Path $PSScriptRoot 'Mai论文结果同步.py')
    if ($LASTEXITCODE -ne 0) { throw '英文稿数据同步失败。' }
    Push-Location -LiteralPath (Join-Path $projectRoot 'latex\AI')
    try {
        foreach ($paper in @('B题英文论文.tex','B题中文论文.tex')) {
            for ($pass=1; $pass -le 2; $pass++) {
                & xelatex -interaction=nonstopmode -halt-on-error $paper
                if ($LASTEXITCODE -ne 0) { throw "论文编译失败：$paper" }
            }
        }
    } finally { Pop-Location }
    Write-Host '固定题给场景已复现。若修改模型或输入，须重新核对中文正文和结论；本脚本不翻译或自动提交。'
} finally { Pop-Location }
