'use strict';
// Optional browser QA. Runtime app has no npm dependency.
// npm install playwright (in an isolated test directory), or set PLAYWRIGHT_MODULE.
const fs=require('node:fs'),path=require('node:path'),os=require('node:os');
const assert=require('node:assert/strict');
const {pathToFileURL}=require('node:url');
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'playwright');
const source=path.resolve(__dirname),packageRoot=path.resolve(source,'../..');
const artifacts=fs.mkdtempSync(path.join(os.tmpdir(),'UPC 网页验收 '));
const channels=(process.env.TEST_BROWSER_CHANNELS||'msedge,chrome').split(',');
const report={schemaVersion:'browser-acceptance-1',date:new Date().toISOString(),node:process.version,platform:process.platform,offline:true,checks:[],browsers:[],passed:false};
const record=(name,detail={})=>report.checks.push({name,passed:true,...detail});
async function download(page,id,name){const promise=page.waitForEvent('download');await page.locator('#'+id).click();const d=await promise;const target=path.join(artifacts,name);await d.saveAs(target);assert(fs.statSync(target).size>0);return target;}
async function recalc(p){const prev=await p.evaluate(()=>window.UPCApp.run?.runId);await p.locator('#runBtn').click();await p.waitForFunction(id=>window.UPCApp.run?.runId!==id,prev,{timeout:120000});return p.evaluate(()=>window.UPCApp.run);}
async function imported(p,c){const f=path.join(artifacts,'配置重载.json');fs.writeFileSync(f,JSON.stringify(c));await p.locator('#configFile').setInputFiles(f);await p.waitForFunction(()=>document.querySelector('#status').textContent.includes('已导入'));return recalc(p);}
(async()=>{
 for(const channel of channels){
  const browser=await chromium.launch({channel,headless:true});
  try{
   const relocated=path.join(artifacts,channel+' 中文和空格');fs.cpSync(source,relocated,{recursive:true});
   const context=await browser.newContext({viewport:{width:1440,height:1100},offline:true,acceptDownloads:true});
   const p=await context.newPage(),errors=[],network=[];p.on('pageerror',e=>errors.push(e.message));p.on('request',r=>{if(/^https?:/.test(r.url()))network.push(r.url());});
   const start=Date.now();await p.goto(pathToFileURL(path.join(relocated,'index.html')).href);
   await p.waitForFunction(()=>window.UPCApp?.run?.prediction?.tGoal>0,null,{timeout:120000});
   const baseline=await p.evaluate(()=>window.UPCApp.run),initialComputeMs=Date.now()-start;
   assert(Math.abs(baseline.prediction.yGoal-2.5994405492)<1e-4);
   record(channel+': A01/A10 offline relocated startup',{initialComputeMs});
   if(channel===channels[0])await p.screenshot({path:path.join(packageRoot,'img/网页首版.png'),fullPage:true});
   for(let j=0;j<6;j++){
    await p.locator('#officialBtn').click();await p.locator('#obs'+Math.floor(j/3)+(j%3)).fill(String(baseline.config.observations[Math.floor(j/3)].position[j%3]+.005));
    assert(await p.locator('#staleOverlay').isVisible());const r=await recalc(p);
    assert(Math.hypot(r.prediction.yGoal-baseline.prediction.yGoal,r.prediction.zGoal-baseline.prediction.zGoal)>1e-4);
   }
   record(channel+': A02 six observation coordinates refit');
   await p.locator('#officialBtn').click();await recalc(p);
   for(const model of ['M0','M1','M2','M3']){await p.locator('#model').selectOption(model);const r=await recalc(p);assert.equal(r.config.model,model);assert.equal(r.prediction.status,'crossed_goal_plane');}
   record(channel+': A03 M0–M3 actual solver');
   await p.locator('[data-tab="uncertainty"]').click();await p.locator('#sampleCount').fill('200');const t=Date.now();await p.locator('#uncertaintyBtn').click();
   await p.waitForFunction(()=>window.UPCApp.statistics?.samples?.length===200,null,{timeout:120000});
   const statistics=await p.evaluate(()=>window.UPCApp.statistics);assert.equal(statistics.summary.total,200);
   record(channel+': A06 bootstrap200',{elapsedMs:Date.now()-t,summary:statistics.summary});
   await p.locator('[data-tab="export"]').click();const sampleFile=await download(p,'samplesCsvBtn',channel+'-samples.csv');assert.equal(fs.readFileSync(sampleFile,'utf8').trim().split(/\r?\n/).length,201);
   await p.locator('[data-tab="uncertainty"]').click();await p.locator('#uncertaintyMethod').selectOption('scenario');await p.locator('#uncertaintyBtn').click();
   await p.waitForFunction(()=>window.UPCApp.statistics?.method==='scenario');assert.equal(await p.evaluate(()=>window.UPCApp.statistics.samples.length),12);assert.equal(await p.evaluate(()=>window.UPCApp.statistics.summary.conditionalEllipse),null);
   record(channel+': A06 deterministic 12-scenario mode has no probability ellipse');
   await p.locator('[data-tab="comparison"]').click();
   for(const mode of ['refit','fixed']){await p.locator('#compareMode').selectOption(mode);await p.locator('#compareBtn').click();await p.waitForFunction(m=>window.UPCApp.comparisons?.comparison===m,mode);assert.equal(await p.evaluate(()=>window.UPCApp.comparisons.runs.length),4);}
   await download(p,'exportPngBtn',channel+'-comparison.png');await p.locator('[data-tab="export"]').click();const note=await download(p,'figureNoteBtn',channel+'-caption.md');assert(fs.readFileSync(note,'utf8').includes('"model": "M3"'));
   record(channel+': model comparisons fixed/refit and caption metadata');
   const configFile=await download(p,'saveConfigBtn',channel+'-config.json');const saved=JSON.parse(fs.readFileSync(configFile,'utf8'));
   await p.close();const restored=await context.newPage();restored.on('pageerror',e=>errors.push(e.message));
   await restored.goto(pathToFileURL(path.join(relocated,'index.html')).href);await restored.waitForFunction(()=>window.UPCApp?.run?.prediction?.tGoal>0);
   const roundtrip=await imported(restored,saved);assert(Math.abs(roundtrip.prediction.yGoal-baseline.prediction.yGoal)<1e-6);
   record(channel+': A08 close page, import and rerun');
   const a=restored;
   const forward={...saved,mode:'forward',model:'M3',velocity:[30.4,9.4,6.9],b:[.2,-.1,-.3],goalHalfWidth:2.5,goalHeight:2.1};
   const changed=await imported(a,forward);assert.equal(changed.config.b[0],.2);assert.equal(changed.config.goalHalfWidth,2.5);assert.equal(changed.config.goalHeight,2.1);
   await download(a,'exportPngBtn',channel+'-custom-goal.png');record(channel+': imported forward bx and custom goal retained');
   await a.locator('#officialBtn').click();await a.locator('#mode').selectOption('forward');await a.locator('#model').selectOption('M4');
   const m4=await recalc(a);assert.equal(m4.prediction.status,'crossed_goal_plane');assert(m4.trajectory[50].q.length===4);
   await a.locator('#timeline').fill('750');await a.locator('#timeline').dispatchEvent('input');assert(!(await a.locator('#timeReadout').innerText()).startsWith('0.000'));
   await a.locator('#forceToggle').check();await a.locator('[data-tab="curves"]').click();await a.locator('#curveType').selectOption('omega');
   if(channel===channels[0])await a.screenshot({path:path.join(packageRoot,'img/网页真实自旋.png'),fullPage:true});
   await a.locator('[data-tab="export"]').click();const csv=await download(a,'trajectoryCsvBtn',channel+'-M4.csv');const lines=fs.readFileSync(csv,'utf8').replace(/^\uFEFF/,'').trim().split(/\r?\n/),header=lines[0].split(',').map(s=>s.replaceAll('"',''));
   assert(header.includes('omega_z_rad_s')&&header.includes('q_w')&&header.includes('torque_z_N_m')&&header.includes('relativeAirVelocity_x_m_s'));
   assert.equal(lines.length,m4.trajectory.length+1);assert(lines.every(l=>l.split(',').length===header.length));
   const resultFile=await download(a,'resultJsonBtn',channel+'-M4.json');const exported=JSON.parse(fs.readFileSync(resultFile,'utf8'));assert.equal(exported.run.config.model,'M4');assert(exported.run.provenance.verification.includes('未完成实验'));
   const sheet=await download(a,'sequenceExportBtn',channel+'-12frames.png');
   const pngBytes=fs.readFileSync(sheet);assert(pngBytes.includes(Buffer.from('UPC-Metadata')));
   assert.equal(pngBytes.readUInt32BE(16),2400);assert.equal(pngBytes.readUInt32BE(20),1200);
   const sheetMetaFile=await download(a,'sequenceMetadataBtn',channel+'-12frames.json');
   const sheetMeta=JSON.parse(fs.readFileSync(sheetMetaFile,'utf8'));assert.equal(sheetMeta.frames,12);assert.equal(sheetMeta.config.model,'M4');assert.equal(sheetMeta.times_s.length,12);assert(sheetMeta.coreVersion);
   if(channel===channels[0])fs.copyFileSync(sheet,path.join(packageRoot,'img/交互端12帧演示.png'));
   await a.locator('#sequenceExportBtn').click();await a.locator('#cancelSequenceBtn').click();await a.waitForFunction(()=>document.querySelector('#sequenceExportBtn').disabled===false);assert((await a.locator('#sequenceStatus').innerText()).includes('取消'));
   record(channel+': A04/A05/A09 M4 time, curves, complete export and cancellable frame contact sheet');
   await a.locator('#dragModel').selectOption('paper2018');await a.locator('#liftModel').selectOption('teamgeist2018');const empirical=await recalc(a);assert.equal(empirical.prediction.status,'crossed_goal_plane');
   await a.locator('#v0').fill('10');const invalid=await recalc(a);assert.equal(invalid.prediction.status,'outside_model_domain');assert.equal(await a.locator('#status').getAttribute('class'),'error');
   record(channel+': experimental relation and model-domain failure visible');
   await a.locator('#officialBtn').click();await recalc(a);await a.locator('details').filter({has:a.locator('#airTemp')}).evaluate(e=>e.open=true);await a.locator('#airTemp').fill('20');await a.locator('#airRH').fill('50');await a.locator('#applyAirBtn').click();assert(Math.abs(Number(await a.locator('#rho').inputValue())-1.198869387)<1e-8);
   record(channel+': moist-air assumption helper');
   await a.locator('#officialBtn').click();await recalc(a);await a.locator('[data-tab="uncertainty"]').click();await a.locator('#sigma').fill('0');await recalc(a);await a.locator('#uncertaintyBtn').click();await a.waitForFunction(()=>document.querySelector('#status').className==='error');assert.equal(await a.evaluate(()=>window.UPCApp.statistics),null);
   record(channel+': A07 zero assumed noise disables probability sampling');
   await a.locator('#officialBtn').click();await recalc(a);await a.locator('#sampleCount').fill('600');await a.locator('#uncertaintyBtn').click();await a.locator('#cancelBtn').click();await a.waitForTimeout(100);assert.equal(await a.evaluate(()=>window.UPCApp.statistics),null);
   await a.locator('#uncertaintyBtn').click();await a.locator('#obs00').fill('1.51');assert((await a.locator('#status').innerText()).includes('取消'));const newest=await recalc(a);await a.waitForTimeout(250);assert.equal((await a.evaluate(()=>window.UPCApp.run)).runId,newest.runId);
   record(channel+': A11 cancellation and changed-input stale-job isolation');
   await a.locator('#officialBtn').click();await a.locator('#obs00').fill('20');await a.locator('#obs10').fill('40');const incompatible=await recalc(a);assert.equal(incompatible.fit.diagnostics.compatibility.flag,'incompatible_with_declared_noise');assert((await a.locator('#resultBadge').innerText()).includes('残差过大'));
   record(channel+': A12 incompatible data is visibly conditional');
   await a.setViewportSize({width:390,height:844});assert(!(await a.evaluate(()=>document.documentElement.scrollWidth>innerWidth+2)));record(channel+': mobile layout at 390px');
   assert.deepEqual(errors,[]);assert.deepEqual(network,[]);report.browsers.push({channel,version:await browser.version(),noPageErrors:true,noNetworkRequests:true});
  } finally {await browser.close();}
 }
 report.passed=true;fs.writeFileSync(path.join(packageRoot,'data/验证基准/网页验收.json'),JSON.stringify(report,null,2));console.log(JSON.stringify({passed:true,checks:report.checks.length,browsers:report.browsers,artifacts},null,2));
})().catch(e=>{console.error(e);process.exitCode=1;});
