/* Optional Node.js entry point for reproducible experiments. Browser users open index.html. */
'use strict';
const fs=require('node:fs');
const {createPhysics}=require('./physics.js'),{createDensity}=require('./density.js');
const {createStudy}=require('./study.js'),{createMapping}=require('./mapping.js');
const {createExperiments}=require('./experiments.js');
const P=createPhysics(),D=createDensity(),S=createStudy(P,D),M=createMapping(P,D,S),E=createExperiments(P,D,S,M);
const [kind,input,output]=process.argv.slice(2);
if(!['study','virtual','batch','default-batch'].includes(kind)||(!input&&kind!=='default-batch')){
 console.log('Usage: node run-lab.cjs study|virtual|batch input.json [output.json]\n       node run-lab.cjs default-batch [output.json]');process.exitCode=1;
}else{
 try{
  let result,destination=output;
  if(kind==='default-batch'){result=E.defaultConfig();destination=input;}
  else{const cfg=JSON.parse(fs.readFileSync(input,'utf8'));
   if(kind==='batch')result=E.runBatch(cfg);
   else if(kind==='virtual')result=E.runExperiment(cfg.reference,cfg.study,{observationSigma:cfg.study.sigma,seed:cfg.study.seed,observationsOverride:cfg.manual?cfg.study.base.observations:undefined,regionLevel:.95});
   else result=M.withMapping(S.runStudy(cfg.study||cfg));
  }
  const text=JSON.stringify(result,null,2)+'\n';
  if(destination)fs.writeFileSync(destination,text);else process.stdout.write(text);
 }catch(e){process.stderr.write(e.message+'\n');process.exitCode=1;}
}
