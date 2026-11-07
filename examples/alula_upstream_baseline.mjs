// SPDX-License-Identifier: GPL-2.0-or-later
// Independent API baseline for integration verification, never a production CLI.
import {readFileSync,writeFileSync} from 'node:fs';
import {resolve,join} from 'node:path';
import {pathToFileURL} from 'node:url';
const [rootPath,inputPath,outputPath]=process.argv.slice(2);
const root=resolve(rootPath),input=JSON.parse(readFileSync(inputPath,'utf8'));
const {conditions:c,moment_reference:r,settings:s}=input;
const args={elements:input.elements.map(e=>({name:e.name,points:e.points})),alpha:c.alpha_deg,mach:c.mach,
  referenceChord:1,momentReference:{x:r.x_over_c,y:r.y_over_c},boundaryCondition:s.boundary_condition};
let raw;
if(input.mode==='inviscid'){
  const {solveInviscid}=await import(pathToFileURL(join(root,'src/inviscid/linear-vortex.js')));
  raw=solveInviscid(args);
}else if(input.mode==='viscous_incompressible'){
  const {solveViscousAssembly}=await import(pathToFileURL(join(root,'src/viscous/result.js')));
  const v=s.viscous;
  raw=solveViscousAssembly({...args,reynolds:c.reynolds,ncrit:v.ncrit,trips:v.trips,
    maxIterations:v.max_iterations,wakeLength:v.wake_length,wakeCount:v.wake_count,
    initialization:args.elements.length>1?'auto':'native',jacobianMode:'analytic'},
    {tolerance:v.tolerance,wakeTolerance:v.wake_tolerance,maxWakeIterations:v.max_wake_iterations,wakeRelaxation:v.wake_relaxation});
}else throw new Error('Unsupported baseline mode');
const out=Object.fromEntries(['status','cl','cd','cm','coefficientStatus','diagnostics','elements','solverSettings','warnings','limitations'].filter(k=>raw[k]!==undefined).map(k=>[k,raw[k]]));
writeFileSync(outputPath,JSON.stringify(out,(_,value)=>typeof value==='number'&&!Number.isFinite(value)?String(value):value,2)+'\n',{flag:'wx'});
