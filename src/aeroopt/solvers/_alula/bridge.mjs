// SPDX-License-Identifier: GPL-2.0-or-later
// AeroOpt integration, 2026. No Alula numerical algorithms are copied here.
import {readFileSync, writeFileSync} from 'node:fs';
import {createHash} from 'node:crypto';
import {resolve, dirname, join} from 'node:path';
import {fileURLToPath, pathToFileURL} from 'node:url';
import {performance} from 'node:perf_hooks';

export const VERSION='0.4.0';
const here=dirname(fileURLToPath(import.meta.url));
const lockBytes=readFileSync(join(here,'upstream-lock.json'));
const lock=JSON.parse(lockBytes);
const sha=data=>createHash('sha256').update(data).digest('hex');
const fail=message=>{throw new Error(message);};
const finite=(value,label)=>{if(typeof value!=='number'||!Number.isFinite(value))fail(`${label} must be finite`);return value;};
const integer=(value,min,max,label)=>{finite(value,label);if(!Number.isInteger(value)||value<min||value>max)fail(`Invalid ${label}`);};
const range=(value,min,max,label)=>{finite(value,label);if(value<min||value>max)fail(`Invalid ${label}`);};
function keys(value,names,label){
  if(!value||typeof value!=='object'||Array.isArray(value)||Object.keys(value).sort().join('|')!==[...names].sort().join('|'))fail(`Invalid ${label} fields`);
}
function point(value){keys(value,['x','y'],'point');finite(value.x,'x');finite(value.y,'y');}
const same=(a,b)=>a.x===b.x&&a.y===b.y;
const equalPoints=(a,b)=>Array.isArray(a)&&a.length===b.length&&a.every((p,i)=>same(p,b[i]));

export function validateInput(input){
  keys(input,['schema_version','case_id','mode','elements','conditions','moment_reference','settings','runtime'],'input');
  if(input.schema_version!=='aeroopt.alula-input.v1'||!/^[a-f0-9]{64}$/.test(input.case_id))fail('Invalid schema or case ID');
  if(!['inviscid','viscous_incompressible'].includes(input.mode))fail('Unsupported mode');
  keys(input.conditions,['alpha_deg','reynolds','mach'],'conditions');
  const c=input.conditions,viscous=input.mode==='viscous_incompressible';
  range(c.alpha_deg,viscous?-20:-90,viscous?20:90,'alpha');
  finite(c.mach,'Mach');if(c.mach!==0)fail('Mach must equal zero');
  finite(c.reynolds,'Reynolds');if(c.reynolds<=0)fail('Reynolds must be positive');
  if(viscous)range(c.reynolds,5e4,5e7,'Reynolds');
  keys(input.moment_reference,['x_over_c','y_over_c','positive_sense','normalization'],'moment reference');
  const ref=input.moment_reference;finite(ref.x_over_c,'moment x');finite(ref.y_over_c,'moment y');
  if(!['clockwise','counterclockwise'].includes(ref.positive_sense)||ref.normalization!=='q_inf*c_main^2')fail('Unsupported moment convention');
  keys(input.runtime,['node_version','upstream_commit'],'runtime');
  if(input.runtime.node_version!==process.version||input.runtime.upstream_commit!==lock.commit)fail('Runtime identity mismatch');
  const settings=input.settings;
  keys(settings,['boundary_condition','viscous'],'settings');
  if(!['normal-velocity','streamfunction'].includes(settings.boundary_condition))fail('Unknown boundary condition');
  if(!Array.isArray(input.elements)||input.elements.length<1||input.elements.length>6)fail('Supply 1-6 elements');
  let total=0;
  for(const [index,e] of input.elements.entries()){
    keys(e,['name','role','points','closure'],'element');
    if(typeof e.name!=='string'||!e.name||e.role!==(index===0?'main':index===1?'flap':`element_${index+1}`))fail('Element order/role mismatch');
    if(!Array.isArray(e.points)||e.points.length<9)fail('At least nine supplied points required');
    e.points.forEach(point);
    const closed=same(e.points[0],e.points.at(-1)),n=e.points.length-1;
    if(n+(closed?0:1)>600)fail('At most 600 solid panels per element');
    if(e.closure!==(closed?'sharp_repeated_endpoint':'straight_TE_cap'))fail('Closure declaration differs from supplied endpoints');
    if(!closed&&settings.boundary_condition!=='streamfunction')fail('Finite bases require explicitly selected streamfunction formulation');
    const p0=e.points[0];let area=0;
    for(let i=0;i<e.points.length;i++){
      const a=e.points[i],b=e.points[(i+1)%e.points.length];
      area+=(a.x-p0.x)*(b.y-p0.y)-(b.x-p0.x)*(a.y-p0.y);
      if(i<e.points.length-1&&same(a,b))fail('Repeated successive points');
    }
    if(!(area>0))fail('Counterclockwise contours required; no automatic reordering');
    if(viscous&&n<40)fail('Viscous mode requires at least 40 wetted panels per element');
    total+=n;
  }
  if(total>700)fail('At most 700 total wetted panels');
  if(viscous){
    if(settings.boundary_condition!=='streamfunction')fail('Viscous mode uses the streamfunction outer solution');
    const s=settings.viscous;
    keys(s,['ncrit','trips','max_iterations','wake_length','wake_count','tolerance','wake_tolerance','max_wake_iterations','wake_relaxation'],'viscous settings');
    range(s.ncrit,0,14,'Ncrit');
    if(!Array.isArray(s.trips)||s.trips.length!==2)fail('Two explicit transition trips required');
    s.trips.forEach(t=>{range(t,Number.MIN_VALUE,1,'trip');});
    integer(s.max_iterations,1,200,'Newton iteration count');integer(s.wake_count,4,200,'wake count');
    integer(s.max_wake_iterations,0,200,'wake iteration count');
    for(const k of ['wake_length','tolerance','wake_tolerance']){finite(s[k],k);if(s[k]<=0)fail(`${k} must be positive`);}
    range(s.wake_relaxation,Number.MIN_VALUE,1,'wake relaxation');
  }else if(settings.viscous!==null)fail('Inviscid input cannot request boundary-layer parameters');
  return input;
}

export async function solve(input,upstream){
  const bridgeStart=performance.now();
  validateInput(input);
  for(const [name,expected] of Object.entries(lock.files)){
    if(sha(readFileSync(join(upstream,name)))!==expected)fail(`Pinned Alula source mismatch: ${name}`);
  }
  const importStart=performance.now();
  const integrityMs=importStart-bridgeStart;
  const {prepareAirfoilElement}=await import(pathToFileURL(join(upstream,'src/geometry/airfoil-element.js')));
  const module=await import(pathToFileURL(join(upstream,lock.entry_points[input.mode])));
  const importMs=performance.now()-importStart;
  const elements=input.elements.map(e=>({name:e.name,points:structuredClone(e.points)}));
  const prepared=elements.map(prepareAirfoilElement);
  prepared.forEach((e,i)=>{
    const p=elements[i].points,expected=same(p[0],p.at(-1))?p:[...p,p[0]];
    if(!equalPoints(e.points,expected))fail('Upstream preparation changed supplied coordinates/order');
  });
  const ref=input.moment_reference,c=input.conditions,s=input.settings.viscous;
  const args={elements,alpha:c.alpha_deg,mach:c.mach,referenceChord:1,
    momentReference:{x:ref.x_over_c,y:ref.y_over_c},boundaryCondition:input.settings.boundary_condition};
  const start=performance.now();
  const raw=input.mode==='inviscid'?module.solveInviscid(args):module.solveViscousAssembly({
    ...args,reynolds:c.reynolds,ncrit:s.ncrit,trips:s.trips,maxIterations:s.max_iterations,
    wakeLength:s.wake_length,wakeCount:s.wake_count,initialization:elements.length>1?'auto':'native',jacobianMode:'analytic'
  },{tolerance:s.tolerance,wakeTolerance:s.wake_tolerance,maxWakeIterations:s.max_wake_iterations,wakeRelaxation:s.wake_relaxation});
  const solverMs=performance.now()-start;
  return translate(input,raw,{validation_integrity_ms:integrityMs,import_ms:importMs,solver_ms:solverMs,
    bridge_ms:performance.now()-bridgeStart,node_uptime_ms:process.uptime()*1000});
}

export function translate(input,raw,timings){
  const inv=input.mode==='inviscid',d=raw.diagnostics??{},s=input.settings.viscous;
  const warnings=[...(raw.warnings??[])];
  const limits=['RESEARCH_MODEL_NOT_PHYSICALLY_VALIDATED'];
  if(input.elements.some(e=>e.closure==='straight_TE_cap'))limits.push('ALULA_FINITE_BASE_RESEARCH');
  const gate=(key,limit)=>typeof d[key]==='number'&&Number.isFinite(d[key])&&d[key]>=0&&d[key]<=limit;
  const extraInviscidKeys=input.settings.boundary_condition==='normal-velocity'?['normalVelocityResidual']:
    ['nodalStreamfunctionResidual','interiorTEResidual'];
  let numerical=(raw.status==='solved'||raw.status==='outside-model')&&(inv?gate('linearResidual',1e-9)&&gate('kuttaResidual',1e-9)&&extraInviscidKeys.every(k=>gate(k,1e-9)):
    gate('equationResidual',s.tolerance)&&gate('wakeResidual',s.wake_tolerance));
  const modelAdmissible=inv||gate('maxRelativeDisplacement',.1);
  const provisional=raw.coefficientStatus!==undefined&&raw.coefficientStatus!=='available'&&raw.coefficientStatus!=='converged';
  const required=inv?['cl','cm']:['cl','cd','cm'];
  for(const key of ['cl','cd','cm'])if(raw[key]!==null&&raw[key]!==undefined)finite(raw[key],`upstream ${key}`);
  if(inv&&raw.cd!==null&&raw.cd!==undefined)fail('Inviscid drag cannot be a physical coefficient');
  const present=required.every(k=>typeof raw[k]==='number'&&Number.isFinite(raw[k]));
  if(!Array.isArray(raw.elements)||raw.elements.length!==input.elements.length)fail('Upstream element count mismatch');
  for(const [i,e] of raw.elements.entries()){
    const p=input.elements[i].points,closed=same(p[0],p.at(-1))?p:[...p,p[0]];
    if(!equalPoints(e.points,closed))fail('Returned upstream solid geometry differs from input');
    if(!Array.isArray(e.cp)||numerical&&(!e.cp.length||e.cp.some(q=>!Number.isFinite(q.cp)||!Number.isFinite(q.x)||!Number.isFinite(q.y))))fail('Missing/invalid pressure data');
  }
  const nonfinite=[];
  function clean(value,path=''){
    if(typeof value==='number'&&!Number.isFinite(value)){nonfinite.push({path,value:String(value)});return null;}
    if(Array.isArray(value)||ArrayBuffer.isView(value))return Array.from(value,(v,i)=>clean(v,`${path}/${i}`));
    if(value&&typeof value==='object')return Object.fromEntries(Object.entries(value).filter(([,v])=>v!==undefined).map(([k,v])=>[k,clean(v,`${path}/${k}`)]));
    return value;
  }
  const upstream=clean({model:raw.model,status:raw.status,coefficientStatus:raw.coefficientStatus,
    coefficients:{cl:raw.cl??null,cd:raw.cd??null,cm_clockwise:raw.cm??null},diagnostics:d,
    elements:raw.elements,boundaryLayer:raw.boundaryLayer,solverSettings:raw.solverSettings,
    history:raw.history,wakeHistory:raw.wakeHistory,limitations:raw.limitations});
  if(nonfinite.length){numerical=false;warnings.push('Nonfinite upstream diagnostics retained as null with explicit path/value labels.');}
  const eligible=numerical&&modelAdmissible&&raw.status==='solved'&&present&&!provisional;
  if(!eligible)warnings.push('Coefficients withheld: convergence, availability and eligibility are separate requirements.');
  const factor=input.moment_reference.positive_sense==='counterclockwise'?-1:1;
  return {schema_version:'aeroopt.alula-output.v1',case_id:input.case_id,mode:input.mode,
    conditions:input.conditions,moment_reference:input.moment_reference,
    status:{upstream:raw.status,numerical_convergence:numerical?'converged':'not_converged',
      model_admissible:modelAdmissible,coefficient_status:provisional?'provisional':present?'available':'missing',eligible},
    coefficients:{cl:eligible?raw.cl:null,cd:eligible&&!inv?raw.cd:null,cm:eligible?factor*raw.cm:null},
    geometry:input.elements,warnings,limitations:limits,nonfinite_diagnostics:nonfinite,upstream,
    metadata:{bridge_version:VERSION,bridge_sha256:sha(readFileSync(fileURLToPath(import.meta.url))),
      upstream_commit:lock.commit,lock_sha256:sha(lockBytes),node_version:process.version,
      numerical_dependencies:lock.runtime_dependencies,reference_chord:1,upstream_moment_positive:'clockwise',
      timings,memory:process.memoryUsage(),process_versions:process.versions,numerical_eligibility_only:true,
      comparison_scope:'controlled_same_model_comparison_requires_refinement',aerodynamic_ranking_qualified:false,
      physical_validation:false}};
}

async function main(){
  const args=process.argv.slice(2);
  if(args.length!==6||args[0]!=='--upstream'||args[2]!=='--input'||args[4]!=='--output')fail('Expected --upstream PATH --input JSON --output JSON');
  const data=readFileSync(resolve(args[3]));
  if(data.length>10_000_000)fail('Input exceeds 10 MB');
  const input=JSON.parse(data.toString('utf8'));
  const output=await solve(input,resolve(args[1]));
  output.input_file_sha256=sha(data);
  writeFileSync(resolve(args[5]),JSON.stringify(output,null,2)+'\n',{flag:'wx'});
}
if(process.argv[1]&&resolve(process.argv[1])===fileURLToPath(import.meta.url)){
  main().catch(error=>{console.error(`ALULA_BRIDGE_FAILED: ${error.stack??error}`);process.exitCode=2;});
}
