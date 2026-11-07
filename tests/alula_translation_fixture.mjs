// SPDX-License-Identifier: GPL-2.0-or-later
// Fault injection into a real, previously solved upstream result. Test only.
import {readFileSync,writeFileSync} from 'node:fs';
import {resolve} from 'node:path';
import {pathToFileURL} from 'node:url';
const [bridge,inputPath,outputPath,mutation,resultPath]=process.argv.slice(2);
const {translate}=await import(pathToFileURL(resolve(bridge)));
const input=JSON.parse(readFileSync(inputPath,'utf8'));
const saved=JSON.parse(readFileSync(outputPath,'utf8'));
const raw={...saved.upstream,...saved.upstream.coefficients,cm:saved.upstream.coefficients.cm_clockwise};
if(mutation==='provisional')raw.coefficientStatus='unconverged';
else if(mutation==='missing')delete raw.cl;
else if(mutation==='missing_diagnostic')delete raw.diagnostics.linearResidual;
else if(mutation==='nonfinite')raw.cl=Infinity;
else if(mutation==='outside_model'){
  input.mode='viscous_incompressible';input.settings.boundary_condition='streamfunction';
  input.settings.viscous={tolerance:1e-8,wake_tolerance:1e-6};
  raw.status='outside-model';raw.diagnostics={equationResidual:1e-10,wakeResidual:1e-8,maxRelativeDisplacement:.2};
  raw.cl=raw.cd=raw.cm=null;
}else throw new Error('Unknown test mutation');
writeFileSync(resultPath,JSON.stringify(translate(input,raw,{})),{flag:'wx'});
