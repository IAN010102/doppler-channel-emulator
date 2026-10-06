// Sentence 2 of item 10 (withSignal SMI, delta_theta = 0, L = 100, K = 20 dB): mean +- SE of the three runs that give -5.5 ... -5.7 dB (E3 seeds 3000+, E4 seeds 4000+, T17d seeds 1717+6000+), and the effect of the modulation (QPSK / 16-QAM) with the same seeds. node docs/diagnostics/withsignal_se_check.js
const Core=require('../../core.js'),X=require('../../experiments.js');
function run(over,seed0,n){const a=[];for(let t=0;t<n;t++){Core.setSeed(seed0+t);const s=Core.createSys();s.calZ=new Array(16).fill(0);Object.assign(s,X.BASE,over,{freshRealization:true});s.snaps=[];s.computeMath();a.push(s.sinrDb);}
 const m=a.reduce((x,y)=>x+y)/n,sd=Math.sqrt(a.reduce((x,y)=>x+(y-m)*(y-m),0)/(n-1));return m.toFixed(2)+' ± '+(sd/Math.sqrt(n)).toFixed(2)+' (sd '+sd.toFixed(2)+')';}
const o={algo:'SMI',trainMode:'withSignal',pointErrDeg:0};
console.log('E3 setting (QAM16, seeds 3000+):',run(o,3000,1000));
console.log('E4 setting (QAM16, seeds 4000+):',run(o,4000,1000));
console.log('T17d setting (QPSK, seeds 1717+6000+):',run(Object.assign({},o,{mod:'QPSK'}),7717,1000));
console.log('QAM16 with the T17d seeds:',run(o,7717,1000));
console.log('QPSK with the E3 seeds:',run(Object.assign({},o,{mod:'QPSK'}),3000,1000));
