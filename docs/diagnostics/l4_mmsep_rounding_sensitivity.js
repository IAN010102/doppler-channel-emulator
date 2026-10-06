// L = 4 realisations (unified, signalFree, K = 20 dB): sensitivity of the pinv MMSE-P weights to a change of the summation order of R_hat (~1e-16 relative) versus kappa of the non-zero part
const Core=require(require('path').join(__dirname,'..','..','core.js')); const {Cplx}=Core;
function pinvC(H,eps){return Core.pinvHermitian(H,eps).pinv;}
const rows=[];
for(let t=0;t<400;t++){Core.setSeed(31000+t);const s=Core.createSys();s.calZ=new Array(16).fill(0);Object.assign(s,{model:'unified',algo:'MMSEP',trainMode:'signalFree',N:8,L:4,aoaT:0,aoaJ:40,snr:20,sir:-10,v:0,latMs:0,kDb:20,calDeg:0,mod:'QPSK',freshRealization:true});s.snaps=[];s.computeMath();
 const N=8,L=4,x=s.snaps.map(sn=>Array.from({length:N},(_,i)=>new Cplx(sn.rr[i]+sn.tr[i],sn.ri[i]+sn.ti[i]))),sy=s.snaps.map(sn=>new Cplx(sn.s1r,sn.s1i));
 const form=order=>Array.from({length:N},(_,m)=>Array.from({length:N},(_,q)=>{let a=new Cplx(0,0);for(const l of order)a=Cplx.add(a,Cplx.mul(x[l][m],Cplx.conj(x[l][q])));return new Cplx(a.r/L,a.i/L);}));
 const r=Array.from({length:N},(_,m)=>{let a=new Cplx(0,0);for(let l=0;l<L;l++)a=Cplx.add(a,Cplx.mul(x[l][m],Cplx.conj(sy[l])));return new Cplx(a.r/L,a.i/L);});
 const Rf=form([0,1,2,3]),Rb=form([3,2,1,0]);
 const M=Array.from({length:16},()=>new Float64Array(16));for(let i=0;i<N;i++)for(let j=0;j<N;j++){const a=Rf[i][j].r,b=Rf[i][j].i;M[i][j]=a;M[i+N][j+N]=a;M[i][j+N]=-b;M[i+N][j]=b;}
 const ev=Core.eigSymDecomp(M).vals.filter((_,k)=>k%2===0).slice(0,4),kap=ev[0]/ev[3];
 const w1=Core.matMulVec(pinvC(Rf,1e-10),r),w2=Core.matMulVec(pinvC(Rb,1e-10),r);
 const dn=Math.sqrt(w1.reduce((a,c,i)=>a+(c.r-w2[i].r)**2+(c.i-w2[i].i)**2,0))/Math.sqrt(w1.reduce((a,c)=>a+c.r*c.r+c.i*c.i,0));
 rows.push({kap,dn});}
const lk=rows.map(r=>Math.log10(r.kap)),ld=rows.map(r=>Math.log10(Math.max(r.dn,1e-18)));
const m=a=>a.reduce((x,y)=>x+y)/a.length,mk=m(lk),md=m(ld);let sxy=0,sxx=0,syy=0;for(let i=0;i<rows.length;i++){sxy+=(lk[i]-mk)*(ld[i]-md);sxx+=(lk[i]-mk)**2;syy+=(ld[i]-md)**2;}
const sorted=rows.map(r=>r.dn).sort((a,b)=>a-b);
console.log('400 realisations: kappa_nz median',rows.map(r=>r.kap).sort((a,b)=>a-b)[200].toExponential(2),' (min',Math.min(...rows.map(r=>r.kap)).toExponential(1),'max',Math.max(...rows.map(r=>r.kap)).toExponential(1),')');
console.log('relative change of w from reordering the sum: median',sorted[200].toExponential(2),' p90',sorted[360].toExponential(2),' max',sorted[399].toExponential(2));
console.log('log-log: slope d log(err)/d log(kappa) =',(sxy/sxx).toFixed(2),' correlation',(sxy/Math.sqrt(sxx*syy)).toFixed(2),'; share of realisations with change > 1e-9:',(rows.filter(r=>r.dn>1e-9).length/rows.length).toFixed(3));
