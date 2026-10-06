// read-only investigation of the L < N cases (the full-matrix condition number is meaningless there: the dropped eigenvalues are rounding noise) of data/matlab_cases.json (MMSE-P): nothing in the repo is modified
const path=require('path').join(__dirname,'..','..')+'/';
const Core=require(path+'core.js'), J=require(path+'data/matlab_cases.json');
const cz=z=>z.re.map((r,i)=>({r,i:z.im[i]}));   // check the layout below
function mat(z){ // z.re, z.im are N x L arrays (row major) possibly nested
  return z.re.map((row,i)=>row.map((r,j)=>({r,i:z.im[i][j]})));
}
function vec(z){ return z.re.map((r,i)=>({r,i:z.im[i]})); }
const mul=(a,b)=>({r:a.r*b.r-a.i*b.i,i:a.r*b.i+a.i*b.r}), cj=a=>({r:a.r,i:-a.i}), add=(a,b)=>({r:a.r+b.r,i:a.i+b.i}), sub=(a,b)=>({r:a.r-b.r,i:a.i-b.i});
const abs2=a=>a.r*a.r+a.i*a.i;
function dot(u,v){ let s={r:0,i:0}; for(let k=0;k<u.length;k++) s=add(s,mul(cj(u[k]),v[k])); return s; } // u^H v
function norm(u){ return Math.sqrt(u.reduce((s,c)=>s+abs2(c),0)); }
for(const c of J.cases){
  const P=c.params; if(!(P.L<P.N)) continue;
  const N=P.N,L=P.L, Xm=mat(c.X), s=c.s.re.map((r,i)=>({r,i:c.s.im[i]})), a=vec(c.a_assumed), h=vec(c.h), aJ=vec(c.a_J);
  const cols=[]; for(let l=0;l<L;l++) cols.push(Xm.map(row=>row[l]));        // columns of X (N-vectors)
  // Rw = X X^H / L
  const Rw=Array.from({length:N},(_,m)=>Array.from({length:N},(_,q)=>{let t={r:0,i:0};for(let l=0;l<L;l++)t=add(t,mul(cols[l][m],cj(cols[l][q])));return {r:t.r/L,i:t.i/L};}));
  const M=Array.from({length:2*N},()=>new Float64Array(2*N));
  for(let i=0;i<N;i++)for(let j=0;j<N;j++){const x=Rw[i][j].r,y=Rw[i][j].i;M[i][j]=x;M[i+N][j+N]=x;M[i][j+N]=-y;M[i+N][j]=y;}
  const ev=Core.eigSymDecomp(M).vals.filter((_,k)=>k%2===0), lmax=ev[0];
  const kept=ev.filter(v=>v>P.epsRank*lmax), dropped=ev.filter(v=>v<=P.epsRank*lmax);
  // reference: minimum-norm solution of X^H w = conj(s)  via modified Gram-Schmidt QR of X (applied twice): w = Q R^-H conj(s)
  const Q=[],R=Array.from({length:L},()=>Array.from({length:L},()=>({r:0,i:0})));
  for(let l=0;l<L;l++){let v=cols[l].map(c=>({...c}));for(let pass=0;pass<2;pass++)for(let k=0;k<l;k++){const p=dot(Q[k],v);R[k][l]=add(R[k][l],p);v=v.map((c,i)=>sub(c,mul(Q[k][i],p)));}const nv=norm(v);R[l][l]={r:nv,i:0};Q.push(v.map(c=>({r:c.r/nv,i:c.i/nv})));}
  // solve R^H y = conj(s) (forward substitution); y_l
  const sc=s.slice(0,L).map(cj), y=[];
  for(let l=0;l<L;l++){let t=sc[l];for(let k=0;k<l;k++)t=sub(t,mul(cj(R[k][l]),y[k]));const d=cj(R[l][l]);y.push({r:(t.r*d.r+t.i*d.i)/abs2(d),i:(t.i*d.r-t.r*d.i)/abs2(d)});}
  let wRef=Array.from({length:N},()=>({r:0,i:0})); for(let l=0;l<L;l++)for(let i=0;i<N;i++)wRef[i]=add(wRef[i],mul(Q[l][i],y[l]));
  // singular values of X: sqrt of the kept eigenvalues of X X^H (L of them) -> kappa(X)
  const kapX=Math.sqrt(kept[0]/kept[kept.length-1]);
  const nrm=w=>{const d=dot(a,w);const q={r:d.r/abs2(d),i:-d.i/abs2(d)};return w.map(x=>mul(x,q));};  // w / (a^H w)
  const wWeb=vec(c.algorithms['MMSE-P']?c.algorithms['MMSE-P'].w:c.algorithms.MMSE_P.w), nW=nrm(wWeb), nR=nrm(wRef);
  const eW=norm(nW.map((x,i)=>sub(x,nR[i])))/norm(nR);
  const sinr=w=>{const wh=dot(w,h),wa=dot(w,aJ);return 10*Math.log10(abs2(wh)/(P.P_j*abs2(wa)+P.sigma2*w.reduce((t,x)=>t+abs2(x),0)));};
  console.log(c.id.padEnd(26),'| eigs kept:',kept.map(x=>x.toExponential(2)).join(' '),'| dropped:',dropped.map(x=>x.toExponential(2)).join(' '));
  console.log('   kappa(R_nz)=',(kept[0]/kept[kept.length-1]).toExponential(2),' kappa(X)=',kapX.toExponential(2),' smallest kept/lmax=',(kept[kept.length-1]/lmax).toExponential(2),' largest dropped/lmax=',(dropped[0]/lmax).toExponential(2),'  threshold eps=',P.epsRank,
    '\n   web w vs QR reference: rel err (w^H a=1)',eW.toExponential(2),'  SINR web',c.algorithms['MMSE-P']? '' : '', 'SINR(web w)-SINR(QR w)=',(sinr(wWeb)-sinr(wRef)).toExponential(2),' dB; json sinr',(c.algorithms.MMSE_P||c.algorithms['MMSE-P']).sinr_dB, ' (kappa_nz*2.2e-16 =',(kept[0]/kept[kept.length-1]*2.2e-16).toExponential(2),')');
}
