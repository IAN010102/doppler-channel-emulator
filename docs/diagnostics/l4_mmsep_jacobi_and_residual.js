// does the Jacobi stopping tolerance of core.pinvHermitian explain the ~1e-10 difference? (scratch copy of the routine with a tolerance argument; repo untouched)
const path=require('path').join(__dirname,'..','..')+'/'; const Core=require(path+'core.js'), J=require(path+'data/matlab_cases.json'); const {Cplx}=Core;
function pinvTol(H,eps,tol,maxSweep){const n=H.length,m=2*n,A=Array.from({length:m},()=>new Float64Array(m)),V=Array.from({length:m},(_,i)=>{const r=new Float64Array(m);r[i]=1;return r;});
 for(let i=0;i<n;i++)for(let j=0;j<n;j++){const a=H[i][j].r,b=H[i][j].i;A[i][j]=a;A[i+n][j+n]=a;A[i][j+n]=-b;A[i+n][j]=b;}
 let scale=0;for(let i=0;i<m;i++)scale+=A[i][i]*A[i][i];scale=scale||1;let sw=0;
 for(;sw<maxSweep;sw++){let off=0;for(let p=0;p<m;p++)for(let q=p+1;q<m;q++)off+=A[p][q]*A[p][q];if(off<tol*scale)break;
  for(let p=0;p<m-1;p++)for(let q=p+1;q<m;q++){if(Math.abs(A[p][q])<1e-300)continue;const th=(A[q][q]-A[p][p])/(2*A[p][q]),t=(th>=0?1:-1)/(Math.abs(th)+Math.sqrt(th*th+1)),c=1/Math.sqrt(t*t+1),s=t*c;
   for(let k=0;k<m;k++){const x=A[k][p],y=A[k][q];A[k][p]=c*x-s*y;A[k][q]=s*x+c*y;}for(let k=0;k<m;k++){const x=A[p][k],y=A[q][k];A[p][k]=c*x-s*y;A[q][k]=s*x+c*y;}for(let k=0;k<m;k++){const x=V[k][p],y=V[k][q];V[k][p]=c*x-s*y;V[k][q]=s*x+c*y;}}}
 let lmax=0;for(let i=0;i<m;i++)lmax=Math.max(lmax,A[i][i]);const keep=[];for(let i=0;i<m;i++)if(A[i][i]>eps*lmax)keep.push(i);
 return {P:Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>{let re=0,im=0;for(const k of keep){re+=V[i][k]*V[j][k]/A[k][k];im+=V[i+n][k]*V[j][k]/A[k][k];}return new Cplx(re,im);})),sweeps:sw};}
const cj=a=>({r:a.r,i:-a.i}),mul=(a,b)=>({r:a.r*b.r-a.i*b.i,i:a.r*b.i+a.i*b.r}),add=(a,b)=>({r:a.r+b.r,i:a.i+b.i}),abs2=a=>a.r*a.r+a.i*a.i;
const c=J.cases.find(x=>x.id==='static_L4_signalFree_d3'),P=c.params,N=P.N,L=P.L;
const X=c.X.re.map((row,i)=>row.map((r,j)=>new Cplx(r,c.X.im[i][j]))), s=c.s.re.map((r,i)=>new Cplx(r,c.s.im[i])), a=c.a_assumed.re.map((r,i)=>new Cplx(r,c.a_assumed.im[i])), h=c.h.re.map((r,i)=>new Cplx(r,c.h.im[i])), aJ=c.a_J.re.map((r,i)=>new Cplx(r,c.a_J.im[i]));
const Rw=Array.from({length:N},(_,m)=>Array.from({length:N},(_,q)=>{let t=new Cplx(0,0);for(let l=0;l<L;l++)t=Cplx.add(t,Cplx.mul(X[m][l],Cplx.conj(X[q][l])));return new Cplx(t.r/L,t.i/L);}));
const r=Array.from({length:N},(_,m)=>{let t=new Cplx(0,0);for(let l=0;l<L;l++)t=Cplx.add(t,Cplx.mul(X[m][l],Cplx.conj(s[l])));return new Cplx(t.r/L,t.i/L);});
const wWeb=c.algorithms['MMSE-P']?c.algorithms['MMSE-P']:c.algorithms.MMSE_P, ww=wWeb.w.re.map((x,i)=>new Cplx(x,wWeb.w.im[i]));
const nrm=w=>{let d=new Cplx(0,0);for(let i=0;i<N;i++)d=Cplx.add(d,Cplx.mul(Cplx.conj(a[i]),w[i]));const q=new Cplx(d.r/(d.r*d.r+d.i*d.i),-d.i/(d.r*d.r+d.i*d.i));return w.map(x=>Cplx.mul(x,q));};
const dn=(u,v)=>Math.sqrt(u.reduce((t,x,i)=>t+abs2({r:x.r-v[i].r,i:x.i-v[i].i}),0))/Math.sqrt(v.reduce((t,x)=>t+abs2(x),0));
const vd=(u,v)=>{let t=new Cplx(0,0);for(let i=0;i<N;i++)t=Cplx.add(t,Cplx.mul(Cplx.conj(u[i]),v[i]));return t;};
const sinr=w=>{const wh=vd(w,h),wa=vd(w,aJ);return 10*Math.log10((wh.r*wh.r+wh.i*wh.i)/(P.P_j*(wa.r*wa.r+wa.i*wa.i)+P.sigma2*w.reduce((t,x)=>t+abs2(x),0)));};
const ref=Core.matMulVec(pinvTol(Rw,1e-10,1e-300,400).P,r);   // converged to machine precision
console.log('json w vs the converged Jacobi solution: rel err',dn(nrm(ww),nrm(ref)).toExponential(2),' SINR diff',(sinr(ww)-sinr(ref)).toExponential(2));
for(const tol of [1e-30,1e-34,1e-40,1e-60,1e-300]){const o=pinvTol(Rw,1e-10,tol,400),w=Core.matMulVec(o.P,r);console.log('tol',tol,'sweeps',o.sweeps,' rel err vs json',dn(nrm(w),nrm(ww)).toExponential(2),' vs converged',dn(nrm(w),nrm(ref)).toExponential(2),' SINR diff to converged',(sinr(w)-sinr(ref)).toExponential(2));}

// ---- residual of the consistent system X^H w = conj(s) (exact solution has residual 0): which candidate is closer to the exact minimum-norm solution?
function qrRef(X,s,N,L){const cols=[];for(let l=0;l<L;l++)cols.push(Array.from({length:N},(_,i)=>X[i][l]));
 const Q=[],R=Array.from({length:L},()=>Array.from({length:L},()=>new Cplx(0,0)));
 for(let l=0;l<L;l++){let v=cols[l].map(c=>new Cplx(c.r,c.i));for(let pass=0;pass<2;pass++)for(let k=0;k<l;k++){let p=new Cplx(0,0);for(let i=0;i<N;i++)p=Cplx.add(p,Cplx.mul(Cplx.conj(Q[k][i]),v[i]));R[k][l]=Cplx.add(R[k][l],p);v=v.map((c,i)=>Cplx.sub(c,Cplx.mul(Q[k][i],p)));}
  const nv=Math.sqrt(v.reduce((t,x)=>t+abs2(x),0));R[l][l]=new Cplx(nv,0);Q.push(v.map(c=>new Cplx(c.r/nv,c.i/nv)));}
 const y=[];for(let l=0;l<L;l++){let t=Cplx.conj(s[l]);for(let k=0;k<l;k++)t=Cplx.sub(t,Cplx.mul(Cplx.conj(R[k][l]),y[k]));const d=Cplx.conj(R[l][l]);y.push(new Cplx((t.r*d.r+t.i*d.i)/abs2(d),(t.i*d.r-t.r*d.i)/abs2(d)));}
 const w=Array.from({length:N},()=>new Cplx(0,0));for(let l=0;l<L;l++)for(let i=0;i<N;i++)w[i]=Cplx.add(w[i],Cplx.mul(Q[l][i],y[l]));return w;}
function resid(w){let t=0,tn=0;for(let l=0;l<L;l++){let u=new Cplx(0,0);for(let i=0;i<N;i++)u=Cplx.add(u,Cplx.mul(Cplx.conj(X[i][l]),w[i]));const e=Cplx.sub(u,Cplx.conj(s[l]));t+=abs2(e);tn+=abs2(s[l]);}return Math.sqrt(t/tn);}
const wQ=qrRef(X,s,N,L);
console.log('residual |X^H w - conj(s)| / |s|:  web json',resid(ww).toExponential(2),'| converged Jacobi pinv',resid(ref).toExponential(2),'| QR (Gram-Schmidt x2)',resid(wQ).toExponential(2));
console.log('web vs QR rel err (w^H a = 1):',dn(nrm(ww),nrm(wQ)).toExponential(2),' ; converged-Jacobi vs QR:',dn(nrm(ref),nrm(wQ)).toExponential(2),' ; ||w|| web',Math.sqrt(ww.reduce((t,x)=>t+abs2(x),0)).toFixed(4),' QR',Math.sqrt(wQ.reduce((t,x)=>t+abs2(x),0)).toFixed(4));
// is R_hat (formed from X) the same as the matrix the web used? the web stores R_hat_train_web (signalFree: from X_sf); compare Rw formed here with a second summation order
const Rw2=Array.from({length:N},(_,m)=>Array.from({length:N},(_,q)=>{let t=new Cplx(0,0);for(let l=L-1;l>=0;l--)t=Cplx.add(t,Cplx.mul(X[m][l],Cplx.conj(X[q][l])));return new Cplx(t.r/L,t.i/L);}));
let dm=0,nm=0;for(let m=0;m<N;m++)for(let q=0;q<N;q++){dm=Math.max(dm,Math.hypot(Rw[m][q].r-Rw2[m][q].r,Rw[m][q].i-Rw2[m][q].i));nm=Math.max(nm,Math.hypot(Rw[m][q].r,Rw[m][q].i));}
console.log('R_hat summation order (forward vs backward) changes entries by',(dm/nm).toExponential(2),'relative');
const w3=Core.matMulVec(pinvTol(Rw2,1e-10,1e-300,400).P,r);console.log('pinv solution from the backward-summed R_hat vs forward-summed: rel err',dn(nrm(w3),nrm(ref)).toExponential(2),' SINR diff',(sinr(w3)-sinr(ref)).toExponential(2));
