
// ===================== Organic molecules =====================
const ORG_VAL={B:[3],C:[4],N:[3,5],O:[2],P:[3,5],S:[2,4,6],F:[1],Cl:[1],Br:[1],I:[1]};
function parseSMILES(src){
  const s=String(src).trim();if(!s||/\s/.test(s))throw chemErr('Not a SMILES string.');
  const atoms=[],bonds=[],ring={},stack=[];let prev=-1,bo=null,i=0;
  const mk=(a,b,o)=>[a,b,o===null?(atoms[a].arom&&atoms[b].arom?'ar':1):o];
  const add=a=>{atoms.push(a);const k=atoms.length-1;if(prev>=0)bonds.push(mk(prev,k,bo));prev=k;bo=null};
  const bad=m=>chemErr(m||'Not a valid SMILES string.');
  while(i<s.length){const c=s[i];
    if(c==='('){if(prev<0)throw bad();stack.push(prev);i++;continue}
    if(c===')'){if(!stack.length)throw bad();prev=stack.pop();i++;continue}
    if(c==='-'||c==='/'||c==='\\'){bo=1;i++;continue}
    if(c==='='){bo=2;i++;continue}if(c==='#'){bo=3;i++;continue}if(c===':'){bo='ar';i++;continue}
    if(c==='.'){prev=-1;bo=null;i++;continue}
    if(/\d/.test(c)||c==='%'){let n;if(c==='%'){n=s.substr(i+1,2);i+=3}else{n=c;i++}
      if(prev<0)throw bad();
      if(ring[n]){const r=ring[n];delete ring[n];bonds.push(mk(r.atom,prev,bo!==null?bo:r.bo));bo=null}else{ring[n]={atom:prev,bo};bo=null}
      continue}
    if(c==='['){const j=s.indexOf(']',i);if(j<0)throw bad();const t=s.slice(i+1,j);i=j+1;
      const m=t.match(/^(\d*)(se|as|[cnospb]|[A-Z][a-z]?)(@*)(?:(H)(\d*))?([+-]+\d*)?(?::\d+)?$/);
      if(!m)throw bad();let sym=m[2];const arom=/^[a-z]/.test(sym);if(arom)sym=sym[0].toUpperCase()+sym.slice(1);
      if(!EL[sym]){if(sym.length===2&&EL[sym[0]])throw bad();throw bad()}
      const h=m[4]?(m[5]?+m[5]:1):0;let fc=0;
      if(m[6]){const ch=m[6];if(/^[+-]\d+$/.test(ch))fc=(ch[0]==='-'?-1:1)*+ch.slice(1);else fc=[...ch].reduce((a,x)=>a+(x==='+'?1:-1),0)}
      add({el:sym,fc,arom,h});continue}
    const two=s.substr(i,2);
    if(two==='Cl'||two==='Br'){add({el:two,fc:0,arom:false,h:null});i+=2;continue}
    if('BCNOPSFI'.includes(c)){add({el:c,fc:0,arom:false,h:null});i++;continue}
    if('bcnosp'.includes(c)){add({el:c.toUpperCase(),fc:0,arom:true,h:null});i++;continue}
    throw bad();
  }
  if(Object.keys(ring).length||stack.length||!atoms.length)throw bad();
  const deg=atoms.map(()=>0);bonds.forEach(([a,b,o])=>{const v=o==='ar'?1:o;deg[a]+=v;deg[b]+=v});
  atoms.forEach((a,k)=>{if(a.h!==null)return;const sum=deg[k]+(a.arom?1:0);const v=(ORG_VAL[a.el]||[]).find(x=>x>=sum);a.h=v===undefined?0:v-sum});
  const TG={C:4,N:3,O:2,S:2,P:3,B:3};
  const need=atoms.map((a,k)=>{if(!a.arom)return 0;const t=(TG[a.el]||4)+(a.el==='C'?-Math.abs(a.fc):a.fc);return t-deg[k]-a.h===1?1:0});
  const matched=new Set();
  const solve=()=>{const u=atoms.findIndex((a,k)=>need[k]&&!matched.has(k));if(u<0)return true;
    for(const b of bonds){if(b[2]!=='ar')continue;const v=b[0]===u?b[1]:b[1]===u?b[0]:-1;if(v<0||!need[v]||matched.has(v))continue;
      b[2]=2;matched.add(u);matched.add(v);if(solve())return true;b[2]='ar';matched.delete(u);matched.delete(v)}
    return false};
  if(!solve())throw chemErr('I could not work out the bonds in the aromatic ring. Try writing it with alternating = bonds, e.g. C1=CC=CC=C1.');
  bonds.forEach(b=>{if(b[2]==='ar')b[2]=1});
  const g={atoms:atoms.map(a=>({el:a.el,fc:a.fc})),bonds:bonds.map(b=>b.slice()),heavy:atoms.length};
  atoms.forEach((a,k)=>{for(let t=0;t<a.h;t++){g.atoms.push({el:'H',fc:0});g.bonds.push([k,g.atoms.length-1,1])}});
  return g;
}
function graphComp(g){const c={};g.atoms.forEach(a=>c[a.el]=(c[a.el]||0)+1);return c}
// put heavy atoms first (H last); permute coordinate arrays the same way
function heavyFirst(g,...coordSets){const order=[...g.atoms.keys()].filter(i=>g.atoms[i].el!=='H').concat([...g.atoms.keys()].filter(i=>g.atoms[i].el==='H'));
  const inv=[];order.forEach((o,n)=>inv[o]=n);
  const ng={atoms:order.map(i=>g.atoms[i]),bonds:g.bonds.map(([a,b,o])=>[inv[a],inv[b],o]),heavy:g.atoms.filter(a=>a.el!=='H').length};
  return [ng,...coordSets.map(cs=>cs?order.map(i=>cs[i]):null)];
}
function orgON(g){
  const on=g.atoms.map(a=>a.fc),det=g.atoms.map(()=>({}));
  g.bonds.forEach(([a,b,o])=>{const A=g.atoms[a].el,B=g.atoms[b].el,ea=EL[A].en,eb=EL[B].en;const s=A===B||ea===eb?0:(eb>ea?1:-1);
    on[a]+=s*o;on[b]-=s*o;
    const put=(x,el,sg)=>{const k=el+'|'+o;const d=det[x][k]=det[x][k]||{el,o,c:0,s:sg};d.c++};put(a,B,s);put(b,A,-s)});
  return {on,det};
}
// ---------- 2D layout (heavy atoms) ----------
function findRings(N,adj){
  const rings=[],seen=new Set();
  for(let a=0;a<N;a++)for(const[b]of adj[a]){if(b<a)continue;
    // shortest path a->b avoiding edge a-b
    const prev=Array(N).fill(-2);prev[a]=-1;const q=[a];
    while(q.length){const u=q.shift();if(u===b)break;for(const[v]of adj[u]){if(u===a&&v===b)continue;if(prev[v]===-2){prev[v]=u;q.push(v)}}}
    if(prev[b]===-2)continue;const path=[];for(let u=b;u!==-1;u=prev[u])path.push(u);
    if(path.length>8)continue;const key=[...path].sort((x,y)=>x-y).join(',');if(seen.has(key))continue;seen.add(key);rings.push(path)}
  return rings.sort((x,y)=>y.length-x.length===0?0:x.length-y.length);
}
function rot2(v,a){const c=Math.cos(a),s=Math.sin(a);return[v[0]*c-v[1]*s,v[0]*s+v[1]*c]}
function layout2D(g){
  const N=g.heavy,adj=[...Array(N)].map(()=>[]);
  g.bonds.forEach(([a,b,o])=>{if(a<N&&b<N){adj[a].push([b,o]);adj[b].push([a,o])}});
  const pos=Array(N).fill(null),rings=findRings(N,adj),inRing=atom=>rings.filter(r=>r.includes(atom));
  const center=Array(N).fill(null),inDir=Array(N).fill(null),side=Array(N).fill(1);
  const placeRing=(R,anchor,apos,out)=>{const n=R.length,rad=1/(2*Math.sin(Math.PI/n));const c=[apos[0]+out[0]*rad,apos[1]+out[1]*rad];
    const k0=R.indexOf(anchor);const a0=Math.atan2(apos[1]-c[1],apos[0]-c[0]);
    for(let k=0;k<n;k++){const at=R[(k0+k)%n];if(pos[at])continue;const a=a0+k*2*Math.PI/n;pos[at]=[c[0]+rad*Math.cos(a),c[1]+rad*Math.sin(a)];center[at]=c}};
  const isSp=u=>adj[u].some(([,o])=>o===3)||adj[u].filter(([,o])=>o===2).length>=2;
  const queue=[];
  const comps=[];const done=new Set();
  for(let s0=0;s0<N;s0++){if(pos[s0])continue;
    const off=comps.length?comps.reduce((m,c)=>Math.max(m,c),-1e9)+2.5:0;
    const R=inRing(s0)[0];
    if(R){placeRing(R,s0,[off,0],[0,1])}else{pos[s0]=[off,0]}
    const startQ=R?R.slice():[s0];queue.push(...startQ);
    while(queue.length){const u=queue.shift();if(done.has(u))continue;done.add(u);
      const kids=adj[u].filter(([v])=>!pos[v]).map(([v,o])=>[v,o]);if(!kids.length)continue;
      let dirs=[];
      if(center[u]&&inRing(u).length){let b=[pos[u][0]-center[u][0],pos[u][1]-center[u][1]];const l=Math.hypot(...b);b=[b[0]/l,b[1]/l];
        dirs=kids.length===1?[b]:kids.length===2?[rot2(b,0.5236),rot2(b,-0.5236)]:[rot2(b,0.9),b,rot2(b,-0.9)]}
      else if(inDir[u]){const d=inDir[u];
        if(kids.length===1)dirs=[isSp(u)?d:rot2(d,side[u]*Math.PI/3)];
        else if(kids.length===2)dirs=[rot2(d,Math.PI/3),rot2(d,-Math.PI/3)];
        else dirs=[rot2(d,Math.PI/2),d,rot2(d,-Math.PI/2)];}
      else{const n=kids.length;const base=n===1?[Math.cos(Math.PI/6),Math.sin(Math.PI/6)]:[1,0];dirs=[...Array(n)].map((_,k)=>rot2(base,k*2*Math.PI/Math.max(n,isSp(u)?2:n)));
        if(n===2&&!isSp(u))dirs=[rot2([1,0],Math.PI/6),rot2([1,0],Math.PI/6+2*Math.PI/3)]}
      kids.forEach(([v],k)=>{const d=dirs[k%dirs.length];const p=[pos[u][0]+d[0],pos[u][1]+d[1]];
        const Rv=inRing(v).find(r=>r.every(x=>!pos[x]));
        if(Rv){placeRing(Rv,v,p,d);if(!pos[v])pos[v]=p;queue.push(...Rv)}else{pos[v]=p;center[v]=null}
        inDir[v]=d;side[v]=-side[u];queue.push(v)});
    }
    comps.push(Math.max(...pos.filter(Boolean).map(p=>p[0])));
  }
  // relax: bond springs + overlap repulsion
  for(let it=0;it<250;it++){const F=pos.map(()=>[0,0]);
    g.bonds.forEach(([a,b])=>{if(a>=N||b>=N)return;const dx=pos[b][0]-pos[a][0],dy=pos[b][1]-pos[a][1],d=Math.hypot(dx,dy)||1e-3,f=(d-1)*0.3;F[a][0]+=f*dx/d;F[a][1]+=f*dy/d;F[b][0]-=f*dx/d;F[b][1]-=f*dy/d});
    for(let a=0;a<N;a++)for(let b=a+1;b<N;b++){const dx=pos[b][0]-pos[a][0],dy=pos[b][1]-pos[a][1],d=Math.hypot(dx,dy)||1e-3;if(d<1.25){const f=(1.25-d)*0.25;F[a][0]-=f*dx/d;F[a][1]-=f*dy/d;F[b][0]+=f*dx/d;F[b][1]+=f*dy/d}}
    pos.forEach((p,i)=>{p[0]+=F[i][0];p[1]+=F[i][1]})}
  return pos;
}
// ---------- 3D model from 2D layout + quick force-field cleanup ----------
function embed3D(g,p2){
  const n=g.atoms.length,N=g.heavy,adj=[...Array(n)].map(()=>[]);g.bonds.forEach(([a,b,o])=>{adj[a].push([b,o]);adj[b].push([a,o])});
  const P=Array(n);let seed=7;const rnd=()=>{seed=(seed*9301+49297)%233280;return seed/233280-0.5};
  for(let i=0;i<N;i++)P[i]=[p2[i][0]*1.45,p2[i][1]*1.45,rnd()*0.05];
  const hyb=i=>{const os=adj[i].map(([,o])=>o);if(os.includes(3)||os.filter(o=>o===2).length>=2)return 1;if(os.includes(2))return 2;
    if(['N','O'].includes(g.atoms[i].el)&&adj[i].some(([j])=>j<N&&adj[j].some(([,o])=>o===2)))return 2;return 3};
  for(let i=0;i<N;i++){const hs=adj[i].filter(([j])=>j>=N).map(([j])=>j);if(!hs.length)continue;
    const nb=adj[i].filter(([j])=>j<N).map(([j])=>V3.norm(V3.sub(P[j],P[i])));
    let a=nb.length?V3.mul(nb.reduce((s,v)=>V3.add(s,v),[0,0,0]),-1):[1,0,0];if(V3.len(a)<0.2)a=nb.length?V3.cross(nb[0],[0,0,1]):[1,0,0];a=V3.norm(a);
    let p=V3.cross(a,[0,0,1]);if(V3.len(p)<0.1)p=[1,0,0];p=V3.norm(p);const z=V3.norm(V3.cross(a,p));
    const h=hyb(i),m=hs.length;let dirs;
    if(m===1)dirs=[a];
    else if(m===2)dirs=h===2?[V3.add(V3.mul(a,.5),V3.mul(p,.87)),V3.add(V3.mul(a,.5),V3.mul(p,-.87))]:[V3.add(V3.mul(a,.58),V3.mul(z,.82)),V3.add(V3.mul(a,.58),V3.mul(z,-.82))];
    else if(m===3)dirs=[V3.add(V3.mul(a,.33),V3.mul(p,.94)),V3.add(V3.add(V3.mul(a,.33),V3.mul(p,-.47)),V3.mul(z,.82)),V3.add(V3.add(V3.mul(a,.33),V3.mul(p,-.47)),V3.mul(z,-.82))];
    else dirs=IDEAL[4];
    hs.forEach((j,k)=>{P[j]=V3.add(P[i],V3.mul(V3.norm(dirs[k%dirs.length]),bondLen(g.atoms[i].el,'H',1)))})}
  for(let i=0;i<n;i++)if(!P[i])P[i]=[rnd()*3,rnd()*3,rnd()*3];
  // restraints
  const R=[];const pairKey=new Set();const key=(a,b)=>a<b?a+'_'+b:b+'_'+a;
  g.bonds.forEach(([a,b,o])=>{R.push([a,b,bondLen(g.atoms[a].el,g.atoms[b].el,o),1]);pairKey.add(key(a,b))});
  for(let j=0;j<n;j++){const nb=adj[j];if(nb.length<2)continue;const h=hyb(j);const th=h===1?Math.PI:h===2?2*Math.PI/3:1.9106;
    for(let x=0;x<nb.length;x++)for(let y=x+1;y<nb.length;y++){const[a,oa]=nb[x],[b,ob]=nb[y];const la=bondLen(g.atoms[j].el,g.atoms[a].el,oa),lb=bondLen(g.atoms[j].el,g.atoms[b].el,ob);
      R.push([a,b,Math.sqrt(la*la+lb*lb-2*la*lb*Math.cos(th)),0.5]);pairKey.add(key(a,b))}}
  g.bonds.forEach(([b,c,o])=>{if(o!==2)return;const L=bondLen(g.atoms[b].el,g.atoms[c].el,2);
    for(const[a]of adj[b])if(a!==c)for(const[d]of adj[c])if(d!==b){if(pairKey.has(key(a,d)))continue;
      const la=V3.len(V3.sub(P[a],P[b])),ld=V3.len(V3.sub(P[d],P[c]));
      const A=[-la/2,la*.866],Dc=[L+ld/2,ld*.866],Dt=[L+ld/2,-ld*.866];const cis=Math.hypot(Dc[0]-A[0],Dc[1]-A[1]),trans=Math.hypot(Dt[0]-A[0],Dt[1]-A[1]);
      const cur=V3.len(V3.sub(P[a],P[d]));R.push([a,d,Math.abs(cur-cis)<Math.abs(cur-trans)?cis:trans,0.3]);pairKey.add(key(a,d))}});
  for(let it=0;it<1500;it++){const F=P.map(()=>[0,0,0]);
    for(const[a,b,d0,k]of R){const v=V3.sub(P[b],P[a]),d=V3.len(v)||1e-3,f=k*(d-d0)/d;F[a]=V3.add(F[a],V3.mul(v,f));F[b]=V3.sub(F[b],V3.mul(v,f))}
    for(let a=0;a<n;a++)for(let b=a+1;b<n;b++){if(pairKey.has(key(a,b)))continue;const v=V3.sub(P[b],P[a]),d=V3.len(v)||1e-3;
      const lim=(g.atoms[a].el==='H'||g.atoms[b].el==='H')?2.3:2.9;if(d<lim){const f=-0.15*(lim-d)/d;F[a]=V3.add(F[a],V3.mul(v,f));F[b]=V3.sub(F[b],V3.mul(v,f))}}
    const step=it<1000?0.12:0.06;for(let i=0;i<n;i++)P[i]=V3.add(P[i],V3.mul(F[i],step))}
  return P;
}
// ---------- PubChem lookup for organics ----------
function parseSDFfull(t){const L=t.split(/\r?\n/);const na=+L[3].slice(0,3),nb=+L[3].slice(3,6);const atoms=[],xyz=[],bonds=[];
  const CC={1:3,2:2,3:1,5:-1,6:-2,7:-3};
  for(let i=0;i<na;i++){const l=L[4+i];xyz.push([+l.slice(0,10),+l.slice(10,20),+l.slice(20,30)]);atoms.push({el:l.slice(31,34).trim(),fc:CC[+l.slice(36,39)]||0})}
  for(let i=0;i<nb;i++){const l=L[4+na+i];bonds.push([+l.slice(0,3)-1,+l.slice(3,6)-1,+l.slice(6,9)])}
  for(let i=4+na+nb;i<L.length;i++){if(L[i].startsWith('M  CHG')){const f=L[i].slice(6).trim().split(/\s+/).map(Number);for(let k=1;k<f.length;k+=2)atoms[f[k]-1].fc=f[k+1]}if(L[i].startsWith('M  END'))break}
  return {g:{atoms,bonds,heavy:atoms.filter(a=>a.el!=='H').length},xyz,cid:+L[0].trim()||null};
}
async function fetchT(url,ms=8000){const ctl=new AbortController();const t=setTimeout(()=>ctl.abort(),ms);try{const r=await fetch(url,{signal:ctl.signal});return r.ok?await r.text():null}catch(e){return null}finally{clearTimeout(t)}}
async function pubchemOrganic(q,ns){
  const base='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/';
  const t2=await fetchT(`${base}${ns}/${encodeURIComponent(q)}/SDF?record_type=2d`);if(!t2)return null;
  let s2;try{s2=parseSDFfull(t2)}catch(e){return null}
  let p3=null;if(s2.cid){const t3=await fetchT(`${base}cid/${s2.cid}/SDF?record_type=3d`);if(t3){try{const s3=parseSDFfull(t3);if(s3.g.atoms.length===s2.g.atoms.length&&s3.g.atoms.every((a,i)=>a.el===s2.g.atoms[i].el))p3=s3.xyz}catch(e){}}}
  let title=null;if(s2.cid){const tj=await fetchT(`${base}cid/${s2.cid}/property/Title/JSON`);try{title=JSON.parse(tj).PropertyTable.Properties[0].Title}catch(e){}}
  let[g,c2,c3]=heavyFirst(s2.g,s2.xyz.map(p=>[p[0],p[1]]),p3);
  // normalise 2D scale to bond length 1
  const bl=g.bonds.filter(([a,b])=>a<g.heavy&&b<g.heavy).map(([a,b])=>Math.hypot(c2[a][0]-c2[b][0],c2[a][1]-c2[b][1]));const m=bl.length?bl.reduce((x,y)=>x+y)/bl.length:1;
  c2=c2.map(p=>[p[0]/m,p[1]/m]);
  return {g,p2:c2.slice(0,g.heavy),p3:c3,title,cid:s2.cid};
}
// ---------- Organic analysis ----------
function hillOf(c){const els=Object.keys(c);const o=c.C?['C',...(c.H?['H']:[]),...els.filter(e=>e!=='C'&&e!=='H').sort()]:els.sort();return o.map(e=>e+(c[e]>1?c[e]:'')).join('')}
function orgAnalysis(g,p2,p3,name,src){
  const {on,det}=orgON(g);const n=g.atoms.length;
  const charge=g.atoms.reduce((s,a)=>s+a.fc,0);
  const tag=Array(n).fill(null);let cn=0;g.atoms.forEach((a,i)=>{if(a.el==='C')tag[i]='C'+(++cn)});
  const hCount=i=>g.bonds.filter(([a,b])=>(a===i&&g.atoms[b].el==='H')||(b===i&&g.atoms[a].el==='H')).length;
  const BN={1:'single',2:'double',3:'triple'};const desc=i=>{const parts=Object.values(det[i]).sort((x,y)=>y.o-x.o).map(d=>{const v=d.s*d.o*d.c;const nm=d.o===1?'':BN[d.o]+' ';return `${d.c} ${nm}bond${d.c>1?'s':''} to ${d.el} (${d.c>1?d.c+' × '+Q.of(d.s*d.o).str()+' = ':''}${Q.of(v).str()})`});
    if(g.atoms[i].fc)parts.push(`formal charge ${Q.of(g.atoms[i].fc).str()}`);return parts.length?parts.join(', '):'no bonds'};
  const slots=[],envSlot={};
  g.atoms.forEach((a,i)=>{if(a.el!=='C')return;const v=Q.of(on[i]);
    slots.push({key:tag[i],el:'C',label:tag[i],count:1,value:v,rule:'org',why:`${tag[i]} (${uniF('C'+(hCount(i)?'H'+(hCount(i)>1?hCount(i):''):''))}): ${desc(i)} → ${v.str()}.`});envSlot['a'+i+'|C']=tag[i]});
  const order=[];g.atoms.forEach(a=>{if(a.el!=='C'&&!order.includes(a.el))order.push(a.el)});
  for(const el of order){const idx=[...g.atoms.keys()].filter(i=>g.atoms[i].el===el);const same=idx.every(i=>on[i]===on[idx[0]]);
    const rf=v=>{const q=Q.of(v);if(el==='O'&&q.eq(-2))return 'O';if(el==='H'&&q.eq(1))return 'H';if(el==='H'&&q.n<0)return 'hydride';if((HALOGENS.includes(el)||el==='F')&&q.eq(-1))return el==='F'?'F':'X';if(isMetal(el)&&!idx.some(i=>Object.keys(det[i]).length))return monoRule(el,q,false);return 'org'};
    if(same){const v=Q.of(on[idx[0]]);const ds=[...new Set(idx.map(desc))];
      slots.push({key:el,el,label:el,count:idx.length,value:v,rule:rf(on[idx[0]]),why:`${idx.length>1?'Each '+el:el}: ${ds.join('; or ')} → ${v.str()}.`});idx.forEach(i=>envSlot['a'+i+'|'+el]=el)}
    else idx.forEach((i,k)=>{const key=el+(k+1);tag[i]=key;const v=Q.of(on[i]);slots.push({key,el,label:key,count:1,value:v,rule:rf(on[i]),why:`${key}: ${desc(i)} → ${v.str()}.`});envSlot['a'+i+'|'+el]=key})}
  slots.forEach(s=>s.fraction=false);
  const c=graphComp(g),hill=hillOf(c);
  const cs=slots.filter(s=>s.el==='C');let avg=new Q(0);cs.forEach(s=>avg=avg.add(s.value));avg=avg.div(cs.length||1);
  return {isOrg:true,name,src,body:hill,charge,html:htmlF(hill,charge),text:uniF(hill)+chargeSup(charge),slots,envSlot,avgC:avg,maxC:cs.reduce((m,s)=>Math.max(m,s.value.val()),-99),
    org:{g,p2,p3,tag,hCount:g.atoms.map((_,i)=>hCount(i))},parts:[{kind:'org'}]};
}
async function analyzeOrganic(input,allowNet=true){
  const q=String(input).trim();if(!q)throw chemErr('Please type a name or a SMILES string.');
  const db=ORG_DB.find(x=>x.name.toLowerCase()===q.toLowerCase()||x.smi===q||(x.alt||[]).includes(q.toLowerCase()));
  let g,p2=null,p3=null,src='built-in model',name=q;
  if(db){g=parseSMILES(db.smi);name=db.name}
  else{let parsed=null;try{parsed=parseSMILES(q)}catch(e){}
    if(parsed){const pc=allowNet?await pubchemOrganic(q,'smiles'):null;
      if(pc&&keyOf(graphComp(pc.g))===keyOf(graphComp(parsed))){({g,p2,p3}=pc);name=pc.title||q;src='PubChem CID '+pc.cid}else{g=parsed;src='approximate model'}}
    else{const pc=allowNet?await pubchemOrganic(q,'name'):null;
      if(!pc)throw chemErr(`I couldn't find "${q}". Check the spelling, or enter a SMILES string (e.g. CCO for ethanol). Looking up names needs an internet connection.`);
      ({g,p2,p3}=pc);name=pc.title||q;src='PubChem CID '+pc.cid}}
  if(!g.atoms.some(a=>a.el==='C'))throw chemErr('That molecule has no carbon. Use the Formula option for inorganic compounds.');
  if(g.atoms.length>120)throw chemErr('That molecule is too large for this tool.');
  if(!p2)p2=layout2D(g);if(!p3)p3=embed3D(g,p2);
  return orgAnalysis(g,p2,p3,name,src);
}
function orgScene(an){const o=an.org,{on}=orgON(o.g);
  const cen=o.p3.reduce((a,p)=>V3.add(a,p),[0,0,0]).map(v=>v/o.p3.length);
  return{atoms:o.g.atoms.map((a,i)=>({el:a.el,p:V3.sub(o.p3[i],cen),on:Q.of(on[i]),env:'a'+i,tag:a.el==='C'?o.tag[i]:null})),bonds:o.g.bonds.map(b=>b.slice()),missing:false,truncated:false,sources:[an.src]}}

const ORG_DB=[
 {smi:'C',name:'methane'},{smi:'CC',name:'ethane'},{smi:'CCC',name:'propane'},{smi:'CC(C)C',name:'2-methylpropane',alt:['isobutane']},
 {smi:'C=C',name:'ethene',alt:['ethylene']},{smi:'C=CC',name:'propene',alt:['propylene']},{smi:'C#C',name:'ethyne',alt:['acetylene']},
 {smi:'CO',name:'methanol'},{smi:'CCO',name:'ethanol'},{smi:'CC(O)C',name:'2-propanol',alt:['isopropanol','isopropyl alcohol']},{smi:'OCCO',name:'ethylene glycol'},
 {smi:'C=O',name:'formaldehyde',alt:['methanal']},{smi:'CC=O',name:'acetaldehyde',alt:['ethanal']},{smi:'O=Cc1ccccc1',name:'benzaldehyde'},
 {smi:'CC(=O)C',name:'acetone',alt:['propanone']},{smi:'OC=O',name:'formic acid',alt:['methanoic acid']},{smi:'CC(=O)O',name:'acetic acid',alt:['ethanoic acid']},
 {smi:'OC(=O)C(=O)O',name:'oxalic acid'},{smi:'CC(=O)OC',name:'methyl acetate'},{smi:'COC',name:'dimethyl ether'},{smi:'CCOCC',name:'diethyl ether'},
 {smi:'CN',name:'methylamine'},{smi:'CN(C)C',name:'trimethylamine'},{smi:'CC(N)=O',name:'acetamide'},{smi:'NC(N)=O',name:'urea'},{smi:'CC#N',name:'acetonitrile'},
 {smi:'CCl',name:'chloromethane'},{smi:'ClCCl',name:'dichloromethane'},{smi:'ClC(Cl)Cl',name:'chloroform',alt:['trichloromethane']},{smi:'CCBr',name:'bromoethane'},
 {smi:'c1ccccc1',name:'benzene'},{smi:'Cc1ccccc1',name:'toluene'},{smi:'Oc1ccccc1',name:'phenol'},{smi:'OC(=O)c1ccccc1',name:'benzoic acid'},{smi:'c1ccncc1',name:'pyridine'},
 {smi:'C[N+](=O)[O-]',name:'nitromethane'},{smi:'CC(=O)[O-]',name:'acetate ion'},{smi:'OCC(O)C(O)C(O)C(O)C=O',name:'glucose (open-chain form)',alt:['glucose']},
 {smi:'O=C=O',name:'carbon dioxide'}];
LEVELS[6]={name:'Organic compounds',desc:'Each carbon by bond counting, plus oxidation series'};
const OS=(s,n)=>[s,n,'org'];
PROBLEMS[6]=[...ORG_DB.filter(x=>!['O=C=O'].includes(x.smi)).map(x=>OS(x.smi,x.name)),
 {series:'Oxidation series: methane → carbon dioxide',items:[OS('C','methane'),OS('CO','methanol'),OS('C=O','formaldehyde'),OS('OC=O','formic acid'),OS('O=C=O','carbon dioxide')]},
 {series:'Oxidation series: ethane → acetic acid',items:[OS('CC','ethane'),OS('CCO','ethanol'),OS('CC=O','acetaldehyde'),OS('CC(=O)O','acetic acid')]}];
if(typeof module!=='undefined')Object.assign(module.exports,{parseSMILES,analyzeOrganic,orgON,orgScene,ORG_DB,layout2D});
