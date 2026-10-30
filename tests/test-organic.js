const C=require('./all.js');
(async()=>{let bad=0;
for(const x of C.ORG_DB){try{const a=await C.analyzeOrganic(x.smi,false);
 // sum check
 let t=new C.Q(0);a.slots.forEach(s=>t=t.add(s.value.mul(s.count)));
 const sc=C.orgScene(a);
 // geometry sanity: bond lengths and clashes
 const A=sc.atoms;let maxDev=0;sc.bonds.forEach(([i,j])=>{const d=Math.hypot(...[0,1,2].map(k=>A[i].p[k]-A[j].p[k]));maxDev=Math.max(maxDev,Math.abs(d-1.3))});
 let clash=0;const bs=new Set(sc.bonds.map(([i,j])=>i+'_'+j));for(let i=0;i<A.length;i++)for(let j=i+1;j<A.length;j++){if(bs.has(i+'_'+j)||bs.has(j+'_'+i))continue;const d=Math.hypot(...[0,1,2].map(k=>A[i].p[k]-A[j].p[k]));if(d<1.4)clash++}
 console.log(x.name.padEnd(26),a.text.padEnd(10),a.slots.map(s=>s.label+'='+s.value.str()).join(' '),'| sum',t.str(),'avgC',a.avgC.str(),'clash',clash);
 if(!t.eq(a.charge)||clash)bad++;
}catch(e){bad++;console.log('ERR',x.name,e.message)}}
console.log('bad',bad);
for(const s of ['c1ccc2ccccc2c1','c1cc[nH]c1','CCCCCCCC','C1CCCCC1','not a molecule','[Na+].CC(=O)[O-]']){try{const a=await C.analyzeOrganic(s,false);console.log(s,a.slots.map(x=>x.label+'='+x.value.str()).join(' '))}catch(e){console.log(s,'ERR',e.message)}}
const a=await C.analyzeOrganic('ethanol',false);a.slots.forEach(s=>console.log(s.label,s.rule,s.why));
})();
