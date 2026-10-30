const C=require('./all.js');
(async()=>{
let bad=0;
for(const lv of [1,2,3,4,5])for(const [f,n] of C.PROBLEMS[lv]){
 try{const a=C.analyze(f);
  const s=a.slots.map(s=>`${s.label}=${s.value.str()}[${s.rule}]`).join('  ');
  const sc=await C.buildScene(a,false);
  console.log(`L${lv} ${a.text.padEnd(22)} ${s}  | 3D:${sc?sc.atoms.length+' atoms'+(sc.missing?' MISSING':''):'NONE'}`);
  if(!sc||sc.missing)bad++;
 }catch(e){bad++;console.log('ERR',f,e.message)}
}
console.log('bad',bad);
for(const t of ['FeSO4','Fe(NO3)3','KAl(SO4)2','NaSCN','HCN','CuSO4.5H2O','SO4 2-','Cr2O7²⁻','Xx','Fe(NO3','[Co(NH3)5Br]SO4','KClO4','BF3','C2H6','NaBH4']){try{const a=C.analyze(t);console.log(t,'->',a.slots.map(s=>`${s.label}=${s.value.str()}`).join(' '))}catch(e){console.log(t,'ERR',e.message)}}
const a=C.analyze('Fe2(SO4)3');a.slots.forEach(s=>console.log(s.label,s.why));
const b=C.analyze('Na2O2');b.slots.forEach(s=>console.log(s.label,s.why));
const c=C.analyze('[Co(NH3)5Cl]Cl2');c.slots.forEach(s=>console.log(s.label,s.why));
for(const x of ['+3','-2','8/3','2.67','-0.5','3+','abc','',' -1/2']) console.log(JSON.stringify(x), C.parseAnswer(x));
})();
