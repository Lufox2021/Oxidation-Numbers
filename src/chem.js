// ===================== Oxidation Number Engine =====================
// ---------- Exact fractions ----------
function gcd(a,b){a=Math.abs(a);b=Math.abs(b);while(b){const t=a%b;a=b;b=t}return a||1}
class Q{
  constructor(n,d=1){if(d===0)throw new Error('Division by zero');if(d<0){n=-n;d=-d}const g=gcd(n,d);this.n=n/g;this.d=d/g}
  static of(x){return x instanceof Q?x:new Q(x,1)}
  add(o){o=Q.of(o);return new Q(this.n*o.d+o.n*this.d,this.d*o.d)}
  sub(o){o=Q.of(o);return new Q(this.n*o.d-o.n*this.d,this.d*o.d)}
  mul(o){o=Q.of(o);return new Q(this.n*o.n,this.d*o.d)}
  div(o){o=Q.of(o);return new Q(this.n*o.d,this.d*o.n)}
  neg(){return new Q(-this.n,this.d)}
  eq(o){o=Q.of(o);return this.n===o.n&&this.d===o.d}
  isInt(){return this.d===1}
  val(){return this.n/this.d}
  str(){if(this.n===0)return '0';const s=this.n>0?'+':'−';const a=Math.abs(this.n);return s+(this.d===1?a:a+'/'+this.d)}
}
function chemErr(m){const e=new Error(m);e.chem=true;return e}

// ---------- Elements: symbol group EN covalent-radius ----------
const ELDATA=`H 1 2.20 .31|He 18 0 .28|Li 1 .98 1.28|Be 2 1.57 .96|B 13 2.04 .84|C 14 2.55 .76|N 15 3.04 .71|O 16 3.44 .66|F 17 3.98 .57|Ne 18 0 .58|Na 1 .93 1.66|Mg 2 1.31 1.41|Al 13 1.61 1.21|Si 14 1.90 1.11|P 15 2.19 1.07|S 16 2.58 1.05|Cl 17 3.16 1.02|Ar 18 0 1.06|K 1 .82 2.03|Ca 2 1.00 1.76|Sc 3 1.36 1.70|Ti 4 1.54 1.60|V 5 1.63 1.53|Cr 6 1.66 1.39|Mn 7 1.55 1.39|Fe 8 1.83 1.32|Co 9 1.88 1.26|Ni 10 1.91 1.24|Cu 11 1.90 1.32|Zn 12 1.65 1.22|Ga 13 1.81 1.22|Ge 14 2.01 1.20|As 15 2.18 1.19|Se 16 2.55 1.20|Br 17 2.96 1.20|Kr 18 3.00 1.16|Rb 1 .82 2.20|Sr 2 .95 1.95|Y 3 1.22 1.90|Zr 4 1.33 1.75|Nb 5 1.6 1.64|Mo 6 2.16 1.54|Tc 7 1.9 1.47|Ru 8 2.2 1.46|Rh 9 2.28 1.42|Pd 10 2.20 1.39|Ag 11 1.93 1.45|Cd 12 1.69 1.44|In 13 1.78 1.42|Sn 14 1.96 1.39|Sb 15 2.05 1.39|Te 16 2.1 1.38|I 17 2.66 1.39|Xe 18 2.6 1.40|Cs 1 .79 2.44|Ba 2 .89 2.15|La 3 1.1 2.07|Hf 4 1.3 1.75|Ta 5 1.5 1.70|W 6 2.36 1.62|Re 7 1.9 1.51|Os 8 2.2 1.44|Ir 9 2.2 1.41|Pt 10 2.28 1.36|Au 11 2.54 1.36|Hg 12 2.0 1.32|Tl 13 1.62 1.45|Pb 14 2.33 1.46|Bi 15 2.02 1.48`;
const EL={};ELDATA.split('|').forEach(r=>{const[s,g,en,rc]=r.trim().split(/\s+/);EL[s]={sym:s,group:+g,en:+en,rc:+rc}});
const NONMETALS=new Set(['H','He','B','C','N','O','F','Ne','Si','P','S','Cl','Ar','Ge','As','Se','Br','Kr','Sb','Te','I','Xe']);
const isMetal=e=>!NONMETALS.has(e);
const HALOGENS=['Cl','Br','I'];
function valenceE(e){if(e==='H')return 1;const g=EL[e].group;return g<=12?g:g-10}
function typicalAnion(e){if(e==='H')return -1;const g=EL[e].group;return g>=14&&g<=17?g-18:null}

// ---------- Text helpers ----------
const SUBD='₀₁₂₃₄₅₆₇₈₉',SUPD='⁰¹²³⁴⁵⁶⁷⁸⁹';
function uniF(f){return f.replace(/([A-Za-z\)\]])(\d+)/g,(m,a,d)=>a+[...d].map(c=>SUBD[+c]).join(''))}
function chargeSup(q){q=Q.of(q);if(q.n===0)return '';const a=Math.abs(q.n);return (a===1?'':[...String(a)].map(c=>SUPD[+c]).join(''))+(q.n>0?'⁺':'⁻')}
function chargeTxt(q){q=Q.of(q);if(q.n===0)return '0';const a=Math.abs(q.n);return (a===1?'':a)+(q.n>0?'+':'−')}
function htmlF(f,q){let s=f.replace(/([A-Za-z\)\]])(\d+)/g,'$1<sub>$2</sub>');if(q&&Q.of(q).n!==0)s+='<sup>'+chargeTxt(q)+'</sup>';return s}

// ---------- Input normalization ----------
function normalizeInput(raw){
  let s=String(raw).trim().replace(/[−–—]/g,'-').replace(/⁺/g,'+').replace(/⁻/g,'-');
  s=s.replace(/[₀-₉]/g,c=>SUBD.indexOf(c));
  s=s.replace(/[⁰-⁹¹²³]+/g,m=>'^'+[...m].map(c=>SUPD.indexOf(c)).join(''));
  if(!s)throw chemErr('Please type a formula.');
  let charge=0,m;
  if((m=s.match(/\^\s*(\d*)\s*([+-])\s*$/))){charge=(m[2]==='-'?-1:1)*(m[1]?+m[1]:1);s=s.slice(0,m.index)}
  else if((m=s.match(/\^\s*([+-])\s*(\d*)\s*$/))){charge=(m[1]==='-'?-1:1)*(m[2]?+m[2]:1);s=s.slice(0,m.index)}
  else if((m=s.match(/\s+(\d*)([+-])\s*$/))){charge=(m[2]==='-'?-1:1)*(m[1]?+m[1]:1);s=s.slice(0,m.index)}
  else if((m=s.match(/([+-]+)\s*$/))){charge=[...m[1]].reduce((a,c)=>a+(c==='+'?1:-1),0);s=s.slice(0,m.index)}
  s=s.replace(/\s+/g,'');
  if(!s)throw chemErr('Please type a formula.');
  return {body:s,charge};
}

// ---------- Tokenizer ----------
const ABBR={en:{C:2,H:8,N:2},ox:{C:2,O:4}};
function readNum(s,i){let j=i;while(j<s.length&&/\d/.test(s[j]))j++;return [j>i?+s.slice(i,j):1,j]}
function readSym(s,i){const two=s.substr(i,2);if(/^[A-Z][a-z]$/.test(two)&&EL[two])return two;if(EL[s[i]])return s[i];throw chemErr(`"${/^[A-Z][a-z]$/.test(two)?two:s[i]}" is not an element symbol I know.`)}
function tokenize(s){
  const t=[];let i=0;
  while(i<s.length){const c=s[i];
    if(c==='('||c==='['){let d=0,j=i;for(;j<s.length;j++){if(s[j]==='('||s[j]==='[')d++;else if(s[j]===')'||s[j]===']'){d--;if(d===0)break}}
      if(j>=s.length)throw chemErr('A bracket or parenthesis is not closed.');
      const inner=s.slice(i+1,j);if(!inner)throw chemErr('Empty parentheses.');const[n,k]=readNum(s,j+1);
      t.push({type:c==='('?'group':'bracket',inner,count:n,raw:s.slice(i,k)});i=k}
    else if(s.startsWith('en',i)||s.startsWith('ox',i)){const ab=s.substr(i,2);const[n,k]=readNum(s,i+2);t.push({type:'abbr',abbr:ab,count:n,raw:s.slice(i,k)});i=k}
    else if(/[A-Z]/.test(c)){const sym=readSym(s,i);const[n,k]=readNum(s,i+sym.length);t.push({type:'el',el:sym,count:n,raw:s.slice(i,k)});i=k}
    else if(c===')'||c===']')throw chemErr('A closing bracket has no matching opening bracket.');
    else throw chemErr(`I don't understand the character "${c}".`);
  }
  return t;
}
function comp(tokens,mult=1,acc={}){
  for(const t of tokens){const m=mult*t.count;
    if(t.type==='el')acc[t.el]=(acc[t.el]||0)+m;
    else if(t.type==='abbr'){for(const[e,n]of Object.entries(ABBR[t.abbr]))acc[e]=(acc[e]||0)+n*m}
    else comp(tokenize(t.inner),m,acc)}
  return acc;
}
const compOf=f=>comp(tokenize(f));
function keyOf(c){return Object.keys(c).sort().map(e=>e+(c[e]>1?c[e]:'')).join('')}
function atomCount(c){return Object.values(c).reduce((a,b)=>a+b,0)}

// ---------- Known ions & ligands ----------
const ANION_LIST=[['OH','hydroxide',-1],['NO3','nitrate',-1],['NO2','nitrite',-1],['SO4','sulfate',-2],['SO3','sulfite',-2],['HSO4','hydrogen sulfate',-1],['HSO3','hydrogen sulfite',-1],['PO4','phosphate',-3],['HPO4','hydrogen phosphate',-2],['H2PO4','dihydrogen phosphate',-1],['CO3','carbonate',-2],['HCO3','hydrogen carbonate',-1],['ClO4','perchlorate',-1],['ClO3','chlorate',-1],['ClO2','chlorite',-1],['ClO','hypochlorite',-1],['OCl','hypochlorite',-1],['BrO3','bromate',-1],['BrO','hypobromite',-1],['IO3','iodate',-1],['IO4','periodate',-1],['MnO4','permanganate',-1],['CrO4','chromate',-2],['Cr2O7','dichromate',-2],['C2O4','oxalate',-2],['CN','cyanide',-1],['SCN','thiocyanate',-1],['S2O3','thiosulfate',-2],['S4O6','tetrathionate',-2],['N3','azide',-1],['AsO4','arsenate',-3],['SiO3','silicate',-2],['BO3','borate',-3],['AlH4','tetrahydridoaluminate',-1],['BH4','tetrahydridoborate',-1],['O2','peroxide',-2]];
const ANIONS={};ANION_LIST.forEach(([f,name,q])=>{if(f==='O2')return;const k=keyOf(compOf(f));if(!ANIONS[k])ANIONS[k]={f,name,charge:q}});
const CATIONS={};[['NH4','ammonium',1],['H3O','hydronium',1]].forEach(([f,name,q])=>{CATIONS[keyOf(compOf(f))]={f,name,charge:q}});
const LIGAND_LIST=[
  ['NH3','ammine',0,'N(H)(H)H','mono'],['H2O','aqua',0,'O(H)H','mono'],['CO','carbonyl',0,'[C-]#[O+]','mono'],
  ['CN','cyano',-1,'[C-]#N','mono'],['OH','hydroxo',-1,'[O-]H','mono'],['Cl','chloro',-1,'[Cl-]','mono'],['Br','bromo',-1,'[Br-]','mono'],
  ['F','fluoro',-1,'[F-]','mono'],['I','iodo',-1,'[I-]','mono'],['NO2','nitro',-1,'[N+](=O)[O-]','mono'],['SCN','thiocyanato',-1,'N=C=[S-]','mono'],
  ['C2O4','oxalato (ox)',-2,null,'ox'],['C2H8N2','ethylenediamine (en)',0,null,'en']];
const LIGANDS={};LIGAND_LIST.forEach(([f,name,q,smi,type])=>{LIGANDS[keyOf(compOf(f))]={f:f==='C2H8N2'?'en':f==='C2O4'?'C2O4':f,name,charge:q,smi,type}});

// ---------- Rules ----------
const RULES={
  1:{title:'Free elements',text:'An atom in an uncombined (free) element has an oxidation number of 0.',ex:'Na, Fe, O₂, Cl₂, P₄, S₈ → every atom is 0'},
  2:{title:'Monatomic ions',text:'For a monatomic ion, the oxidation number equals the charge of the ion.',ex:'Fe³⁺ → +3;  S²⁻ → −2;  in Na₂S, S = −2'},
  3:{title:'Sum rule',text:'The oxidation numbers of all atoms add up to 0 in a neutral compound, or to the charge of the ion for an ion.',ex:'SO₄²⁻: S + 4(−2) = −2 → S = +6'},
  4:{title:'Group 1 metals',text:'Group 1 (alkali) metals are always +1 in compounds.',ex:'Na in NaCl → +1;  K in KMnO₄ → +1'},
  5:{title:'Group 2 metals and aluminum',text:'Group 2 (alkaline earth) metals are always +2 in compounds; aluminum is always +3.',ex:'Mg in MgO → +2;  Al in Al₂O₃ → +3'},
  6:{title:'Fluorine',text:'Fluorine is always −1 in compounds (it is the most electronegative element).',ex:'F in HF, OF₂, SF₆ → −1'},
  7:{title:'Hydrogen',text:'Hydrogen is +1 when combined with nonmetals and −1 when combined with metals (metal hydrides).',ex:'H in HCl, H₂O → +1;  H in NaH, CaH₂ → −1'},
  8:{title:'Oxygen',text:'Oxygen is usually −2. Exceptions: peroxides (−1), superoxides (−½), and compounds with fluorine (positive, e.g. +2 in OF₂).',ex:'H₂O → −2;  H₂O₂ → −1;  KO₂ → −1/2;  OF₂ → +2'},
  9:{title:'Chlorine, bromine, iodine',text:'Cl, Br and I are usually −1, except when combined with oxygen or with a more electronegative halogen; then they are positive.',ex:'NaCl → −1;  ClO₄⁻ → Cl = +7;  ICl → I = +1'},
  10:{title:'More electronegative element',text:'In a compound of nonmetals not covered above, the more electronegative element takes its usual negative value (its charge as a monatomic anion: N −3, S −2, C −4 …).',ex:'CS₂ → S = −2, C = +4;  CN⁻ → N = −3, C = +2'},
  11:{title:'Ionic compounds with polyatomic ions or variable-charge metals',text:'Split the compound into its ions. Use the known charge of the anion (or polyatomic ion) to find the metal\'s charge, then treat each polyatomic ion with the sum rule.',ex:'Fe₂(SO₄)₃: 3 SO₄²⁻ = −6 → each Fe = +3, and S = +6'},
  12:{title:'Coordination compounds',text:'Find the charge of the complex ion from the counter-ions. Subtract the total ligand charge to get the metal. Neutral ligands: NH₃, H₂O, CO, en. Anionic: Cl⁻, Br⁻, CN⁻, OH⁻, NO₂⁻, SCN⁻ (−1), C₂O₄²⁻ (−2).',ex:'K₄[Fe(CN)₆]: complex = −4; 6 CN⁻ = −6 → Fe = +2'},
  14:{title:'Carbon in organic molecules',text:'Find each carbon separately by counting its bonds: +1 for each bond to a more electronegative atom (O, N, S, halogens; a double bond counts 2), −1 for each bond to H, and 0 for each bond to another carbon. Add any formal charge. The same bond counting works for every atom in the molecule.',ex:'Ethanol CH₃–CH₂–OH: CH₃ = 3(−1) = −3;  CH₂OH = 2(−1) + 1(+1) = −1;  average C = −2'},
  13:{title:'Fractional oxidation numbers',text:'When identical atoms are in different environments, the sum rule gives their average oxidation number, which can be a fraction.',ex:'Fe₃O₄ → Fe = +8/3 (one Fe²⁺ and two Fe³⁺);  S₄O₆²⁻ → S = +5/2'}
};
const RULE_NUM={free:1,monatomic:2,sum:3,g1:4,g2:5,al:5,F:6,H:7,hydride:7,O:8,peroxide:8,superoxide:8,OF:8,X:9,Xpos:9,EN:10,ionic:11,complex:12,org:14};
function hintFor(slot,species){
  const e=slot.el,ch=Q.of(species.charge);
  switch(slot.rule){
    case 'free':return `Is ${e} combined with any other element here?`;
    case 'monatomic':return `${e} is present as a simple (monatomic) ion. What charge does that ion carry?`;
    case 'sum':return `Assign the other atoms first using the rules you're sure of. Then all oxidation numbers must add up to ${ch.n===0?'0 (the compound is neutral)':ch.str()+' (the charge on the ion or unit)'}.`;
    case 'g1':case 'g2':case 'al':return `Where is ${e} on the periodic table?`;
    case 'F':return 'Fluorine is the most electronegative element. What does that mean for its oxidation number?';
    case 'H':case 'hydride':return 'Is hydrogen combined with nonmetals here, or only with a metal?';
    case 'O':return 'Oxygen has a usual oxidation number. Is there anything in this formula that would change it?';
    case 'peroxide':case 'superoxide':case 'OF':return 'Assign everything else first, then check: can oxygen really be −2 here?';
    case 'X':case 'Xpos':return `Halogens are usually negative. Is ${e} combined with oxygen or with a more electronegative halogen?`;
    case 'EN':return `Which element is more electronegative? That one takes its usual negative value.`;
    case 'ionic':return `Identify the anion and its charge first. The ${e} ions must balance that charge.`;
    case 'org':return `Look at every bond from ${slot.label}. Each bond to a more electronegative atom counts +1 (a double bond counts +2), each bond to a less electronegative atom (like H, for carbon) counts −1, and each bond to the same element counts 0.`;
    case 'complex':return 'Find the charge of the complex ion (in the brackets) from the counter-ions, then subtract the charges of the ligands.';
  }
  return 'Use the sum rule.';
}

// ---------- Solve one neutral molecule or one ion ----------
function monoRule(e,v,multi){
  if(e==='O'){if(v.eq(-2))return 'O';if(v.eq(-1))return 'peroxide';if(v.eq(new Q(-1,2)))return 'superoxide';return 'sum'}
  if(e==='H')return v.n<0?'hydride':'H';
  if(e==='F')return 'F';
  if(HALOGENS.includes(e))return v.eq(-1)?'X':'Xpos';
  const g=EL[e].group;
  if(g===1&&v.eq(1))return 'g1';if(g===2&&v.eq(2))return 'g2';if(e==='Al'&&v.eq(3))return 'al';
  return multi?'sum':'monatomic';
}
function solvePart(c,charge,prefix=''){
  const els=Object.keys(c),total=Q.of(charge),res={};
  const set=(e,v,rule,why)=>{res[e]={value:Q.of(v),rule,why:prefix+why}};
  if(els.length===1){const e=els[0],n=c[e],v=total.div(n);
    if(total.n===0)set(e,0,'free',`${e} is not combined with any other element, so it is a free element: 0.`);
    else if(n===1)set(e,v,monoRule(e,v,false),`${e} is a monatomic ion with charge ${total.str()}, so its oxidation number is ${v.str()}.`);
    else set(e,v,monoRule(e,v,true),`${prefix?'':'The '+uniF(e+n)+' unit carries a charge of '+total.str()+'. '}Shared over ${n} identical atoms: ${n} × ${e} = ${total.str()} → ${e} = ${v.str()}.`);
    return res;}
  const unknown=new Set(els),has=e=>e in c;
  const assign=(e,v,rule,why)=>{set(e,v,rule,why);unknown.delete(e)};
  if(has('F'))assign('F',-1,'F','Fluorine is always −1 in compounds.');
  for(const e of els){if(e==='H')continue;const g=EL[e].group;
    if(g===1)assign(e,1,'g1',`${e} is a Group 1 metal: always +1 in compounds.`);
    else if(g===2)assign(e,2,'g2',`${e} is a Group 2 metal: always +2 in compounds.`);
    else if(e==='Al')assign(e,3,'al','Aluminum is always +3 in compounds.');}
  if(unknown.has('H')){const others=els.filter(e=>e!=='H');
    if(others.every(e=>isMetal(e)||e==='B'))assign('H',-1,'hydride',`H is combined only with a metal (${others.join(', ')}), so it is a hydride: −1.`);
    else assign('H',1,'H','H is combined with nonmetals here, so it is +1.');}
  if(unknown.has('O')&&unknown.size>1)assign('O',-2,'O','Oxygen takes its usual value, −2 (no peroxide, superoxide or fluorine here).');
  for(const X of HALOGENS){if(!unknown.has(X)||unknown.size<=1||has('O')||has('F'))continue;
    if(HALOGENS.some(Y=>Y!==X&&has(Y)&&EL[Y].en>EL[X].en))continue;
    assign(X,-1,'X',`${X} is combined with neither oxygen nor a more electronegative halogen, so it is −1.`);}
  while(unknown.size>1){const cand=[...unknown].filter(e=>!isMetal(e)&&typicalAnion(e)!==null).sort((a,b)=>EL[b].en-EL[a].en);
    if(!cand.length)break;const e=cand[0],others=[...unknown].filter(x=>x!==e),t=typicalAnion(e);
    assign(e,t,'EN',`${e} is more electronegative than ${others.join(' and ')}, so it takes its usual negative value, ${Q.of(t).str()}.`);}
  if(unknown.size>1)throw chemErr(`The rules can't pin down ${[...unknown].join(' and ')} separately. If the compound contains a polyatomic ion, try writing it in parentheses, e.g. Fe(NO3)3.`);
  let known=new Q(0);const terms=[];
  for(const e of els)if(!unknown.has(e)){known=known.add(res[e].value.mul(c[e]));terms.push(`${c[e]>1?c[e]+'('+res[e].value.str()+')':res[e].value.str()} for ${e}`)}
  if(unknown.size===1){const e=[...unknown][0],n=c[e],x=total.sub(known).div(n);
    let rule='sum';
    if(e==='O'){if(has('F')&&x.n>-2)rule='OF';else rule=monoRule('O',x,true)}
    else if(HALOGENS.includes(e))rule=x.eq(-1)?'X':'Xpos';
    const tgt=total.n===0?'0 (neutral)':total.str()+' (the charge)';
    set(e,x,rule,`Known: ${terms.join(', ')} → total ${known.str()}. Everything must add to ${tgt}, so ${n>1?n+' × ':''}${e} = ${total.str()} − (${known.str()}) = ${total.sub(known).str()}${n>1?' → '+e+' = '+x.str():''}.`);}
  else if(!known.eq(total))throw chemErr(`Using the rules gives a total of ${known.str()}, but the overall charge is ${total.str()}. Check the formula and the charge.`);
  return res;
}

// ---------- Species analysis ----------
// entry: {el,count,env,envLabel,value,rule,why}
// part (3D): {kind:'unit',comp,charge,mult,sol:{el:Q},env} | {kind:'complex',...}
function addEntries(entries,c,mult,sol,env,envLabel){for(const e of Object.keys(c))entries.push({el:e,count:c[e]*mult,env,envLabel,value:sol[e].value,rule:sol[e].rule,why:sol[e].why})}
function solToVals(sol){const o={};for(const e in sol)o[e]=sol[e].value;return o}
function ionLabel(f,q){return uniF(f)+chargeSup(q)}

function tryIonic(tokens){
  if(tokens.length<2)return null;
  let cat,rest;const t0=tokens[0];
  if(t0.type==='group'){const cc=compOf(t0.inner),k=keyOf(cc);if(!CATIONS[k])return null;cat={comp:cc,charge:CATIONS[k].charge,fixed:true,mult:t0.count,label:CATIONS[k].f,poly:true,name:CATIONS[k].name};rest=tokens.slice(1)}
  else if(t0.type==='el'&&t0.el==='N'&&t0.count===1&&tokens[1].type==='el'&&tokens[1].el==='H'&&tokens[1].count===4&&tokens.length>2){cat={comp:{N:1,H:4},charge:1,fixed:true,mult:1,label:'NH4',poly:true,name:'ammonium'};rest=tokens.slice(2)}
  else if(t0.type==='el'&&isMetal(t0.el)){const g=EL[t0.el].group,fx=g===1?1:g===2?2:t0.el==='Al'?3:null;cat={comp:{[t0.el]:1},charge:fx,fixed:fx!==null,mult:t0.count,label:t0.el,poly:false};rest=tokens.slice(1)}
  else return null;
  if(!rest.length)return null;
  const restComp=comp(rest);if(Object.keys(restComp).some(isMetal)){const rk=rest.length===1&&rest[0].type==='group'?keyOf(compOf(rest[0].inner)):keyOf(restComp);if(!ANIONS[rk])return null}
  let an;
  if(rest.length===1&&rest[0].type==='group')an={comp:compOf(rest[0].inner),mult:rest[0].count};
  else an={comp:restComp,mult:1};
  if(Object.keys(an.comp).length===1){const e=Object.keys(an.comp)[0];an={single:true,el:e,count:an.comp[e]*an.mult}}
  const out={cation:cat,anions:[]};let anTotal;
  if(an.single){const t=typicalAnion(an.el);
    if(cat.fixed){const tot=-cat.charge*cat.mult,per=new Q(tot,an.count);
      if(t!==null&&per.eq(t))out.anions.push({comp:{[an.el]:1},charge:Q.of(t),mult:an.count,label:an.el,fromCation:true});
      else out.anions.push({comp:{[an.el]:an.count},charge:Q.of(tot),mult:1,label:an.el+an.count,fromCation:true,unit:true});
      anTotal=Q.of(tot);}
    else{if(t===null)return null;out.anions.push({comp:{[an.el]:1},charge:Q.of(t),mult:an.count,label:an.el,typical:true});anTotal=Q.of(t*an.count)}
  }else{const info=ANIONS[keyOf(an.comp)];
    if(cat.fixed){const per=new Q(-cat.charge*cat.mult,an.mult);out.anions.push({comp:an.comp,charge:per,mult:an.mult,label:info?info.f:keyOf(an.comp),info,fromCation:true,poly:true});anTotal=per.mul(an.mult)}
    else{if(!info)return null;out.anions.push({comp:an.comp,charge:Q.of(info.charge),mult:an.mult,label:info.f,info,poly:true});anTotal=Q.of(info.charge*an.mult)}
  }
  cat.q=cat.fixed?Q.of(cat.charge):anTotal.neg().div(cat.mult);
  out.anTotal=anTotal;return out;
}

function analyzeSimple(body,charge,entries,parts,notes,envBase){
  const tokens=tokenize(body);
  if(tokens.some(t=>t.type==='bracket'))throw chemErr('Square brackets are only for coordination compounds, e.g. K4[Fe(CN)6].');
  const sp=charge===0?tryIonic(tokens):null;
  if(sp){
    const cat=sp.cation;
    const catEnv=envBase+'cat',catLabel=cat.poly?ionLabel(cat.label,cat.charge):cat.label+' (cation)';
    let catSol;
    if(cat.poly)catSol=solvePart(cat.comp,cat.charge,`In the ${cat.name} ion ${ionLabel(cat.label,cat.charge)}: `);
    else{const e=cat.label;const s=solvePart({[e]:1},cat.q);
      if(cat.fixed){const g=EL[e].group;s[e].why=e==='Al'?'Aluminum is always +3 in compounds.':`${e} is a Group ${g} metal: always +${g} in compounds.`}
      if(!cat.fixed){const a=sp.anions[0];const anDesc=a.poly?`${a.mult>1?a.mult+' × ':''}${ionLabel(a.label,a.charge)}${a.info?' ('+a.info.name+')':''}`:`${a.mult} ${ionLabel(a.label,a.charge)}`;
        s[e]={value:cat.q,rule:'ionic',why:`The anions are ${anDesc}, totaling ${sp.anTotal.str()}. The ${cat.mult>1?cat.mult+' '+e+' ions':e+' ion'} must total ${sp.anTotal.neg().str()}, so each ${e} = ${cat.q.str()}.`};}
      catSol=s;}
    addEntries(entries,cat.comp,cat.mult,catSol,catEnv,catLabel);
    parts.push({kind:'unit',comp:cat.comp,charge:cat.poly?Q.of(cat.charge):cat.q,mult:cat.mult,sol:solToVals(catSol),env:catEnv});
    sp.anions.forEach((a,i)=>{
      const env=envBase+'an'+i;let prefix='';
      if(a.fromCation){const tot=a.charge.mul(a.mult).neg();prefix=`The ${cat.mult>1?cat.mult+' ':''}${cat.poly?ionLabel(cat.label,cat.charge):cat.label} contribute${cat.mult>1?'':'s'} ${tot.str()}, so ${a.mult>1?'each ':'the '}${uniF(a.label)} ${a.unit||a.poly?'unit':'ion'} carries ${a.charge.str()}. `}
      else if(a.info)prefix=`${ionLabel(a.label,a.charge)} is the ${a.info.name} ion. `;
      const sol=solvePart(a.comp,a.charge,a.poly||a.unit?prefix:(a.fromCation?prefix:''));
      const lab=a.poly?ionLabel(a.label,a.charge):a.label+' (anion)';
      addEntries(entries,a.comp,a.mult,sol,env,lab);
      parts.push({kind:'unit',comp:a.comp,charge:a.charge,mult:a.mult,sol:solToVals(sol),env});
    });
  }else{
    const c=comp(tokens);const sol=solvePart(c,charge);const env=envBase+'main';
    addEntries(entries,c,1,sol,env,uniF(body)+chargeSup(charge));
    parts.push({kind:'unit',comp:c,charge:Q.of(charge),mult:1,sol:solToVals(sol),env});
  }
}

function analyzeComplex(body,charge,entries,parts,notes){
  const m=body.match(/^([^\[]*)\[(.+)\](\d*)([^\]]*)$/);
  if(!m)throw chemErr('Write coordination compounds with one set of square brackets, e.g. [Co(NH3)5Cl]Cl2 or K4[Fe(CN)6].');
  const[,pre,inner,cnum,suf]=m;const cmult=cnum?+cnum:1;
  let outside=new Q(0);const outsideDesc=[];
  // counter-cations before the bracket
  if(pre){const tk=tokenize(pre);for(const t of tk){
    if(t.type==='el'&&isMetal(t.el)&&[1,2].includes(EL[t.el].group)){const q=EL[t.el].group;const sol=solvePart({[t.el]:1},q);sol[t.el].rule=q===1?'g1':'g2';sol[t.el].why=`${t.el} is a Group ${q} metal counter-ion: +${q}.`;
      addEntries(entries,{[t.el]:1},t.count,sol,'pre'+t.el,t.el+' (counter-ion)');parts.push({kind:'unit',comp:{[t.el]:1},charge:Q.of(q),mult:t.count,sol:solToVals(sol),env:'pre'+t.el});
      outside=outside.add(q*t.count);outsideDesc.push(`${t.count} ${t.el}⁺${q===2?'²':''}`.replace('⁺²','²⁺'));}
    else if(t.type==='group'&&CATIONS[keyOf(compOf(t.inner))]){const ci=CATIONS[keyOf(compOf(t.inner))];const cc=compOf(t.inner);const sol=solvePart(cc,ci.charge,`In ${ionLabel(ci.f,ci.charge)}: `);
      addEntries(entries,cc,t.count,sol,'pre'+ci.f,ionLabel(ci.f,ci.charge));parts.push({kind:'unit',comp:cc,charge:Q.of(ci.charge),mult:t.count,sol:solToVals(sol),env:'pre'+ci.f});
      outside=outside.add(ci.charge*t.count);outsideDesc.push(`${t.count} ${ionLabel(ci.f,ci.charge)}`);}
    else throw chemErr(`I don't recognize the counter-ion "${t.raw}" before the bracket. Use Group 1/2 metals or (NH4).`);}}
  // counter-anions after the bracket
  if(suf){const tk=tokenize(suf);let ac,amult=1,mono=false;
    if(tk.length===1&&tk[0].type==='el'){ac={[tk[0].el]:1};amult=tk[0].count;mono=true}
    else if(tk.length===1&&tk[0].type==='group'){ac=compOf(tk[0].inner);amult=tk[0].count}
    else ac=comp(tk);
    let q,lab,info;
    if(mono){const e=Object.keys(ac)[0];q=typicalAnion(e);if(q===null)throw chemErr(`I don't recognize the counter-ion "${suf}".`);lab=e}
    else{info=ANIONS[keyOf(ac)];if(!info)throw chemErr(`I don't recognize the counter-ion "${suf}" after the bracket.`);q=info.charge;lab=info.f}
    const sol=solvePart(ac,q,info?`${ionLabel(lab,q)} is the ${info.name} ion. `:'');
    if(mono)sol[lab].why=`${lab} is a counter-ion, ${ionLabel(lab,q)}: ${Q.of(q).str()}.`;
    addEntries(entries,ac,amult,sol,'suf',ionLabel(lab,q)+' (counter-ion)');parts.push({kind:'unit',comp:ac,charge:Q.of(q),mult:amult,sol:solToVals(sol),env:'suf'});
    outside=outside.add(q*amult);outsideDesc.push(`${amult} ${ionLabel(lab,q)}`);}
  const cq=Q.of(charge).sub(outside).div(cmult);
  if(!cq.isInt())throw chemErr('The charges of the counter-ions don\'t give a whole-number charge for the complex. Check the formula.');
  // inside
  const tk=tokenize(inner);if(!tk.length||tk[0].type!=='el'||!isMetal(tk[0].el)||tk[0].count!==1)throw chemErr('Start the bracket with the metal, e.g. [Co(NH3)6].');
  const metal=tk[0].el;const ligs=[];
  for(const t of tk.slice(1)){let c;
    if(t.type==='group'){const it=tokenize(t.inner);c=comp(it)}else if(t.type==='abbr')c={...ABBR[t.abbr]};else if(t.type==='el')c={[t.el]:1};else throw chemErr('Nested brackets are not supported.');
    const L=LIGANDS[keyOf(c)];if(!L)throw chemErr(`I don't recognize the ligand "${t.raw.replace(/\d+$/,'')}". Known ligands: NH3, H2O, CO, en, Cl, Br, F, I, CN, OH, NO2, SCN, C2O4 (ox).`);
    const ex=ligs.find(l=>l.L===L);if(ex)ex.count+=t.count;else ligs.push({L,comp:c,count:t.count});}
  if(!ligs.length)throw chemErr('No ligands found inside the bracket.');
  let ligTot=new Q(0);const ligDesc=[];
  for(const l of ligs){ligTot=ligTot.add(l.L.charge*l.count);ligDesc.push(`${l.count} ${l.L.f==='en'?'en':ionLabel(l.L.f,l.L.charge)} (${l.L.charge===0?'neutral':Q.of(l.L.charge*l.count).str()})`)}
  const mv=cq.sub(ligTot);
  const cqWhy=outsideDesc.length?`Counter-ions: ${outsideDesc.join(' + ')} = ${outside.str()}, so ${cmult>1?'each of the '+cmult+' complexes':'the complex ion'} carries ${cq.str()}.`:`The complex ion carries the overall charge, ${cq.str()}.`;
  const msol={[metal]:{value:mv,rule:'complex',why:`${cqWhy} Ligands: ${ligDesc.join(', ')} → ligand total ${ligTot.str()}. Metal = ${cq.str()} − (${ligTot.str()}) = ${mv.str()}.`}};
  addEntries(entries,{[metal]:1},cmult,msol,'metal',metal+' (metal center)');
  const cpart={kind:'complex',metal,metalValue:mv,mult:cmult,charge:cq,env:'metal',ligands:[]};
  for(const l of ligs){const env='lig'+l.L.f;const sol=solvePart(l.comp,l.L.charge,`In the ${l.L.name} ligand ${l.L.f==='en'?'(en, H₂NCH₂CH₂NH₂)':ionLabel(l.L.f,l.L.charge)}: `);
    addEntries(entries,l.comp,l.count*cmult,sol,env,(l.L.f==='en'?'en':ionLabel(l.L.f,l.L.charge))+' ligand');
    cpart.ligands.push({L:l.L,count:l.count,sol:solToVals(sol),env});}
  parts.push(cpart);
}

function analyze(raw){
  const {body,charge}=normalizeInput(raw);
  const pieces=body.split(/[·•*.]/).filter(Boolean);
  if(pieces.length>1&&charge!==0)throw chemErr('Hydrates should be neutral.');
  const entries=[],parts=[],notes=[];
  pieces.forEach((p,i)=>{
    if(i===0){if(p.includes('[')){analyzeComplex(p,charge,entries,parts,notes);const mv=(arr,f)=>{const a=arr.filter(f),b=arr.filter(x=>!f(x));arr.length=0;arr.push(...b,...a)};mv(entries,e=>e.env==='suf');mv(parts,x=>x.env==='suf')}else analyzeSimple(p,charge,entries,parts,notes,'')}
    else{const m=p.match(/^(\d*)(.+)$/);const k=m[1]?+m[1]:1;const sub=[],sp=[];analyzeSimple(m[2],0,sub,sp,notes,'hyd'+i);
      sub.forEach(e=>{e.count*=k;e.envLabel=uniF(m[2])+' (water of hydration)'});sp.forEach(x=>x.mult*=k);entries.push(...sub);parts.push(...sp)}
  });
  // sanity check
  let tot=new Q(0);entries.forEach(e=>tot=tot.add(e.value.mul(e.count)));
  if(!tot.eq(charge))throw chemErr('These oxidation numbers do not add up to the charge — please check the formula.');
  // slots
  const byEl={};entries.forEach(e=>{(byEl[e.el]=byEl[e.el]||[]).push(e)});
  const order=[];entries.forEach(e=>{if(!order.includes(e.el))order.push(e.el)});
  const slots=[],envSlot={};
  for(const el of order){const list=byEl[el];const same=list.every(e=>e.value.eq(list[0].value));
    if(same){const whys=[...new Set(list.map(e=>e.why))];const s={key:el,el,label:el,count:list.reduce((a,e)=>a+e.count,0),value:list[0].value,rule:list[0].rule,why:whys.join(' ')};slots.push(s);list.forEach(e=>envSlot[e.env+'|'+el]=el)}
    else{const envs={};list.forEach(e=>{(envs[e.env]=envs[e.env]||[]).push(e)});
      for(const env in envs){const L=envs[env];const key=el+'@'+env;slots.push({key,el,label:`${el} in ${L[0].envLabel}`,count:L.reduce((a,e)=>a+e.count,0),value:L[0].value,rule:L[0].rule,why:L[0].why});envSlot[env+'|'+el]=key}}}
  slots.forEach(s=>{s.fraction=!s.value.isInt()});
  return {raw,body,charge,html:htmlF(body.replace(/[·•*.]/g,'·'),charge),text:uniF(body.replace(/[·•*.]/g,'·'))+chargeSup(charge),slots,entries,parts,envSlot};
}

// ---------- Parse a student's answer ----------
function parseAnswer(str){
  let s=String(str).trim().replace(/[−–—]/g,'-').replace(/\s+/g,'');if(!s)return null;
  let m;
  if((m=s.match(/^([+-]?)(\d+)\/(\d+)$/))){const v=new Q(+m[2],+m[3]);return m[1]==='-'?v.neg():v}
  if((m=s.match(/^(\d+)([+-])$/)))return Q.of((m[2]==='-'?-1:1)*+m[1]);   // "3+" style
  if((m=s.match(/^([+-]?)(\d*\.\d+|\d+)$/))){const x=parseFloat(m[2])*(m[1]==='-'?-1:1);
    if(Number.isInteger(x))return Q.of(x);
    for(let d=2;d<=12;d++){const n=Math.round(x*d);if(Math.abs(n/d-x)<0.011)return new Q(n,d)}
    return {approx:x};}
  return undefined;
}
function answerMatches(ans,val){if(!ans)return false;if(ans instanceof Q)return ans.eq(val);if(ans.approx!==undefined)return Math.abs(ans.approx-val.val())<0.011;return false}

// ===================== 3D structure building =====================
function parseSmi(s){
  const atoms=[],bonds=[];let i=0,prev=-1,order=1;const stack=[];
  while(i<s.length){const c=s[i];
    if(c==='('){stack.push(prev);i++}else if(c===')'){prev=stack.pop();i++}
    else if(c==='='){order=2;i++}else if(c==='#'){order=3;i++}else if(c==='-'&&s[i+1]!==']'){order=1;i++}
    else{let el,fc=0;
      if(c==='['){const j=s.indexOf(']',i);const t=s.slice(i+1,j);const mm=t.match(/^([A-Z][a-z]?)(.*)$/);el=mm[1];const ch=mm[2];
        if(ch){const a=ch.match(/^([+-])(\d*)$/);if(a)fc=(a[1]==='-'?-1:1)*(a[2]?+a[2]:1);else{const b=ch.match(/^(\d+)([+-])$/);fc=(b[2]==='-'?-1:1)*+b[1]}}
        i=j+1}
      else{el=readSym(s,i);i+=el.length}
      atoms.push({el,fc});const idx=atoms.length-1;if(prev>=0)bonds.push([prev,idx,order]);prev=idx;order=1;}
  }
  return {atoms,bonds};
}
const V3={add:(a,b)=>[a[0]+b[0],a[1]+b[1],a[2]+b[2]],sub:(a,b)=>[a[0]-b[0],a[1]-b[1],a[2]-b[2]],mul:(a,k)=>[a[0]*k,a[1]*k,a[2]*k],dot:(a,b)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2],
  cross:(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]],len:a=>Math.hypot(a[0],a[1],a[2]),norm:a=>{const l=Math.hypot(a[0],a[1],a[2])||1;return[a[0]/l,a[1]/l,a[2]/l]}};
function rotAxis(ax,ang){const[x,y,z]=V3.norm(ax),c=Math.cos(ang),s=Math.sin(ang),t=1-c;return[[t*x*x+c,t*x*y-s*z,t*x*z+s*y],[t*x*y+s*z,t*y*y+c,t*y*z-s*x],[t*x*z-s*y,t*y*z+s*x,t*z*z+c]]}
const mv=(M,v)=>[M[0][0]*v[0]+M[0][1]*v[1]+M[0][2]*v[2],M[1][0]*v[0]+M[1][1]*v[1]+M[1][2]*v[2],M[2][0]*v[0]+M[2][1]*v[1]+M[2][2]*v[2]];
const mm=(A,B)=>A.map((r,i)=>[0,1,2].map(j=>A[i][0]*B[0][j]+A[i][1]*B[1][j]+A[i][2]*B[2][j]));
const I3=[[1,0,0],[0,1,0],[0,0,1]];
function alignM(a,b){a=V3.norm(a);b=V3.norm(b);const c=V3.dot(a,b);if(c>0.99999)return I3;if(c<-0.99999){let p=V3.cross(a,[1,0,0]);if(V3.len(p)<0.1)p=V3.cross(a,[0,1,0]);return rotAxis(p,Math.PI)}return rotAxis(V3.cross(a,b),Math.acos(c))}
const T=1/Math.sqrt(3);
const IDEAL={1:[[1,0,0]],2:[[1,0,0],[-1,0,0]],3:[[1,0,0],[-.5,.8660254,0],[-.5,-.8660254,0]],4:[[T,T,T],[T,-T,-T],[-T,T,-T],[-T,-T,T]],
  5:[[0,0,1],[0,0,-1],[1,0,0],[-.5,.8660254,0],[-.5,-.8660254,0]],6:[[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]]};
function bondedSlots(sn,lp){const n=sn-lp;
  if(sn===5){const lpSlots=[2,3,4].slice(0,lp);return[0,1,2,3,4].filter(i=>!lpSlots.includes(i)).slice(0,n)}
  if(sn===6){const lpSlots=[4,5].slice(0,lp);return[0,1,2,3,4,5].filter(i=>!lpSlots.includes(i)).slice(0,n)}
  return[...Array(sn).keys()].slice(0,n)}
function bondLen(a,b,o){return(EL[a].rc+EL[b].rc)*(o===2?0.9:o===3?0.85:1)}
// embed a tree graph. opts: {lp:{idx:n}, rootPos, rootParentDir, ref}
function embed(g,opts={}){
  const n=g.atoms.length,adj=[...Array(n)].map(()=>[]);g.bonds.forEach(([a,b,o])=>{adj[a].push([b,o]);adj[b].push([a,o])});
  const pos=Array(n),lpOf=i=>{if(opts.lp&&opts.lp[i]!==undefined)return opts.lp[i];const a=g.atoms[i];if(isMetal(a.el))return 0;
    const bo=adj[i].reduce((s,[,o])=>s+o,0);return Math.max(0,Math.round((valenceE(a.el)-bo-a.fc)/2))};
  function place(i,parent,pdir,ref){
    const kids=adj[i].filter(([j])=>j!==parent);const hasP=!!pdir;
    if(!kids.length)return;
    const deg=kids.length+(hasP?1:0);let lp=lpOf(i);if(deg+lp>6)lp=6-deg;const sn=Math.max(deg+lp,deg);
    const dirs=IDEAL[Math.min(sn,6)]||IDEAL[6];const bs=bondedSlots(Math.min(sn,6),lp);let R=I3;
    if(hasP){R=alignM(dirs[bs[0]],pdir);
      if(bs.length>1&&ref){const p=V3.norm(pdir);const v=mv(R,dirs[bs[1]]);const want=V3.mul(ref,-1);
        const pv=V3.sub(v,V3.mul(p,V3.dot(v,p))),pw=V3.sub(want,V3.mul(p,V3.dot(want,p)));
        if(V3.len(pv)>1e-3&&V3.len(pw)>1e-3){const a=V3.norm(pv),b=V3.norm(pw);let ang=Math.acos(Math.max(-1,Math.min(1,V3.dot(a,b))));if(V3.dot(V3.cross(a,b),p)<0)ang=-ang;R=mm(rotAxis(p,ang),R)}}}
    const free=bs.slice(hasP?1:0);
    kids.forEach(([j,o],k)=>{const d=V3.norm(mv(R,dirs[free[k%free.length]]));pos[j]=V3.add(pos[i],V3.mul(d,bondLen(g.atoms[i].el,g.atoms[j].el,o)));place(j,i,V3.mul(d,-1),pdir||null)});
  }
  pos[0]=opts.rootPos||[0,0,0];place(0,-1,opts.rootParentDir||null,opts.ref||null);
  return {atoms:g.atoms.map((a,i)=>({el:a.el,p:pos[i]||[0,0,0]})),bonds:g.bonds.map(([a,b,o])=>[a,b,o])};
}

const LIB_SRC=[
 ['H2O',0,'HOH'],['H2O2',0,'HOOH'],['NH3',0,'N(H)(H)H'],['NH4',1,'[N+](H)(H)(H)H'],['H3O',1,'[O+](H)(H)H'],['CH4',0,'C(H)(H)(H)H'],
 ['CO2',0,'O=C=O'],['CO',0,'[C-]#[O+]'],['SO2',0,'O=S=O'],['SO3',0,'O=S(=O)=O'],['NO',0,'N=O'],['NO2',0,'O=[N+][O-]'],['N2O',0,'N#[N+][O-]'],
 ['N2O4',0,'O=[N+]([O-])[N+](=O)[O-]'],['H2S',0,'HSH'],['PH3',0,'P(H)(H)H'],['HCl',0,'HCl'],['HF',0,'HF'],['HCN',0,'HC#N'],
 ['CCl4',0,'C(Cl)(Cl)(Cl)Cl'],['CS2',0,'S=C=S'],['PCl3',0,'P(Cl)(Cl)Cl'],['PCl5',0,'P(Cl)(Cl)(Cl)(Cl)Cl'],['SF6',0,'S(F)(F)(F)(F)(F)F'],
 ['OF2',0,'FOF'],['ClF3',0,'Cl(F)(F)F'],['BrF3',0,'Br(F)(F)F'],['XeF4',0,'[Xe](F)(F)(F)F'],['IF5',0,'I(F)(F)(F)(F)F'],['ICl',0,'ICl'],
 ['Cl2O7',0,'O(Cl(=O)(=O)=O)Cl(=O)(=O)=O'],['ClO2',0,'O=Cl=O'],['Cl2O',0,'ClOCl'],['HClO',0,'HOCl'],['HOF',0,'HOF'],['N2H4',0,'N(H)(H)N(H)H'],
 ['O2',0,'O=O'],['O3',0,'O=[O+][O-]'],['N2',0,'N#N'],['HNO3',0,'HO[N+](=O)[O-]'],['H2SO4',0,'HOS(=O)(=O)OH'],['H3PO4',0,'HOP(=O)(OH)OH'],['HClO4',0,'HOCl(=O)(=O)=O'],
 ['OH',-1,'[O-]H'],['NO3',-1,'[O-][N+](=O)[O-]'],['NO2',-1,'[O-]N=O'],['SO4',-2,'[O-]S(=O)(=O)[O-]'],['SO3',-2,'[O-]S(=O)[O-]'],['HSO4',-1,'HOS(=O)(=O)[O-]'],
 ['PO4',-3,'[O-]P(=O)([O-])[O-]'],['HPO4',-2,'HOP(=O)([O-])[O-]'],['H2PO4',-1,'HOP(=O)(OH)[O-]'],['CO3',-2,'[O-]C(=O)[O-]'],['HCO3',-1,'HOC(=O)[O-]'],
 ['ClO4',-1,'[O-]Cl(=O)(=O)=O'],['ClO3',-1,'[O-]Cl(=O)=O'],['ClO2',-1,'[O-]Cl=O'],['ClO',-1,'[O-]Cl'],['BrO3',-1,'[O-]Br(=O)=O'],['IO3',-1,'[O-]I(=O)=O'],['IO4',-1,'[O-]I(=O)(=O)=O'],
 ['MnO4',-1,'[O-][Mn](=O)(=O)=O'],['CrO4',-2,'[O-][Cr](=O)(=O)[O-]'],['Cr2O7',-2,'[O-][Cr](=O)(=O)O[Cr](=O)(=O)[O-]'],['C2O4',-2,'[O-]C(=O)C(=O)[O-]'],
 ['CN',-1,'[C-]#N'],['SCN',-1,'[S-]C#N'],['S2O3',-2,'[O-]S(=O)(=S)[O-]'],['S4O6',-2,'[O-]S(=O)(=O)SSS(=O)(=O)[O-]'],['N3',-1,'[N-]=[N+]=[N-]'],
 ['O2',-2,'[O-][O-]'],['O2',-1,'[O-]O'],['C2',-2,'[C-]#[C-]'],['AsO4',-3,'[O-][As](=O)([O-])[O-]'],['BH4',-1,'[B-](H)(H)(H)H'],['AlH4',-1,'[Al-](H)(H)(H)H']];
const LIB={};LIB_SRC.forEach(([f,q,s])=>{LIB[keyOf(compOf(f))+'|'+q]=s});
function specialStruct(key){
  if(key==='S8|0'){const atoms=[],bonds=[];for(let k=0;k<8;k++){const a=k*Math.PI/4;atoms.push({el:'S',p:[2.37*Math.cos(a),2.37*Math.sin(a),k%2?0.48:-0.48]});bonds.push([k,(k+1)%8,1])}return{atoms,bonds}}
  if(key==='P4|0'){const s=2.21/(2*Math.SQRT2);const P=[[s,s,s],[s,-s,-s],[-s,s,-s],[-s,-s,s]];const bonds=[];for(let i=0;i<4;i++)for(let j=i+1;j<4;j++)bonds.push([i,j,1]);return{atoms:P.map(p=>({el:'P',p})),bonds}}
  return null;
}
function autoGraph(c,charge){
  const els=Object.keys(c),n=atomCount(c);
  if(n===2){const atoms=[];els.forEach(e=>{for(let k=0;k<c[e];k++)atoms.push({el:e,fc:0})});return{g:{atoms,bonds:[[0,1,1]]},lp:{}}}
  const cands=els.filter(e=>c[e]===1&&e!=='H'&&e!=='F'&&!isMetal(e)).sort((a,b)=>EL[a].en-EL[b].en);
  if(!cands.length)return null;const cen=cands[0];
  const atoms=[{el:cen,fc:0}],bonds=[];const Os=[];
  for(const e of els){if(e===cen||e==='H')continue;for(let k=0;k<c[e];k++){atoms.push({el:e,fc:0});bonds.push([0,atoms.length-1,1]);if(e==='O')Os.push(atoms.length-1)}}
  const nH=c.H||0;let hOnO=0;
  for(let k=0;k<nH;k++){if(cen!=='O'&&Os.length>k&&EL[cen].group>=15&&cen!=='N'||(cen==='N'&&Os.length>=3&&k<Os.length)){atoms.push({el:'H',fc:0});bonds.push([Os[k],atoms.length-1,1]);hOnO++}
    else{atoms.push({el:'H',fc:0});bonds.push([0,atoms.length-1,1])}}
  const ve=els.reduce((s,e)=>s+valenceE(e)*c[e],0)-charge;
  const deg=atoms.map(()=>0);bonds.forEach(([a,b])=>{deg[a]++;deg[b]++});
  let used=2*bonds.length;atoms.forEach((a,i)=>{if(i===0||a.el==='H')return;used+=Math.max(0,8-2*deg[i])});
  const lpC=Math.max(0,Math.round((ve-used)/2));
  return{g:{atoms,bonds},lp:{0:lpC}};
}
function hillFormula(c){const els=Object.keys(c);let order;if(c.C){order=['C',...(c.H?['H']:[]),...els.filter(e=>e!=='C'&&e!=='H').sort()]}else order=els.sort();return order.map(e=>e+(c[e]>1?c[e]:'')).join('')}
function parseSDF(t){const L=t.split(/\r?\n/);const na=+L[3].slice(0,3),nb=+L[3].slice(3,6);const atoms=[],bonds=[];
  for(let i=0;i<na;i++){const l=L[4+i];atoms.push({el:l.slice(31,34).trim(),p:[+l.slice(0,10),+l.slice(10,20),+l.slice(20,30)]})}
  for(let i=0;i<nb;i++){const l=L[4+na+i];bonds.push([+l.slice(0,3)-1,+l.slice(3,6)-1,+l.slice(6,9)])}return{atoms,bonds}}
async function pubchem3D(c){
  const f=hillFormula(c);const ctl=typeof AbortController!=='undefined'?new AbortController():null;const tm=setTimeout(()=>ctl&&ctl.abort(),7000);
  try{const r=await fetch(`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastformula/${f}/cids/JSON?MaxRecords=1`,{signal:ctl&&ctl.signal});if(!r.ok)return null;
    const j=await r.json();const cid=j.IdentifierList&&j.IdentifierList.CID&&j.IdentifierList.CID[0];if(!cid)return null;
    const r2=await fetch(`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/record/SDF?record_type=3d`,{signal:ctl&&ctl.signal});if(!r2.ok)return null;
    const s=parseSDF(await r2.text());const cc={};s.atoms.forEach(a=>cc[a.el]=(cc[a.el]||0)+1);if(keyOf(cc)!==keyOf(c))return null;s.source='PubChem CID '+cid;return s;}
  catch(e){return null}finally{clearTimeout(tm)}
}
// returns {atoms:[{el,p,on,env}],bonds,source} for one unit (without multiplicity)
async function buildUnit(part,allowNet){
  const c=part.comp,q=Q.of(part.charge).val(),key=keyOf(c)+'|'+q;let s=null,src='built-in';
  if(atomCount(c)===1)s={atoms:[{el:Object.keys(c)[0],p:[0,0,0]}],bonds:[]};
  else if(specialStruct(key))s=specialStruct(key);
  else if(LIB[key])s=embed(parseSmi(LIB[key]));
  else{if(allowNet&&q===0&&!Object.keys(c).some(isMetal)){s=await pubchem3D(c);if(s)src=s.source}
    if(!s){const ag=autoGraph(c,q);if(ag){s=embed(ag.g,{lp:ag.lp});src='VSEPR model'}}}
  if(!s)return null;
  s.atoms.forEach(a=>{a.on=part.sol[a.el];a.env=part.env});s.source=src;return s;
}
function buildComplex(part){
  const M=part.metal;const atoms=[{el:M,p:[0,0,0],on:part.metalValue,env:part.env}],bonds=[];
  const inst=[];part.ligands.forEach(l=>{for(let k=0;k<l.count;k++)inst.push(l)});
  const cn=inst.reduce((s,l)=>s+(l.L.type==='mono'?1:2),0);
  let dirs;
  if(cn===4){const sq=['Pt','Pd','Au'].includes(M)||(M==='Ni'&&inst.every(l=>l.L.f==='CN'))||M==='Cu';dirs=sq?[[1,0,0],[0,1,0],[-1,0,0],[0,-1,0]]:IDEAL[4]}
  else dirs=(IDEAL[cn]||IDEAL[6]).map(d=>d.slice());
  const used=new Set();const pairs=[[0,2],[1,4],[3,5]];
  const chel=inst.filter(l=>l.L.type!=='mono'),mono=inst.filter(l=>l.L.type==='mono');
  const add=(el,p,on,env)=>{atoms.push({el,p,on,env});return atoms.length-1};
  chel.forEach((l,k)=>{const[i1,i2]=pairs[k]||[0,2];used.add(i1);used.add(i2);const d1=dirs[i1],d2=dirs[i2];
    const D=l.L.type==='ox'?'O':'N',L0=EL[M].rc+EL[D].rc;const A=V3.mul(d1,L0),B=V3.mul(d2,L0);const m=V3.norm(V3.add(d1,d2)),u=V3.norm(V3.sub(A,B)),w=V3.norm(V3.cross(m,u));
    const on=e=>l.sol[e];const a=add(D,A,on(D),l.env),b=add(D,B,on(D),l.env);bonds.push([0,a,1],[0,b,1]);
    const C1=V3.add(V3.add(A,V3.mul(m,1.1)),V3.mul(u,-0.66)),C2=V3.add(V3.add(B,V3.mul(m,1.1)),V3.mul(u,0.66));
    const c1=add('C',C1,on('C'),l.env),c2=add('C',C2,on('C'),l.env);bonds.push([a,c1,1],[c1,c2,1],[c2,b,1]);
    if(l.L.type==='ox'){const o1=add('O',V3.add(C1,V3.mul(V3.norm(V3.add(V3.mul(m,.6),u)),1.22)),on('O'),l.env),o2=add('O',V3.add(C2,V3.mul(V3.norm(V3.sub(V3.mul(m,.6),u)),1.22)),on('O'),l.env);bonds.push([c1,o1,2],[c2,o2,2])}
    else{for(const[X,sg]of[[a,1],[b,-1]]){for(const s of[1,-1]){const h=add('H',V3.add(atoms[X].p,V3.mul(V3.norm(V3.add(V3.add(V3.mul(m,-.25),V3.mul(u,.55*sg)),V3.mul(w,.8*s))),1.01)),on('H'),l.env);bonds.push([X,h,1])}}
      for(const[X,sg]of[[c1,1],[c2,-1]]){for(const s of[1,-1]){const h=add('H',V3.add(atoms[X].p,V3.mul(V3.norm(V3.add(V3.add(V3.mul(m,.55),V3.mul(u,.3*sg)),V3.mul(w,.8*s))),1.09)),on('H'),l.env);bonds.push([X,h,1])}}}
  });
  const freeSlots=dirs.map((_,i)=>i).filter(i=>!used.has(i));
  mono.forEach((l,k)=>{const d=dirs[freeSlots[k%freeSlots.length]];const g=parseSmi(l.L.smi);const D=g.atoms[0].el,L0=EL[M].rc+EL[D].rc;
    const lp={0:Math.max(0,Math.round((valenceE(D)-g.bonds.filter(b=>b[0]===0||b[1]===0).reduce((s,b)=>s+b[2],0)-g.atoms[0].fc)/2)-1)};
    let ref=V3.cross(d,[0,0,1]);if(V3.len(ref)<0.1)ref=V3.cross(d,[0,1,0]);
    const s=embed(g,{lp,rootPos:V3.mul(d,L0),rootParentDir:V3.mul(d,-1),ref:V3.norm(ref)});const off=atoms.length;
    s.atoms.forEach(a=>atoms.push({el:a.el,p:a.p,on:l.sol[a.el],env:l.env}));s.bonds.forEach(([x,y,o])=>bonds.push([x+off,y+off,o]));bonds.push([0,off,1]);});
  return{atoms,bonds,source:'built-in'};
}
// Build the full scene: all ions/molecules side by side
async function buildScene(an,allowNet){
  const units=[];let missing=false;const sources=new Set();
  if(an.isOrg)return orgScene(an);
  for(const part of an.parts){let s;if(part.kind==='complex')s=buildComplex(part);else s=await buildUnit(part,allowNet);
    if(!s){missing=true;continue}sources.add(s.source);const n=Math.min(part.mult,8);for(let k=0;k<n;k++)units.push(s)}
  if(!units.length)return null;
  const MAX=14;const use=units.slice(0,MAX);
  // center each, radius
  const placed=use.map(s=>{const cen=s.atoms.reduce((a,x)=>V3.add(a,x.p),[0,0,0]).map(v=>v/s.atoms.length);const r=Math.max(0.8,...s.atoms.map(a=>V3.len(V3.sub(a.p,cen))+0.6));return{s,cen,r}});
  const perRow=placed.length<=5?placed.length:Math.ceil(Math.sqrt(placed.length*1.6));
  const atoms=[],bonds=[];let x=0,y=0,rowH=0,col=0;const rows=[[]];
  placed.forEach(u=>{if(col===perRow){rows.push([]);col=0}rows[rows.length-1].push(u);col++});
  let yCur=0;rows.forEach(row=>{const h=Math.max(...row.map(u=>u.r));const width=row.reduce((a,u)=>a+2*u.r,0)+(row.length-1)*1.2;let xc=-width/2;
    row.forEach(u=>{const off=[xc+u.r-u.cen[0],yCur-u.cen[1],-u.cen[2]];const base=atoms.length;u.s.atoms.forEach(a=>atoms.push({...a,p:V3.add(a.p,off)}));u.s.bonds.forEach(([a,b,o])=>bonds.push([a+base,b+base,o]));xc+=2*u.r+1.2});
    yCur-=2*h+1.2});
  const cy=atoms.reduce((a,x)=>a+x.p[1],0)/atoms.length;atoms.forEach(a=>a.p[1]-=cy);
  return{atoms,bonds,missing,truncated:units.length>MAX,sources:[...sources]};
}

// ===================== Problem database =====================
const LEVELS={
 1:{name:'Binary ionic compounds',desc:'Metal + nonmetal, plus free elements and simple ions'},
 2:{name:'Binary covalent compounds',desc:'Two nonmetals'},
 3:{name:'Polyatomic ions',desc:'Polyatomic ions, acids and their salts'},
 4:{name:'Exceptions & fractional',desc:'Peroxides, superoxides, hydrides, O–F and halogen–halogen compounds, average oxidation numbers'},
 5:{name:'Coordination compounds',desc:'Complex ions with ligands and counter-ions'}};
const PROBLEMS={
1:[['NaCl','sodium chloride'],['MgO','magnesium oxide'],['Al2O3','aluminum oxide'],['CaF2','calcium fluoride'],['K2S','potassium sulfide'],['Li3N','lithium nitride'],['FeCl3','iron(III) chloride'],['Fe2O3','iron(III) oxide'],['FeO','iron(II) oxide'],['Cu2O','copper(I) oxide'],['CuCl2','copper(II) chloride'],['SnO2','tin(IV) oxide'],['MnO2','manganese(IV) oxide'],['ZnBr2','zinc bromide'],['Cr2O3','chromium(III) oxide'],['PbI2','lead(II) iodide'],['BaCl2','barium chloride'],['O2','oxygen gas'],['S8','sulfur (S₈)'],['P4','white phosphorus'],['Fe','iron metal'],['Fe^3+','iron(III) ion'],['S^2-','sulfide ion']],
2:[['CO2','carbon dioxide'],['CO','carbon monoxide'],['SO2','sulfur dioxide'],['SO3','sulfur trioxide'],['NO','nitrogen monoxide'],['NO2','nitrogen dioxide'],['N2O','dinitrogen monoxide'],['N2O4','dinitrogen tetroxide'],['NH3','ammonia'],['H2S','hydrogen sulfide'],['H2O','water'],['HCl','hydrogen chloride'],['PCl3','phosphorus trichloride'],['PCl5','phosphorus pentachloride'],['CCl4','carbon tetrachloride'],['CS2','carbon disulfide'],['SF6','sulfur hexafluoride'],['Cl2O7','dichlorine heptoxide'],['ClO2','chlorine dioxide'],['N2H4','hydrazine'],['IF5','iodine pentafluoride'],['XeF4','xenon tetrafluoride'],['CH4','methane'],['PH3','phosphine']],
3:[['SO4^2-','sulfate ion'],['NO3^-','nitrate ion'],['PO4^3-','phosphate ion'],['MnO4^-','permanganate ion'],['Cr2O7^2-','dichromate ion'],['NH4^+','ammonium ion'],['ClO4^-','perchlorate ion'],['CO3^2-','carbonate ion'],['NO2^-','nitrite ion'],['K2Cr2O7','potassium dichromate'],['KMnO4','potassium permanganate'],['Na2SO4','sodium sulfate'],['Fe2(SO4)3','iron(III) sulfate'],['NH4NO3','ammonium nitrate'],['(NH4)2SO4','ammonium sulfate'],['Ca3(PO4)2','calcium phosphate'],['NaHCO3','sodium hydrogen carbonate'],['H2SO4','sulfuric acid'],['HNO3','nitric acid'],['H3PO4','phosphoric acid'],['NaClO3','sodium chlorate'],['KNO2','potassium nitrite'],['CuSO4','copper(II) sulfate'],['Na2CrO4','sodium chromate'],['NaClO','sodium hypochlorite'],['NaH2PO4','sodium dihydrogen phosphate'],['Na2S2O3','sodium thiosulfate'],['HClO4','perchloric acid'],['Mg(OH)2','magnesium hydroxide'],['Cu(NO3)2','copper(II) nitrate'],['CuSO4·5H2O','copper(II) sulfate pentahydrate']],
4:[['H2O2','hydrogen peroxide'],['Na2O2','sodium peroxide'],['BaO2','barium peroxide'],['KO2','potassium superoxide'],['O2^-','superoxide ion'],['OF2','oxygen difluoride'],['HOF','hypofluorous acid'],['NaH','sodium hydride'],['CaH2','calcium hydride'],['LiAlH4','lithium aluminum hydride'],['LiH','lithium hydride'],['ICl','iodine monochloride'],['BrF3','bromine trifluoride'],['ClF3','chlorine trifluoride'],['HClO','hypochlorous acid'],['Cl2O','dichlorine monoxide'],['Fe3O4','iron(II,III) oxide (magnetite)'],['Pb3O4','lead(II,IV) oxide'],['Mn3O4','manganese(II,III) oxide'],['Na2S4O6','sodium tetrathionate'],['S4O6^2-','tetrathionate ion'],['NaN3','sodium azide']],
5:[['[Co(NH3)6]Cl3','hexaamminecobalt(III) chloride'],['[Co(NH3)5Cl]Cl2','pentaamminechloridocobalt(III) chloride'],['[Co(NH3)4Cl2]Cl','tetraamminedichloridocobalt(III) chloride'],['K4[Fe(CN)6]','potassium hexacyanoferrate(II)'],['K3[Fe(CN)6]','potassium hexacyanoferrate(III)'],['[Cu(NH3)4]SO4','tetraamminecopper(II) sulfate'],['[Pt(NH3)2Cl2]','cisplatin'],['K2[PtCl4]','potassium tetrachloroplatinate(II)'],['K2[PtCl6]','potassium hexachloroplatinate(IV)'],['[Ni(CO)4]','tetracarbonylnickel'],['Na2[Zn(OH)4]','sodium tetrahydroxozincate'],['[Cr(H2O)6]Cl3','hexaaquachromium(III) chloride'],['K3[Fe(C2O4)3]','potassium trioxalatoferrate(III)'],['[Co(en)3]Cl3','tris(ethylenediamine)cobalt(III) chloride'],['K2[Ni(CN)4]','potassium tetracyanonickelate(II)'],['[Ag(NH3)2]Cl','diamminesilver(I) chloride'],['K[Ag(CN)2]','potassium dicyanoargentate(I)'],['Na3[Co(NO2)6]','sodium hexanitrocobaltate(III)'],['[Fe(H2O)6]^3+','hexaaquairon(III) ion'],['[Fe(SCN)(H2O)5]^2+','pentaaquathiocyanatoiron(III) ion']]
};
function findName(an){for(const lv in PROBLEMS)for(const[f,n]of PROBLEMS[lv]){try{const a=normalizeInput(f);if(a.body===an.body&&a.charge===an.charge)return n}catch(e){}}return ''}

if(typeof module!=='undefined')module.exports={Q,analyze,parseAnswer,answerMatches,buildScene,PROBLEMS,LEVELS,RULES,RULE_NUM,hintFor};
