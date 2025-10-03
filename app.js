// not pro code, just enough for this project

// tiny fraction with BigInt so we avoid float errors
const ONE = 1n;
function g(a,b){ a=a<0n?-a:a; b=b<0n?-b:b; while(b){ [a,b]=[b,a%b] } return a||1n; }
function l(a,b){ return a/g(a,b)*b; }
function norm([n,d]){ if(d===0n) throw Error("div0"); if(n===0n) return [0n,1n]; if(d<0n){ n=-n; d=-d } const x=g(n,d); return [n/x,d/x]; }
const add=(A,B)=>norm([A[0]*B[1]+B[0]*A[1], A[1]*B[1]]);
const sub=(A,B)=>norm([A[0]*B[1]-B[0]*A[1], A[1]*B[1]]);
const mul=(A,B)=>norm([A[0]*B[0], A[1]*B[1]]);
const divf=(A,B)=>norm([A[0]*B[1], A[1]*B[0]]);
const I = n => [BigInt(n), ONE];

// split on dot for hydrates, then parse each piece
function tokens(formula){
  const parts = formula.split(/[\u00B7\.]/g);
  let total = {};
  for(const p of parts){
    const got = parseOne(p.trim());
    for(const k in got) total[k] = (total[k]||0) + got[k];
  }
  return total;

  function parseOne(s){
    let i=0, lead=0;
    while(i<s.length && /\d/.test(s[i])){ lead = lead*10 + (s.charCodeAt(i)-48); i++; }
    if(lead===0) lead=1;

    const stack=[{}];
    function addK(map, el, n){ map[el] = (map[el]||0) + n; }

    while(i<s.length){
      const c = s[i];
      if(c==='('){
        stack.push({}); i++;
      }else if(c===')'){
        i++;
        let num=0; while(i<s.length && /\d/.test(s[i])){ num = num*10 + (s.charCodeAt(i)-48); i++; }
        if(num===0) num=1;
        const grp = stack.pop();
        const top = stack[stack.length-1];
        for(const el in grp) addK(top, el, grp[el]*num);
      }else if(/[A-Z]/.test(c)){
        let el=c; i++;
        if(i<s.length && /[a-z]/.test(s[i])){ el+=s[i]; i++; }
        let num=0; while(i<s.length && /\d/.test(s[i])){ num = num*10 + (s.charCodeAt(i)-48); i++; }
        if(num===0) num=1;
        addK(stack[stack.length-1], el, num);
      }else if(/\s/.test(c)){
        i++;
      }else{
        throw Error("bad char: "+c);
      }
    }
    const base = stack.pop();
    if(stack.length) throw Error("unmatched (");
    const out={}; for(const el in base) out[el]=base[el]*lead;
    return out;
  }
}

function splitEq(s){
  const m = s.split(/->|=>|→/);
  if(m.length!==2) throw Error("Use exactly one '->'");
  const L = m[0].split('+').map(x=>x.trim()).filter(Boolean);
  const R = m[1].split('+').map(x=>x.trim()).filter(Boolean);
  if(!L.length || !R.length) throw Error("Need stuff on both sides");
  return {L,R};
}

function makeMatrix(L,R){
  const all = [...L, ...R];
  const parsed = all.map(tokens);
  const elems = Array.from(new Set(parsed.flatMap(o=>Object.keys(o)))).sort();

  const A = elems.map(_ => all.map(_ => I(0)));
  elems.forEach((el, r)=>{
    parsed.forEach((obj, c)=>{
      const k = obj[el] || 0;
      const sgn = c < L.length ? 1 : -1;
      A[r][c] = I(sgn * k);
    });
  });
  return {A, all, leftCount: L.length};
}

// row-reduce to get one nullspace vector
function oneNull(A){
  const R=A.length, C=A[0].length;
  const M=A.map(row=>row.map(x=>x));
  let r=0;
  const pivCol = Array(R).fill(-1);

  for(let c=0; c<C && r<R; c++){
    let pr=r; while(pr<R && M[pr][c][0]===0n) pr++;
    if(pr===R) continue;
    if(pr!==r) [M[pr],M[r]]=[M[r],M[pr]];
    const inv = divf(I(1), M[r][c]);
    for(let j=c;j<C;j++) M[r][j]=mul(M[r][j],inv);
    for(let i=0;i<R;i++) if(i!==r && M[i][c][0]!==0n){
      const f=M[i][c];
      for(let j=c;j<C;j++) M[i][j]=sub(M[i][j], mul(f, M[r][j]));
    }
    pivCol[r]=c; r++;
  }

  const used = new Set(pivCol.filter(x=>x>=0));
  const free = [];
  for(let c=0;c<C;c++) if(!used.has(c)) free.push(c);
  if(free.length===0) throw Error("no solution (non-trivial)");

  // set first free = 1, others = 0
  const x = Array(C).fill(I(0));
  x[free[0]] = I(1);
  for(let i=0;i<R;i++){
    const pc=pivCol[i]; if(pc<0) continue;
    // row looks like e_pc + sum a_ij e_j = 0
    let sum=I(0);
    for(const fc of free){
      if(M[i][fc][0]!==0n) sum = add(sum, mul(M[i][fc], x[fc]));
    }
    x[pc] = mul(I(-1), sum);
  }

  // make integers
  let LCM=1n; for(const q of x) LCM = l(LCM, q[1]);
  const ints = x.map(q => q[0]*(LCM/q[1]));
  const neg = ints.some(v=>v<0n);
  const vec = ints.map(v=>neg?-v:v);
  let gAll = vec.reduce((a,b)=>g(a,b), 0n) || 1n;
  return vec.map(v=>v/gAll);
}

function format(all,leftCount,coef){
  const side = arr => arr.map(([k,m]) => (k===1n?"":String(k))+m).join(" + ");
  const L = all.slice(0,leftCount).map((m,i)=>[coef[i],m]);
  const R = all.slice(leftCount).map((m,i)=>[coef[i+leftCount],m]);
  return `${side(L)} -> ${side(R)}`;
}

// UI bits
const $ = s => document.querySelector(s);
const eq = $("#eq");
const out = $("#out");

function go(){
  out.textContent=""; out.className="out";
  try{
    const raw = eq.value.trim().replace(/\s+/g,' ');
    if(!raw){ out.textContent="Enter an equation."; return; }
    const {L,R} = splitEq(raw);
    const {A,all,leftCount} = makeMatrix(L,R);
    const coef = oneNull(A);
    // if most are negative, flip
    if(coef.some(c=>c<=0n)){
      const s = coef.reduce((t,c)=>t+(c<0n?-1:c>0n?1:0),0);
      if(s<0) for(let i=0;i<coef.length;i++) coef[i] = -coef[i];
    }
    out.textContent = format(all,leftCount,coef);
    out.classList.add("ok");
  }catch(e){
    out.textContent = "Error: " + e.message;
    out.classList.add("err");
  }
}

$("#bal").addEventListener("click", go);
$("#clr").addEventListener("click", ()=>{ eq.value=""; out.textContent=""; out.className="out"; });
$("#ex1").addEventListener("click", ()=> eq.value="Ca(OH)2 + H3PO4 -> Ca3(PO4)2 + H2O");
$("#ex2").addEventListener("click", ()=> eq.value="K4Fe(CN)6 + H2SO4 + H2O -> K2SO4 + FeSO4 + (NH4)2SO4 + CO");
$("#ex3").addEventListener("click", ()=> eq.value="C12H22O11 + KClO3 -> KCl + CO2 + H2O");
document.addEventListener("keydown", e=>{ if(e.key==="Enter") go(); });