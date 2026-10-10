const GREEK = { alpha: 'α', beta: 'β', gamma: 'γ', delta: 'δ', epsilon: 'ϵ', varepsilon: 'ε', zeta: 'ζ', eta: 'η', theta: 'θ', vartheta: 'ϑ', iota: 'ι', kappa: 'κ', lambda: 'λ', mu: 'μ', nu: 'ν', xi: 'ξ', pi: 'π', varpi: 'ϖ', rho: 'ρ', varrho: 'ϱ', sigma: 'σ', varsigma: 'ς', tau: 'τ', upsilon: 'υ', phi: 'ϕ', varphi: 'φ', chi: 'χ', psi: 'ψ', omega: 'ω', Gamma: 'Γ', Delta: 'Δ', Theta: 'Θ', Lambda: 'Λ', Xi: 'Ξ', Pi: 'Π', Sigma: 'Σ', Upsilon: 'Υ', Phi: 'Φ', Psi: 'Ψ', Omega: 'Ω' };
const SYMBOLS = { partial: '∂', nabla: '∇', infty: '∞', ell: 'ℓ' };
const OPS = { cdot: '⋅', times: '×', pm: '±', mp: '∓', div: '÷', circ: '∘', approx: '≈', equiv: '≡', sim: '∼', propto: '∝', le: '≤', leq: '≤', ge: '≥', geq: '≥', ne: '≠', neq: '≠', ll: '≪', gg: '≫', lt: '&lt;', gt: '&gt;', in: '∈', to: '→', rightarrow: '→', leftarrow: '←', gets: '←', Rightarrow: '⇒', ldots: '…', cdots: '⋯', prime: '′', '{': '{', '}': '}', '|': '‖', '%': '%', '&': '&amp;' };
const PLAIN = { '+': '+', '-': '−', '=': '=', '<': '&lt;', '>': '&gt;', '(': '(', ')': ')', '[': '[', ']': ']', ',': ',', '.': '.', '/': '/', '|': '|', '!': '!', ':': ':', ';': ';', '*': '∗' };
const SPACES = { ',': '0.1667em', ':': '0.2222em', '>': '0.2222em', ';': '0.2778em', '!': '-0.1667em', ' ': '0.3333em', quad: '1em', qquad: '2em' };
const BIG = { sum: '∑', prod: '∏', int: '∫', iint: '∬', oint: '∮' }, LIMITS = new Set(['sum', 'prod', 'max', 'min', 'lim']);
const FUNCTIONS = new Set(['ln', 'log', 'exp', 'sin', 'cos', 'tan', 'sinh', 'cosh', 'tanh', 'max', 'min', 'lim']);
const ACCENTS = { dot: '˙', ddot: '¨', hat: 'ˆ', bar: '¯', tilde: '˜', vec: '→', overline: '‾', widehat: 'ˆ', widetilde: '˜' }, WIDE = new Set(['overline', 'widehat', 'widetilde']);
const DELIMITERS = { '(': '(', ')': ')', '[': '[', ']': ']', '|': '|', '.': '', '\\{': '{', '\\}': '}', '\\|': '‖', '\\langle': '⟨', '\\rangle': '⟩' };
const FONTS = { mathbf: 'bf', boldsymbol: 'bm', mathrm: 'rm', operatorname: 'rm' };
const FIXED = new Set('()[]{}|‖/'), TIGHT = new Set('|‖/.!′″‴⁗'), BINARY = new Set('+−±∓×⋅÷∘∗');
const LEADING = new Set(['=', '&lt;', '&gt;', ':', '≈', '≡', '∼', '∝', '≤', '≥', '≠', '≪', '≫', '∈', '→', '←', '⇒', '(', '[', '{', ',', ';']);
const LATIN = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz', GREEK_ORDER = 'ΑΒΓΔΕΖΗΘΙΚΛΜΝΞΟΠΡϴΣΤΥΦΧΨΩ∇αβγδεζηθικλμνξοπρςστυφχψω∂ϵϑϰϕϱϖ', DIGITS = '0123456789';
const ALPHANUMERIC = { bold: [0x1d400, 0x1d6a8], 'bold-italic': [0x1d468, 0x1d71c] };
const UPPER_GREEK = /[Α-Ω]/, COMMAND = /\\(?:[a-zA-Z]+|[^])/y, NUMBER = /\d+(?:\.\d+)?/y, WORD = /[a-zA-Z]+/y, PRIMES = /'+/y;

const has = (table, key) => Object.hasOwn(table, key);
const escape = (s) => s.replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' })[c]);
const wrap = (items) => (items.length === 1 ? items[0] : `<mrow>${items.join('')}</mrow>`);
const space = (width) => (width.startsWith('-') ? `<mrow style="margin-left: ${width}"></mrow>` : `<mspace width="${width}"></mspace>`);
const op = (c, tight = TIGHT.has(c)) => `<mo${FIXED.has(c) ? ' stretchy="false"' : ''}${tight ? ' lspace="0em" rspace="0em"' : ''}>${c}</mo>`;
const mo = (c) => ({ ml: op(c), mo: c });
const mathChar = (ch, v) => { const [latin, greek] = ALPHANUMERIC[v], a = LATIN.indexOf(ch), g = GREEK_ORDER.indexOf(ch), d = DIGITS.indexOf(ch); return a >= 0 ? String.fromCodePoint(latin + a) : g >= 0 ? String.fromCodePoint(greek + g) : d >= 0 ? String.fromCodePoint(0x1d7ce + d) : ch; };
const APPLY = '<mo>&#x2061;</mo>', THIN = space('0.1667em');

function parse(src, display) {
  let i = 0, font = null;
  const fail = (what) => { throw new SyntaxError(`${what} at position ${i}`); };
  const sticky = (re) => { re.lastIndex = i; return re.exec(src)?.[0] ?? null; };
  const skip = () => { while (i < src.length && /\s/.test(src[i])) i++; };
  const name = () => sticky(COMMAND)?.slice(1) ?? null;
  const command = () => { const n = name(); if (n === null) fail('lone backslash'); i += n.length + 1; return n; };
  const expect = (c) => { skip(); if (src[i] !== c) fail(i < src.length ? `expected ${c} but found ${src[i]}` : `missing ${c}`); i++; };
  const raw = () => { expect('{'); for (let start = i, depth = 1; i < src.length; i++) { if (src[i] === '\\') i++; else if (src[i] === '{') depth++; else if (src[i] === '}' && !--depth) return src.slice(start, i++); } fail('missing }'); };
  const ended = (bracket = false) => { skip(); const c = src[i]; return c === undefined || c === '}' || c === '&' || (bracket && c === ']') || (c === '\\' && ['\\', 'right', 'end'].includes(name())); };
  const spacing = () => { skip(); if (src[i] === '~') return 1; const n = src[i] === '\\' ? name() : null; return n !== null && (has(SPACES, n) || /^\s$/.test(n)) ? n.length + 1 : 0; };
  const spaced = () => { const at = i; for (let step; (step = spacing()); ) i += step; const c = src[i], n = c === '\\' ? name() : null; i = at; if (c === undefined) return false; if (c !== '\\') return !")]}([,;.:!=<>/&'".includes(c); return has(OPS, n) ? BINARY.has(OPS[n]) : !['\\', 'right', 'end'].includes(n); };
  const bold = () => font === 'bf' || font === 'bm';
  const variant = (ch) => (/[A-Za-z]/.test(ch) ? (font === 'bm' ? 'bold-italic' : 'bold') : /[0-9Α-Ω∇]/.test(ch) ? 'bold' : 'bold-italic');
  const styled = (ch) => { const v = variant(ch); return `<mi mathvariant="${v}">${mathChar(ch, v)}</mi>`; };
  const letter = (ch, upright = UPPER_GREEK.test(ch)) => (font === 'rm' || (!font && upright) ? `<mi mathvariant="normal">${ch}</mi>` : bold() ? styled(ch) : `<mi>${ch}</mi>`);
  const number = (n) => (bold() ? `<mn mathvariant="bold">${[...n].map((d) => mathChar(d, 'bold')).join('')}</mn>` : `<mn>${n}</mn>`);
  const text = (s) => { const t = s.replace(/\s+/g, ' ').replace(/\\([{}%&_#$ ])/g, '$1'); if (t.includes('\\')) fail('commands inside \\text are not supported'); return escape(t.replace(/~/g, ' ').replace(/^ +| +$/g, (m) => ' '.repeat(m.length))); };

  function list(bracket = false, lead = true) {
    const out = [];
    while (!ended(bracket)) {
      const base = atom(lead), unary = lead && BINARY.has(base.mo);
      if (unary) base.ml = op(base.mo, true);
      else if (!out.length && BINARY.has(base.mo)) base.ml = `<mo form="infix">${base.mo}</mo>`;
      if (!base.space) lead = !unary && (BINARY.has(base.mo) || LEADING.has(base.mo) || Boolean(base.lead || base.apply));
      const node = scripts(base);
      out.push(node.ml);
      if (node.apply) out.push(APPLY + (spaced() ? THIN : ''));
    }
    return out;
  }

  function items(lead = true) {
    skip();
    if (src[i] === '{') { i++; const out = list(false, lead); expect('}'); return out; }
    if (ended() || src[i] === '^' || src[i] === '_') fail('missing argument');
    if (/[0-9]/.test(src[i])) return [number(src[i++])];
    if (/[a-zA-Z]/.test(src[i])) return [letter(src[i++])];
    const node = atom();
    return [BINARY.has(node.mo) ? op(node.mo, true) : node.ml];
  }
  const arg = () => wrap(items());

  function scripts(base) {
    let sub = null, sup = null;
    for (skip(); src[i] === '^' || src[i] === '_'; skip()) {
      if (src[i++] === '^') { if (sup !== null) fail('double superscript'); sup = arg(); }
      else { if (sub !== null) fail('double subscript'); sub = arg(); }
    }
    if (sub === null && sup === null) return base;
    const [lower, upper, both] = base.limits && display ? ['munder', 'mover', 'munderover'] : ['msub', 'msup', 'msubsup'], tag = sub === null ? upper : sup === null ? lower : both;
    return { ml: `<${tag}>${base.ml}${sub ?? ''}${sup ?? ''}</${tag}>`, apply: base.apply };
  }

  function atom(lead = true) {
    const c = src[i];
    if (c === '{') return { ml: arg() };
    if (c === '^' || c === '_') return { ml: '<mrow></mrow>' };
    if (c === '\\') return control(command(), lead);
    if (c === "'") { const run = sticky(PRIMES); i += run.length; return { ml: op('′″‴⁗'[Math.min(run.length, 4) - 1]) }; }
    if (c === '~') { i++; return { ml: space(SPACES[' ']), space: true }; }
    if (/[0-9]/.test(c)) { const n = sticky(NUMBER); i += n.length; return { ml: number(n) }; }
    if (has(PLAIN, c)) { i++; return mo(PLAIN[c]); }
    if (font === 'rm' && /[a-zA-Z]/.test(c)) { const run = sticky(WORD); i += run.length; return { ml: run.length > 1 ? `<mi>${run}</mi>` : letter(run) }; }
    const ch = String.fromCodePoint(src.codePointAt(i));
    if (!/\p{L}/u.test(ch)) fail(`unexpected ${ch}`);
    i += ch.length;
    return { ml: letter(ch) };
  }

  function control(n, lead = true) {
    if (has(GREEK, n)) return { ml: letter(GREEK[n]) };
    if (has(SYMBOLS, n)) return { ml: letter(SYMBOLS[n], true) };
    if (has(OPS, n)) return mo(OPS[n]);
    if (has(SPACES, n) || /^\s$/.test(n)) return { ml: space(SPACES[n] ?? SPACES[' ']), space: true };
    if (has(BIG, n)) return { ml: `<mo>${BIG[n]}</mo>`, limits: LIMITS.has(n), lead: true };
    if (FUNCTIONS.has(n)) return { ml: `<mi>${n}</mi>`, limits: LIMITS.has(n), apply: true };
    if (has(ACCENTS, n)) return { ml: `<mover accent="true">${arg()}<mo stretchy="${WIDE.has(n)}">${ACCENTS[n]}</mo></mover>` };
    if (n === 'frac' || n === 'dfrac' || n === 'tfrac') { const fraction = `<mfrac>${arg()}${arg()}</mfrac>`; return { ml: n === 'frac' ? fraction : `<mstyle displaystyle="${n === 'dfrac'}" scriptlevel="0">${fraction}</mstyle>` }; }
    if (n === 'sqrt') { skip(); if (src[i] !== '[') return { ml: `<msqrt>${items().join('')}</msqrt>` }; i++; const index = wrap(list(true)); expect(']'); return { ml: `<mroot>${arg()}${index}</mroot>` }; }
    if (has(FONTS, n)) { const outer = font; font = FONTS[n]; const body = arg(); font = outer; return { ml: body, apply: n === 'operatorname' }; }
    if (n === 'text') return { ml: `<mtext>${text(raw())}</mtext>` };
    if (n === 'class') { const names = raw().trim(); if (!/^[A-Za-z_][\w-]*(?:\s+[A-Za-z_][\w-]*)*$/.test(names)) fail(`bad class name ${names}`); return { ml: `<mrow class="${names.replace(/\s+/g, ' ')}">${items(lead).join('')}</mrow>` }; }
    if (n === 'left') return fenced();
    if (n === 'begin') return aligned();
    fail(`unknown command \\${n}`);
  }

  function delimiter() {
    skip();
    const key = src[i] === '\\' ? `\\${command()}` : src[i++];
    if (!has(DELIMITERS, key ?? '')) fail(`bad delimiter ${key ?? 'at the end'}`);
    return DELIMITERS[key];
  }

  function fenced() {
    const open = delimiter(), body = list();
    if (name() !== 'right') fail('missing \\right');
    command();
    const close = delimiter(), fence = (c, form) => (c ? `<mo fence="true" form="${form}" stretchy="true">${c}</mo>` : '');
    return { ml: `<mrow>${fence(open, 'prefix')}${body.join('')}${fence(close, 'postfix')}</mrow>` };
  }

  function aligned() {
    const env = raw().trim();
    if (env !== 'aligned') fail(`unsupported environment ${env}`);
    const rows = [[]];
    for (;;) {
      const row = rows[rows.length - 1];
      row.push(list(false, row.length % 2 === 0));
      if (src[i] === '&') { i++; continue; }
      const n = name();
      if (n === '\\') { command(); rows.push([]); continue; }
      if (n !== 'end') fail('missing \\end{aligned}');
      command();
      if (raw().trim() !== env) fail('mismatched \\end');
      break;
    }
    const last = rows[rows.length - 1];
    if (rows.length > 1 && last.length === 1 && !last[0].length) rows.pop();
    const columns = Math.max(...rows.map((row) => row.length)), side = (k) => (k % 2 ? 'left' : 'right');
    const body = rows.map((row) => `<mtr>${row.map((cell, k) => `<mtd columnalign="${side(k)}">${k % 2 ? '<mi></mi>' : ''}${cell.join('')}</mtd>`).join('')}</mtr>`).join('');
    return { ml: `<mtable class="tex-aligned" displaystyle="true" columnalign="${Array.from({ length: columns }, (_, k) => side(k)).join(' ')}" columnspacing="0em">${body}</mtable>` };
  }

  const body = list();
  if (i < src.length) fail(src[i] === '}' ? 'unmatched }' : src[i] === '&' ? '& outside aligned' : `unexpected \\${name()}`);
  return body.join('');
}

export function texToMathML(tex, { display = false } = {}) {
  const source = String(tex);
  return `<math display="${display ? 'block' : 'inline'}"><semantics><mrow>${parse(source, display)}</mrow><annotation encoding="application/x-tex">${escape(source)}</annotation></semantics></math>`;
}

export function renderMath(root = document) {
  for (const el of root.querySelectorAll('.tex')) {
    if (el.hasAttribute('data-math')) continue;
    const tex = el.textContent.trim();
    try {
      const template = document.createElement('template');
      template.innerHTML = texToMathML(tex, { display: el.localName === 'div' || el.classList.contains('display') });
      el.replaceChildren(template.content);
      el.setAttribute('data-math', 'rendered');
    } catch (error) {
      el.classList.add('tex-error');
      el.setAttribute('data-math', 'error');
      console.warn(`Could not render “${tex}”: ${error.message}`);
    }
  }
}

export function alignMath(root = document) {
  const cells = [...root.querySelectorAll('.tex .tex-aligned > mtr > mtd:nth-child(odd)')];
  for (const cell of cells) cell.style.paddingLeft = '';
  const pads = cells.map((cell) => {
    if (!cell.children.length) return null;
    const box = cell.getBoundingClientRect(), style = getComputedStyle(cell);
    const right = Math.max(...[...cell.children].map((child) => child.getBoundingClientRect().right));
    const slack = box.right - parseFloat(style.paddingRight) - right;
    return slack > 0.5 ? `${parseFloat(style.paddingLeft) + slack}px` : null;
  });
  cells.forEach((cell, n) => { if (pads[n]) cell.style.paddingLeft = pads[n]; });
}

export function keepAligned(root = document) {
  let frame = 0;
  const schedule = () => { cancelAnimationFrame(frame); frame = requestAnimationFrame(() => alignMath(root)); };
  schedule();
  document.fonts?.ready.then(schedule);
  addEventListener('resize', schedule);
}
