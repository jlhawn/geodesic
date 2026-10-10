import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { texToMathML, renderMath } from '../explain/math.module.js';

const WRAPPER = /^<math display="(block|inline)"><semantics><mrow>([^]*)<\/mrow><annotation encoding="application\/x-tex">([^]*)<\/annotation><\/semantics><\/math>$/;
const body = (tex, options) => WRAPPER.exec(texToMathML(tex, options))[2];
const display = (tex) => body(tex, { display: true });
const THIN = '<mspace width="0.1667em"></mspace>', APPLY = '<mo>&#x2061;</mo>', BACK = '<mrow style="margin-left: -0.1667em"></mrow>';
const unary = (c) => `<mo lspace="0em" rspace="0em">${c}</mo>`;
const unescape = (s) => s.replace(/&(amp|lt|gt|quot|#x([0-9a-f]+));/gi, (_, name, hex) => (hex ? String.fromCodePoint(parseInt(hex, 16)) : { amp: '&', lt: '<', gt: '>', quot: '"' }[name]));

const ELEMENTS = new Set(['math', 'semantics', 'annotation', 'mrow', 'mi', 'mn', 'mo', 'mtext', 'mspace', 'mfrac', 'msqrt', 'mroot', 'msub', 'msup', 'msubsup', 'munder', 'mover', 'munderover', 'mtable', 'mtr', 'mtd', 'mstyle']);
const TOKENS = new Set(['mi', 'mn', 'mo', 'mtext', 'annotation']);
const ARITY = { mfrac: 2, mroot: 2, msub: 2, msup: 2, munder: 2, mover: 2, msubsup: 3, munderover: 3, semantics: 2 };
const ENTITY = /&(?!amp;|lt;|gt;|quot;|#x[0-9a-f]+;)/i;

function tree(markup) {
  const root = { name: '#root', children: [] }, stack = [root], piece = /<(\/?)([a-z]+)((?:\s[a-z-]+="[^"<>]*")*)>|([^<>]+)/gy;
  while (piece.lastIndex < markup.length) {
    const at = piece.lastIndex, m = piece.exec(markup);
    assert.ok(m, `malformed markup at ${at}: ${markup.slice(at, at + 40)}`);
    const [, closing, name, attrs, text] = m, top = stack[stack.length - 1];
    if (text !== undefined) { assert.doesNotMatch(text, ENTITY, `bare & in ${text}`); top.children.push({ text }); continue; }
    if (closing) { assert.equal(top.name, name, `</${name}> closes <${top.name}>`); stack.pop(); continue; }
    assert.doesNotMatch(attrs, ENTITY, `bare & in attributes ${attrs}`);
    const node = { name, attrs: Object.fromEntries([...attrs.matchAll(/([a-z-]+)="([^"]*)"/g)].map(([, k, v]) => [k, v])), children: [] };
    top.children.push(node);
    stack.push(node);
  }
  assert.equal(stack.length, 1, `unclosed <${stack[stack.length - 1].name}>`);
  return root.children;
}

function check(node) {
  if (node.text !== undefined) return;
  assert.ok(ELEMENTS.has(node.name), `<${node.name}> is MathML Core`);
  const elements = node.children.filter((c) => c.text === undefined);
  if (TOKENS.has(node.name)) assert.equal(elements.length, 0, `<${node.name}> holds only text`);
  if (node.name === 'mspace') assert.equal(node.children.length, 0, '<mspace> is empty');
  if (node.name in ARITY) { assert.equal(node.children.length, elements.length, `<${node.name}> has no loose text`); assert.equal(elements.length, ARITY[node.name], `<${node.name}> has ${ARITY[node.name]} children`); }
  if (node.name === 'mtable') assert.ok(elements.every((c) => c.name === 'mtr'), '<mtable> holds rows');
  if (node.name === 'mtr') assert.ok(elements.every((c) => c.name === 'mtd'), '<mtr> holds cells');
  if (['mrow', 'mfrac', 'msqrt', 'mroot', 'msub', 'msup', 'msubsup', 'munder', 'mover', 'munderover', 'mtd', 'mstyle'].includes(node.name)) assert.equal(node.children.length, elements.length, `<${node.name}> has no loose text`);
  elements.forEach(check);
}

function wellFormed(tex, options) {
  const markup = texToMathML(tex, options), nodes = tree(markup);
  assert.equal(nodes.length, 1);
  check(nodes[0]);
  const [, mode, , annotation] = WRAPPER.exec(markup);
  assert.equal(mode, options?.display ? 'block' : 'inline');
  assert.equal(unescape(annotation), tex);
  return nodes[0];
}

const corpus = [
  String.raw`\frac{\partial \pi}{\partial t} + \nabla_\sigma \cdot (\pi \mathbf{u}) + \frac{\partial (\pi \dot\sigma)}{\partial \sigma} = 0`,
  String.raw`\frac{\partial \pi}{\partial t} = -\int_0^1 \nabla_\sigma \cdot (\pi \mathbf{u}) \, \mathrm{d}\sigma`,
  String.raw`\pi \dot\sigma = -\sigma \frac{\partial \pi}{\partial t} - \int_0^\sigma \nabla_\sigma \cdot (\pi \mathbf{u}) \, \mathrm{d}\sigma'`,
  String.raw`p = \sigma \pi, \qquad \sigma = \frac{p}{\pi}, \qquad 0 \le \sigma \le 1`,
  String.raw`\frac{\partial \Phi}{\partial \sigma} = -\frac{R T}{\sigma}`,
  String.raw`\Pi = c_p \left( \frac{p}{p_0} \right)^{\kappa}, \qquad \frac{\partial \Phi}{\partial \Pi} = -\theta`,
  String.raw`\theta = T \left( \frac{p_0}{p} \right)^{R/c_p}, \qquad \kappa = \frac{R}{c_p} \approx 0.286`,
  String.raw`\frac{\partial \mathbf{u}}{\partial t} = -\class{term-a}{(\zeta + f)\, \hat{\mathbf{k}} \times \mathbf{u}} - \class{term-b}{\nabla_\sigma (\Phi + K)} - \class{term-c}{R T \nabla_\sigma \ln \pi} - \class{term-d}{\dot\sigma \frac{\partial \mathbf{u}}{\partial \sigma}} + \class{term-e}{\mathbf{F}}`,
  String.raw`(\mathbf{u} \cdot \nabla) \mathbf{u} = \zeta\, \hat{\mathbf{k}} \times \mathbf{u} + \nabla \left( \tfrac{1}{2} |\mathbf{u}|^2 \right)`,
  String.raw`\zeta = \hat{\mathbf{k}} \cdot \nabla \times \mathbf{u}, \qquad K = \tfrac{1}{2} |\mathbf{u}|^2, \qquad f = 2 \Omega \sin \varphi`,
  String.raw`\frac{\partial (\pi \theta)}{\partial t} + \nabla_\sigma \cdot (\pi \theta \mathbf{u}) + \frac{\partial (\pi \dot\sigma \theta)}{\partial \sigma} = \frac{\pi Q}{c_p} \left( \frac{p_0}{p} \right)^{\kappa}`,
  String.raw`\frac{D \theta}{D t} = \frac{\partial \theta}{\partial t} + \mathbf{u} \cdot \nabla_\sigma \theta + \dot\sigma \frac{\partial \theta}{\partial \sigma} = \frac{\theta}{c_p T} \dot{Q}`,
  String.raw`\begin{aligned} \frac{\mathrm{d} \pi_i}{\mathrm{d} t} &= -\frac{1}{A_i} \sum_{e \in \mathrm{EC}(i)} n_{e,i}\, l_e \sum_{k=1}^{K} \Delta\sigma_k\, \overline{\pi}_e\, u_{e,k} \\ (\pi \dot\sigma)_{i,k+1/2} &= (\pi \dot\sigma)_{i,k-1/2} - \Delta\sigma_k \left( \frac{1}{A_i} \sum_{e \in \mathrm{EC}(i)} n_{e,i}\, l_e\, \overline{\pi}_e\, u_{e,k} + \frac{\mathrm{d} \pi_i}{\mathrm{d} t} \right) \end{aligned}`,
  String.raw`\begin{aligned} \frac{\mathrm{d} u_{e,k}}{\mathrm{d} t} &= \class{term-a}{\sum_{e'} w_{e,e'}\, l_{e'}\, \overline{\pi}_{e'} u_{e',k}\, \frac{q_e + q_{e'}}{2}} \\ &\quad - \class{term-b}{\frac{(\Phi + K)_{c_2(e)} - (\Phi + K)_{c_1(e)}}{d_e}} - \class{term-c}{c_p \overline{\theta}_e \frac{\Pi_{c_2(e)} - \Pi_{c_1(e)}}{d_e}} \end{aligned}`,
  String.raw`K_i = \frac{1}{A_i} \sum_{e \in \mathrm{EC}(i)} \frac{l_e d_e}{4} u_e^2`,
  String.raw`(\nabla \Phi)_e = \frac{\Phi_{c_2(e)} - \Phi_{c_1(e)}}{d_e}`,
  String.raw`\Phi_k = \Phi_{k+1/2} + c_p \theta_k \left( \Pi_{k+1/2} - \Pi_k \right), \qquad \Phi_{K+1/2} = g z_s`,
  String.raw`\Pi_k = \frac{1}{1 + \kappa} \frac{\Pi_{k+1/2} \sigma_{k+1/2} - \Pi_{k-1/2} \sigma_{k-1/2}}{\sigma_{k+1/2} - \sigma_{k-1/2}}`,
  String.raw`\frac{\mathrm{d} \mathbf{y}}{\mathrm{d} t} = \mathbf{F}(\mathbf{y}), \qquad \mathbf{y} = \left( \pi_i,\ \theta_{i,k},\ u_{e,k} \right)`,
  String.raw`\mathbf{y}^{n+1} = \mathbf{y}^n + \Delta t\, \mathbf{F}(\mathbf{y}^n)`,
  String.raw`\mathbf{y}(t + \Delta t) = \mathbf{y}(t) + \Delta t\, \dot{\mathbf{y}}(t) + \frac{\Delta t^2}{2} \ddot{\mathbf{y}}(t) + O(\Delta t^3)`,
  String.raw`y^{n+1} = (1 + i \omega \Delta t)\, y^n, \qquad |1 + i \omega \Delta t| = \sqrt{1 + \omega^2 \Delta t^2} > 1`,
  String.raw`\begin{aligned} \mathbf{k}_1 &= \mathbf{F}(\mathbf{y}^n) \\ \mathbf{k}_2 &= \mathbf{F}\left( \mathbf{y}^n + \tfrac{\Delta t}{2} \mathbf{k}_1 \right) \\ \mathbf{k}_3 &= \mathbf{F}\left( \mathbf{y}^n + \tfrac{\Delta t}{2} \mathbf{k}_2 \right) \\ \mathbf{k}_4 &= \mathbf{F}\left( \mathbf{y}^n + \Delta t\, \mathbf{k}_3 \right) \\ \mathbf{y}^{n+1} &= \mathbf{y}^n + \frac{\Delta t}{6} \left( \mathbf{k}_1 + 2 \mathbf{k}_2 + 2 \mathbf{k}_3 + \mathbf{k}_4 \right) \end{aligned}`,
  String.raw`A(z) = 1 + z + \frac{z^2}{2} + \frac{z^3}{6} + \frac{z^4}{24}, \quad z = \lambda \Delta t, \quad |A(i \omega \Delta t)| \le 1 \ \text{for}\ \omega \Delta t \le 2\sqrt{2}`,
  String.raw`\Delta t \le \frac{\Delta x}{c}, \qquad c = \sqrt{g H} \approx 300\ \mathrm{m\,s^{-1}}, \qquad \Delta t = \frac{86\,400\ \mathrm{s}}{500} \approx 173\ \mathrm{s}`,
  String.raw`\frac{\partial \mathbf{u}}{\partial t} = \cdots - k_v(\sigma)\, \mathbf{u}, \qquad k_v = k_f \max\!\left( 0, \frac{\sigma - \sigma_b}{1 - \sigma_b} \right)`,
  String.raw`T_{\mathrm{eq}} = \max\left\{ 200\,\mathrm{K},\ \left[ 315\,\mathrm{K} - (\Delta T)_y \sin^2 \varphi - (\Delta \theta)_z \log\left( \frac{p}{p_0} \right) \cos^2 \varphi \right] \left( \frac{p}{p_0} \right)^{\kappa} \right\}`,
  String.raw`\frac{\partial u}{\partial t} = -\nu_4 \nabla^4 u, \qquad \nu_4 = \frac{1}{\tau} \left( \frac{\Delta x}{\pi} \right)^4`,
  String.raw`\sum_i A_i \pi_i^{n+1} = \sum_i A_i \pi_i^n \quad \Rightarrow \quad \sum_i A_i \sum_{e \in \mathrm{EC}(i)} n_{e,i}\, F_e = 0`,
  String.raw`\overline{\pi}_e = \tfrac{1}{2} \left( \pi_{c_1(e)} + \pi_{c_2(e)} \right), \qquad F_{e,k} = \overline{\pi}_e\, u_{e,k}\, \Delta\sigma_k`,
  String.raw`\omega = \frac{\mathrm{d} p}{\mathrm{d} t} = \sigma \left( \frac{\partial \pi}{\partial t} + \mathbf{u} \cdot \nabla \pi \right) + \pi \dot\sigma`,
  String.raw`\rho = \frac{p}{R T}, \qquad T = \theta \left( \frac{p}{p_0} \right)^{\kappa}, \qquad \frac{\partial p}{\partial z} = -\rho g`,
];

test('letters are separate italic identifiers', () => {
  assert.equal(body('x'), '<mi>x</mi>');
  assert.equal(body('xy'), '<mi>x</mi><mi>y</mi>');
  assert.equal(body('a b'), '<mi>a</mi><mi>b</mi>');
  assert.equal(body('θ'), '<mi>θ</mi>');
});

test('numbers and decimals are mn', () => {
  assert.equal(body('3.14'), '<mn>3.14</mn>');
  assert.equal(body('12x'), '<mn>12</mn><mi>x</mi>');
  assert.equal(body('1.'), '<mn>1</mn><mo lspace="0em" rspace="0em">.</mo>');
  assert.equal(body('0.286'), '<mn>0.286</mn>');
});

test('plain operators', () => {
  const cases = { '+': '<mo>+</mo>', '-': '<mo>−</mo>', '=': '<mo>=</mo>', '<': '<mo>&lt;</mo>', '>': '<mo>&gt;</mo>', '(': '<mo stretchy="false">(</mo>', ')': '<mo stretchy="false">)</mo>', '[': '<mo stretchy="false">[</mo>', ']': '<mo stretchy="false">]</mo>', ',': '<mo>,</mo>', '.': '<mo lspace="0em" rspace="0em">.</mo>', '/': '<mo stretchy="false" lspace="0em" rspace="0em">/</mo>', '|': '<mo stretchy="false" lspace="0em" rspace="0em">|</mo>', '!': '<mo lspace="0em" rspace="0em">!</mo>', ':': '<mo>:</mo>', ';': '<mo>;</mo>' };
  for (const [tex, mo] of Object.entries(cases)) assert.equal(body(`a${tex}b`), `<mi>a</mi>${mo}<mi>b</mi>`, tex);
  assert.equal(body('a-b'), '<mi>a</mi><mo>\u2212</mo><mi>b</mi>');
  assert.equal(body("f'"), '<mi>f</mi><mo lspace="0em" rspace="0em">′</mo>');
  assert.equal(body("f''"), '<mi>f</mi><mo lspace="0em" rspace="0em">″</mo>');
  assert.equal(body(String.raw`\{x\}`), '<mo stretchy="false">{</mo><mi>x</mi><mo stretchy="false">}</mo>');
  assert.equal(body(String.raw`a \lt b \gt c`), '<mi>a</mi><mo>&lt;</mo><mi>b</mi><mo>&gt;</mo><mi>c</mi>');
});

test('Greek letters: lowercase italic, uppercase upright', () => {
  const lower = { alpha: 'α', beta: 'β', gamma: 'γ', delta: 'δ', epsilon: 'ϵ', varepsilon: 'ε', zeta: 'ζ', eta: 'η', theta: 'θ', vartheta: 'ϑ', iota: 'ι', kappa: 'κ', lambda: 'λ', mu: 'μ', nu: 'ν', xi: 'ξ', pi: 'π', rho: 'ρ', sigma: 'σ', tau: 'τ', upsilon: 'υ', phi: 'ϕ', varphi: 'φ', chi: 'χ', psi: 'ψ', omega: 'ω' };
  const upper = { Gamma: 'Γ', Delta: 'Δ', Theta: 'Θ', Lambda: 'Λ', Xi: 'Ξ', Pi: 'Π', Sigma: 'Σ', Phi: 'Φ', Psi: 'Ψ', Omega: 'Ω' };
  for (const [name, ch] of Object.entries(lower)) assert.equal(body(`\\${name}`), `<mi>${ch}</mi>`, name);
  for (const [name, ch] of Object.entries(upper)) assert.equal(body(`\\${name}`), `<mi mathvariant="normal">${ch}</mi>`, name);
  assert.equal(body(String.raw`\Delta t`), '<mi mathvariant="normal">Δ</mi><mi>t</mi>');
  assert.equal(body(String.raw`\alpha x`), '<mi>α</mi><mi>x</mi>');
  assert.throws(() => texToMathML(String.raw`\alphax`), /unknown command \\alphax/);
});

test('symbols and relations', () => {
  const identifiers = { partial: '∂', nabla: '∇', infty: '∞' };
  const operators = { cdot: '⋅', times: '×', pm: '±', approx: '≈', equiv: '≡', le: '≤', ge: '≥', ne: '≠', to: '→', rightarrow: '→', leftarrow: '←', ldots: '…', cdots: '⋯' };
  for (const [name, ch] of Object.entries(identifiers)) assert.equal(body(`\\${name}`), `<mi mathvariant="normal">${ch}</mi>`, name);
  for (const [name, ch] of Object.entries(operators)) assert.equal(body(`a \\${name} b`), `<mi>a</mi><mo>${ch}</mo><mi>b</mi>`, name);
  assert.equal(body(String.raw`\prime`), '<mo lspace="0em" rspace="0em">′</mo>');
  assert.equal(body(String.raw`x^\prime`), '<msup><mi>x</mi><mo lspace="0em" rspace="0em">′</mo></msup>');
});

test('sums take limits under and over in display, scripts inline; integrals always take scripts', () => {
  const sum = String.raw`\sum_{i=1}^{n} x_i`, limits = '<mrow><mi>i</mi><mo>=</mo><mn>1</mn></mrow><mi>n</mi>';
  assert.equal(display(sum), `<munderover><mo>∑</mo>${limits}</munderover><msub><mi>x</mi><mi>i</mi></msub>`);
  assert.equal(body(sum), `<msubsup><mo>∑</mo>${limits}</msubsup><msub><mi>x</mi><mi>i</mi></msub>`);
  assert.equal(display(String.raw`\sum_e F_e`), '<munder><mo>∑</mo><mi>e</mi></munder><msub><mi>F</mi><mi>e</mi></msub>');
  assert.equal(body(String.raw`\sum_e`), '<msub><mo>∑</mo><mi>e</mi></msub>');
  assert.equal(display(String.raw`\sum^{K}`), '<mover><mo>∑</mo><mi>K</mi></mover>');
  for (const render of [body, display]) assert.equal(render(String.raw`\int_0^1 f`), '<msubsup><mo>∫</mo><mn>0</mn><mn>1</mn></msubsup><mi>f</mi>');
  assert.equal(display(String.raw`\int f`), '<mo>∫</mo><mi>f</mi>');
});

test('fractions', () => {
  assert.equal(body(String.raw`\frac{a}{b}`), '<mfrac><mi>a</mi><mi>b</mi></mfrac>');
  assert.equal(body(String.raw`\frac{a+b}{2}`), '<mfrac><mrow><mi>a</mi><mo>+</mo><mi>b</mi></mrow><mn>2</mn></mfrac>');
  assert.equal(body(String.raw`\frac12`), '<mfrac><mn>1</mn><mn>2</mn></mfrac>');
  assert.equal(body(String.raw`\frac\partial{\partial t}`), '<mfrac><mi mathvariant="normal">∂</mi><mrow><mi mathvariant="normal">∂</mi><mi>t</mi></mrow></mfrac>');
  assert.equal(body(String.raw`\dfrac{a}{b}`), '<mstyle displaystyle="true" scriptlevel="0"><mfrac><mi>a</mi><mi>b</mi></mfrac></mstyle>');
  assert.equal(body(String.raw`\tfrac{1}{2}`), '<mstyle displaystyle="false" scriptlevel="0"><mfrac><mn>1</mn><mn>2</mn></mfrac></mstyle>');
});

test('roots', () => {
  assert.equal(body(String.raw`\sqrt{x}`), '<msqrt><mi>x</mi></msqrt>');
  assert.equal(body(String.raw`\sqrt{g H}`), '<msqrt><mi>g</mi><mi>H</mi></msqrt>');
  assert.equal(body(String.raw`\sqrt2`), '<msqrt><mn>2</mn></msqrt>');
  assert.equal(body(String.raw`\sqrt[3]{x}`), '<mroot><mi>x</mi><mn>3</mn></mroot>');
  assert.equal(body(String.raw`\sqrt[n+1]{x+y}`), '<mroot><mrow><mi>x</mi><mo>+</mo><mi>y</mi></mrow><mrow><mi>n</mi><mo>+</mo><mn>1</mn></mrow></mroot>');
  assert.equal(body(String.raw`\sqrt{[a]}`), '<msqrt><mo stretchy="false">[</mo><mi>a</mi><mo stretchy="false">]</mo></msqrt>');
});

test('subscripts and superscripts', () => {
  assert.equal(body('x^2'), '<msup><mi>x</mi><mn>2</mn></msup>');
  assert.equal(body('x_i'), '<msub><mi>x</mi><mi>i</mi></msub>');
  assert.equal(body('x_i^2'), '<msubsup><mi>x</mi><mi>i</mi><mn>2</mn></msubsup>');
  assert.equal(body('x^2_i'), '<msubsup><mi>x</mi><mi>i</mi><mn>2</mn></msubsup>');
  assert.equal(body('x ^ 2 _ i'), '<msubsup><mi>x</mi><mi>i</mi><mn>2</mn></msubsup>');
  assert.equal(body('x^{10}'), '<msup><mi>x</mi><mn>10</mn></msup>');
  assert.equal(body('x^23'), '<msup><mi>x</mi><mn>2</mn></msup><mn>3</mn>');
  assert.equal(body('x^ab'), '<msup><mi>x</mi><mi>a</mi></msup><mi>b</mi>');
  assert.equal(body('x_{i,k+1/2}'), '<msub><mi>x</mi><mrow><mi>i</mi><mo>,</mo><mi>k</mi><mo>+</mo><mn>1</mn><mo stretchy="false" lspace="0em" rspace="0em">/</mo><mn>2</mn></mrow></msub>');
  assert.equal(body(String.raw`x^\alpha`), '<msup><mi>x</mi><mi>α</mi></msup>');
  assert.equal(body(String.raw`10^{-3}`), `<msup><mn>10</mn><mrow>${unary('−')}<mn>3</mn></mrow></msup>`);
  assert.equal(body('u^+'), `<msup><mi>u</mi>${unary('+')}</msup>`);
  assert.equal(body('^2'), '<msup><mrow></mrow><mn>2</mn></msup>');
  assert.equal(body('{}_a x'), '<msub><mrow></mrow><mi>a</mi></msub><mi>x</mi>');
  assert.equal(body('(a)_y'), '<mo stretchy="false">(</mo><mi>a</mi><msub><mo stretchy="false">)</mo><mi>y</mi></msub>');
  assert.throws(() => texToMathML('x^2^3'), /double superscript/);
  assert.throws(() => texToMathML('x_a_b'), /double subscript/);
});

test('bold, upright and text', () => {
  assert.equal(body(String.raw`\mathbf{u}`), '<mi mathvariant="bold">\u{1d42e}</mi>');
  assert.equal(body(String.raw`\mathbf u`), '<mi mathvariant="bold">\u{1d42e}</mi>');
  assert.equal(body(String.raw`\mathbf{F}`), '<mi mathvariant="bold">\u{1d405}</mi>');
  assert.equal(body(String.raw`\mathbf{uv}`), '<mrow><mi mathvariant="bold">\u{1d42e}</mi><mi mathvariant="bold">\u{1d42f}</mi></mrow>');
  assert.equal(body(String.raw`\mathbf{2}`), '<mn mathvariant="bold">\u{1d7d0}</mn>');
  assert.equal(body(String.raw`\mathbf{\Omega}`), '<mi mathvariant="bold">\u{1d6c0}</mi>');
  assert.equal(body(String.raw`\mathbf{\omega}`), '<mi mathvariant="bold-italic">\u{1d74e}</mi>');
  assert.equal(body(String.raw`\boldsymbol{\omega}`), '<mi mathvariant="bold-italic">\u{1d74e}</mi>');
  assert.equal(body(String.raw`\boldsymbol{\alpha}`), '<mi mathvariant="bold-italic">\u{1d736}</mi>');
  assert.equal(body(String.raw`\boldsymbol{\varphi}`), '<mi mathvariant="bold-italic">\u{1d74b}</mi>');
  assert.equal(body(String.raw`\boldsymbol{\phi}`), '<mi mathvariant="bold-italic">\u{1d753}</mi>');
  assert.equal(body(String.raw`\boldsymbol{\Omega}`), '<mi mathvariant="bold">\u{1d6c0}</mi>');
  assert.equal(body(String.raw`\boldsymbol{x}`), '<mi mathvariant="bold-italic">\u{1d499}</mi>');
  assert.equal(body(String.raw`\boldsymbol{\nabla}`), '<mi mathvariant="bold">\u{1d6c1}</mi>');
  assert.equal(body(String.raw`\mathbf{u}_e`), '<msub><mi mathvariant="bold">\u{1d42e}</mi><mi>e</mi></msub>');
  assert.equal(body(String.raw`\mathrm{d}`), '<mi mathvariant="normal">d</mi>');
  assert.equal(body(String.raw`\mathrm{d}x`), '<mi mathvariant="normal">d</mi><mi>x</mi>');
  assert.equal(body(String.raw`T_{\mathrm{eq}}`), '<msub><mi>T</mi><mi>eq</mi></msub>');
  assert.equal(body(String.raw`\mathrm{m\,s^{-1}}`), `<mrow><mi mathvariant="normal">m</mi>${THIN}<msup><mi mathvariant="normal">s</mi><mrow>${unary('−')}<mn>1</mn></mrow></msup></mrow>`);
  assert.equal(body(String.raw`\mathrm{\mu m}`), '<mrow><mi mathvariant="normal">μ</mi><mi mathvariant="normal">m</mi></mrow>');
  assert.equal(body(String.raw`\operatorname{div} \mathbf{u}`), `<mi>div</mi>${APPLY}${THIN}<mi mathvariant="bold">\u{1d42e}</mi>`);
  assert.equal(body(String.raw`\operatorname{d}`), `<mi mathvariant="normal">d</mi>${APPLY}`);
  assert.equal(body(String.raw`\text{for all } x`), '<mtext>for all\u00a0</mtext><mi>x</mi>');
  assert.equal(body(String.raw`\text{ 50\% }`), '<mtext>\u00a050%\u00a0</mtext>');
  assert.equal(body(String.raw`\text{a {b} c}`), '<mtext>a {b} c</mtext>');
});

test('named functions are upright and followed by a thin space before an operand', () => {
  for (const name of ['ln', 'log', 'exp', 'sin', 'cos', 'tan', 'max', 'min']) assert.equal(body(`\\${name} x`), `<mi>${name}</mi>${APPLY}${THIN}<mi>x</mi>`, name);
  assert.equal(body(String.raw`\sin(x)`), `<mi>sin</mi>${APPLY}<mo stretchy="false">(</mo><mi>x</mi><mo stretchy="false">)</mo>`);
  assert.equal(body(String.raw`\ln \pi`), `<mi>ln</mi>${APPLY}${THIN}<mi>π</mi>`);
  assert.equal(body(String.raw`\sin^2\varphi`), `<msup><mi>sin</mi><mn>2</mn></msup>${APPLY}${THIN}<mi>φ</mi>`);
  assert.equal(body(String.raw`\log_2 x`), `<msub><mi>log</mi><mn>2</mn></msub>${APPLY}${THIN}<mi>x</mi>`);
  assert.equal(display(String.raw`\max_i u_i`), `<munder><mi>max</mi><mi>i</mi></munder>${APPLY}${THIN}<msub><mi>u</mi><mi>i</mi></msub>`);
  assert.equal(body(String.raw`\max_i u_i`), `<msub><mi>max</mi><mi>i</mi></msub>${APPLY}${THIN}<msub><mi>u</mi><mi>i</mi></msub>`);
  assert.equal(body(String.raw`\exp\left(x\right)`), `<mi>exp</mi>${APPLY}${THIN}<mrow><mo fence="true" form="prefix" stretchy="true">(</mo><mi>x</mi><mo fence="true" form="postfix" stretchy="true">)</mo></mrow>`);
  assert.equal(body(String.raw`\sin\,x`), `<mi>sin</mi>${APPLY}${THIN}${THIN}<mi>x</mi>`);
  assert.equal(body(String.raw`\sin\, (x)`), `<mi>sin</mi>${APPLY}${THIN}<mo stretchy="false">(</mo><mi>x</mi><mo stretchy="false">)</mo>`);
  assert.equal(body(String.raw`\sin -x`), `<mi>sin</mi>${APPLY}${THIN}${unary('−')}<mi>x</mi>`);
  assert.equal(body(String.raw`\max\!\left(x\right)`), `<mi>max</mi>${APPLY}${THIN}${BACK}<mrow><mo fence="true" form="prefix" stretchy="true">(</mo><mi>x</mi><mo fence="true" form="postfix" stretchy="true">)</mo></mrow>`);
  assert.equal(body(String.raw`\sin = 1`), `<mi>sin</mi>${APPLY}<mo>=</mo><mn>1</mn>`);
  assert.equal(body(String.raw`{\cos}`), `<mrow><mi>cos</mi>${APPLY}</mrow>`);
});

test('accents', () => {
  const accents = { dot: ['˙', false], ddot: ['¨', false], hat: ['ˆ', false], bar: ['¯', false], overline: ['‾', true], tilde: ['˜', false], vec: ['→', false] };
  for (const [name, [ch, stretchy]] of Object.entries(accents)) assert.equal(body(`\\${name}{u}`), `<mover accent="true"><mi>u</mi><mo stretchy="${stretchy}">${ch}</mo></mover>`, name);
  assert.equal(body(String.raw`\dot\sigma`), '<mover accent="true"><mi>σ</mi><mo stretchy="false">˙</mo></mover>');
  assert.equal(body(String.raw`\overline{\pi}_e`), '<msub><mover accent="true"><mi>π</mi><mo stretchy="true">‾</mo></mover><mi>e</mi></msub>');
  assert.equal(body(String.raw`\hat{\mathbf{k}}`), '<mover accent="true"><mi mathvariant="bold">\u{1d424}</mi><mo stretchy="false">ˆ</mo></mover>');
  assert.equal(body(String.raw`\overline{u v}`), '<mover accent="true"><mrow><mi>u</mi><mi>v</mi></mrow><mo stretchy="true">‾</mo></mover>');
});

test('stretchy fences', () => {
  const fence = (c, form) => `<mo fence="true" form="${form}" stretchy="true">${c}</mo>`, frac = '<mfrac><mi>a</mi><mi>b</mi></mfrac>';
  assert.equal(body(String.raw`\left( \frac{a}{b} \right)`), `<mrow>${fence('(', 'prefix')}${frac}${fence(')', 'postfix')}</mrow>`);
  assert.equal(body(String.raw`\left[ \frac{a}{b} \right]`), `<mrow>${fence('[', 'prefix')}${frac}${fence(']', 'postfix')}</mrow>`);
  assert.equal(body(String.raw`\left| \frac{a}{b} \right|`), `<mrow>${fence('|', 'prefix')}${frac}${fence('|', 'postfix')}</mrow>`);
  assert.equal(body(String.raw`\left\{ \frac{a}{b} \right\}`), `<mrow>${fence('{', 'prefix')}${frac}${fence('}', 'postfix')}</mrow>`);
  assert.equal(body(String.raw`\left. \frac{a}{b} \right|_{0}`), `<msub><mrow>${frac}${fence('|', 'postfix')}</mrow><mn>0</mn></msub>`);
  assert.equal(body(String.raw`\left( x \right.`), `<mrow>${fence('(', 'prefix')}<mi>x</mi></mrow>`);
  assert.equal(body(String.raw`\left( \left[ x \right] \right)^2`), `<msup><mrow>${fence('(', 'prefix')}<mrow>${fence('[', 'prefix')}<mi>x</mi>${fence(']', 'postfix')}</mrow>${fence(')', 'postfix')}</mrow><mn>2</mn></msup>`);
  assert.equal(body(String.raw`a \rightarrow b`), '<mi>a</mi><mo>→</mo><mi>b</mi>');
});

test('spacing', () => {
  const spaces = { ',': '0.1667em', ':': '0.2222em', ';': '0.2778em', quad: '1em', qquad: '2em' };
  for (const [name, width] of Object.entries(spaces)) assert.equal(body(`a\\${name} b`), `<mi>a</mi><mspace width="${width}"></mspace><mi>b</mi>`, name);
  assert.equal(body(String.raw`a\! b`), `<mi>a</mi>${BACK}<mi>b</mi>`);
  assert.equal(body(String.raw`a\ b`), '<mi>a</mi><mspace width="0.3333em"></mspace><mi>b</mi>');
  assert.equal(body('a~b'), '<mi>a</mi><mspace width="0.3333em"></mspace><mi>b</mi>');
});

test('classes wrap a term in an mrow CSS can color', () => {
  assert.equal(body(String.raw`\class{term-a}{x + y}`), '<mrow class="term-a"><mi>x</mi><mo>+</mo><mi>y</mi></mrow>');
  assert.equal(body(String.raw`\class{term-b}{x}`), '<mrow class="term-b"><mi>x</mi></mrow>');
  assert.equal(body(String.raw`\class{term-c wide}{x}^2`), '<msup><mrow class="term-c wide"><mi>x</mi></mrow><mn>2</mn></msup>');
  assert.equal(body(String.raw`-\class{term-d}{\nabla K}`), `${unary('−')}<mrow class="term-d"><mi mathvariant="normal">∇</mi><mi>K</mi></mrow>`);
  assert.throws(() => texToMathML(String.raw`\class{x" onclick="y}{a}`), /bad class name/);
  assert.throws(() => texToMathML(String.raw`\class{}{a}`), /bad class name/);
});

test('aligned becomes a table with alternating right and left columns', () => {
  const cell = (side, content) => `<mtd columnalign="${side}">${side === 'left' ? '<mi></mi>' : ''}${content}</mtd>`;
  const table = (align, rows) => `<mtable class="tex-aligned" displaystyle="true" columnalign="${align}" columnspacing="0em">${rows.map((r) => `<mtr>${r}</mtr>`).join('')}</mtable>`;
  assert.equal(display(String.raw`\begin{aligned} a &= b \\ c &= d + e \end{aligned}`), table('right left', [cell('right', '<mi>a</mi>') + cell('left', '<mo>=</mo><mi>b</mi>'), cell('right', '<mi>c</mi>') + cell('left', '<mo>=</mo><mi>d</mi><mo>+</mo><mi>e</mi>')]));
  assert.equal(display(String.raw`\begin{aligned} a &= b \\ \end{aligned}`), table('right left', [cell('right', '<mi>a</mi>') + cell('left', '<mo>=</mo><mi>b</mi>')]));
  assert.equal(display(String.raw`\begin{aligned} a &= b & c &= d \end{aligned}`), table('right left right left', [cell('right', '<mi>a</mi>') + cell('left', '<mo>=</mo><mi>b</mi>') + cell('right', '<mi>c</mi>') + cell('left', '<mo>=</mo><mi>d</mi>')]));
  assert.equal(display(String.raw`\begin{aligned} a &= b \\ &\quad + c \end{aligned}`), table('right left', [cell('right', '<mi>a</mi>') + cell('left', '<mo>=</mo><mi>b</mi>'), cell('right', '') + cell('left', '<mspace width="1em"></mspace><mo>+</mo><mi>c</mi>')]));
  assert.equal(display(String.raw`\begin{aligned} a &= -b \\ &\quad - c \end{aligned}`), table('right left', [cell('right', '<mi>a</mi>') + cell('left', `<mo>=</mo>${unary('−')}<mi>b</mi>`), cell('right', '') + cell('left', '<mspace width="1em"></mspace><mo>−</mo><mi>c</mi>')]));
  assert.equal(display(String.raw`\begin{aligned} -a &= b \end{aligned}`), table('right left', [cell('right', `${unary('−')}<mi>a</mi>`) + cell('left', '<mo>=</mo><mi>b</mi>')]));
  const sums = display(String.raw`\begin{aligned} x &= \sum_{e} F_e \end{aligned}`);
  assert.match(sums, /<munder><mo>∑<\/mo><mi>e<\/mi><\/munder>/);
});

test('a binary operator is unary at the start and after a relation, opening, comma, binary operator or big operator', () => {
  assert.equal(body('a - b'), '<mi>a</mi><mo>−</mo><mi>b</mi>');
  assert.equal(body('-a'), `${unary('−')}<mi>a</mi>`);
  assert.equal(body('x = -y'), `<mi>x</mi><mo>=</mo>${unary('−')}<mi>y</mi>`);
  assert.equal(body(String.raw`x \le +y`), `<mi>x</mi><mo>≤</mo>${unary('+')}<mi>y</mi>`);
  assert.equal(body(String.raw`\sigma = \pm 1`), `<mi>σ</mi><mo>=</mo>${unary('±')}<mn>1</mn>`);
  assert.equal(body('(-x)'), `<mo stretchy="false">(</mo>${unary('−')}<mi>x</mi><mo stretchy="false">)</mo>`);
  assert.equal(body('a, -b'), `<mi>a</mi><mo>,</mo>${unary('−')}<mi>b</mi>`);
  assert.equal(body('a + -b'), `<mi>a</mi><mo>+</mo>${unary('−')}<mi>b</mi>`);
  assert.equal(body('--a'), `${unary('−')}<mo>−</mo><mi>a</mi>`);
  assert.equal(body(String.raw`\left( -x \right) - y`), `<mrow><mo fence="true" form="prefix" stretchy="true">(</mo>${unary('−')}<mi>x</mi><mo fence="true" form="postfix" stretchy="true">)</mo></mrow><mo>−</mo><mi>y</mi>`);
  assert.equal(body(String.raw`\sum_i -x_i`), `<msub><mo>∑</mo><mi>i</mi></msub>${unary('−')}<msub><mi>x</mi><mi>i</mi></msub>`);
  assert.equal(body(String.raw`= \, -x`), `<mo>=</mo>${THIN}${unary('−')}<mi>x</mi>`);
  assert.equal(body(String.raw`\frac{a}{b} - c`), '<mfrac><mi>a</mi><mi>b</mi></mfrac><mo>−</mo><mi>c</mi>');
  assert.equal(body(String.raw`x \cdots - y`), '<mi>x</mi><mo>⋯</mo><mo>−</mo><mi>y</mi>');
});

test('braces group without adding markup for a single item', () => {
  assert.equal(body('{a+b}^2'), '<msup><mrow><mi>a</mi><mo>+</mo><mi>b</mi></mrow><mn>2</mn></msup>');
  assert.equal(body('{x}'), '<mi>x</mi>');
  assert.equal(body('{}'), '<mrow></mrow>');
  assert.equal(body('a{=}b'), '<mi>a</mi><mo>=</mo><mi>b</mi>');
});

test('the wrapper carries the escaped source', () => {
  assert.equal(texToMathML('x', { display: true }), '<math display="block"><semantics><mrow><mi>x</mi></mrow><annotation encoding="application/x-tex">x</annotation></semantics></math>');
  assert.equal(texToMathML('x'), '<math display="inline"><semantics><mrow><mi>x</mi></mrow><annotation encoding="application/x-tex">x</annotation></semantics></math>');
  assert.equal(texToMathML(''), '<math display="inline"><semantics><mrow></mrow><annotation encoding="application/x-tex"></annotation></semantics></math>');
});

test('<, > and & are escaped everywhere', () => {
  assert.equal(body('a<b>c'), '<mi>a</mi><mo>&lt;</mo><mi>b</mi><mo>&gt;</mo><mi>c</mi>');
  assert.equal(body(String.raw`\text{<b>&amp;"x"</b>}`), '<mtext>&lt;b&gt;&amp;amp;&quot;x&quot;&lt;/b&gt;</mtext>');
  assert.equal(body(String.raw`a \& b`), '<mi>a</mi><mo>&amp;</mo><mi>b</mi>');
  const markup = texToMathML(String.raw`\begin{aligned} a &< b \\ c &> d \end{aligned} `, { display: true });
  assert.match(markup, /<annotation encoding="application\/x-tex">\\begin\{aligned\} a &amp;&lt; b \\\\ c &amp;&gt; d \\end\{aligned\} <\/annotation>/);
  wellFormed(String.raw`\text{</math><script>alert(1)</script>} < b >`);
  assert.doesNotMatch(texToMathML(String.raw`\text{</math><script>}`), /<script>|<\/math><script/);
});

test('bad TeX throws a SyntaxError that says what went wrong', () => {
  const cases = [
    [String.raw`\foo`, /unknown command \\foo/], [String.raw`\toString`, /unknown command/], ['{x', /missing \}/], ['x}', /unmatched \}/],
    [String.raw`\left( x`, /missing \\right/], [String.raw`x \right)`, /unexpected \\right/], ['a & b', /& outside aligned/], [String.raw`a \\ b`, /unexpected \\\\/],
    [String.raw`\begin{matrix} a \end{matrix}`, /unsupported environment matrix/], [String.raw`\begin{aligned} a \end{array}`, /mismatched \\end/], [String.raw`\begin{aligned} a &= b`, /missing \\end\{aligned\}/],
    [String.raw`\left< x \right>`, /bad delimiter </], [String.raw`\left`, /bad delimiter at the end/], [String.raw`\frac{a}`, /missing argument/], ['x^', /missing argument/], ['x^}', /missing argument/],
    ['\\', /lone backslash/], [String.raw`\text{\alpha}`, /commands inside \\text/], [String.raw`\text{a`, /missing \}/], ['#', /unexpected #/], [String.raw`\sqrt[3{x}`, /missing \]/], [String.raw`\frac{a}{b`, /missing \}/],
  ];
  for (const [tex, message] of cases) assert.throws(() => texToMathML(tex), (error) => error instanceof SyntaxError && message.test(error.message), tex);
});

test('the well-formedness check rejects broken markup', () => {
  for (const markup of ['<mi>x</mo>', '<mrow><mi>x</mi>', '<mi>a & b</mi>', '<mi>a > b</mi>', '<mi x=1>a</mi>']) assert.throws(() => tree(markup), markup);
  assert.throws(() => check(tree('<mfrac><mi>a</mi></mfrac>')[0]), /2 children/);
  assert.throws(() => check(tree('<mi><mn>1</mn></mi>')[0]), /only text/);
  assert.throws(() => check(tree('<menclose><mi>a</mi></menclose>')[0]), /MathML Core/);
});

test(`a corpus of ${corpus.length} model formulas is well formed in both modes`, () => {
  assert.ok(corpus.length >= 25);
  for (const tex of corpus) for (const mode of [false, true]) {
    const math = wellFormed(tex, { display: mode });
    assert.equal(math.name, 'math');
  }
  const momentum = texToMathML(corpus[7], { display: true });
  for (const term of ['a', 'b', 'c', 'd', 'e']) assert.match(momentum, new RegExp(`<mrow class="term-${term}">`));
  const rk4 = texToMathML(corpus[22], { display: true });
  assert.equal(rk4.match(/<mtr>/g).length, 5);
  assert.equal(rk4.match(/<mtd columnalign="right">/g).length, 5);
  assert.equal(rk4.match(/<mtd columnalign="left"><mi><\/mi><mo>=<\/mo>/g).length, 5);
  const discrete = texToMathML(corpus[12], { display: true });
  assert.equal(discrete.match(/<munder><mo>∑<\/mo>/g).length, 2);
  assert.equal(discrete.match(/<munderover><mo>∑<\/mo>/g).length, 1);
  assert.equal(texToMathML(corpus[12]).match(/<msubsup><mo>∑<\/mo>/g).length, 1);
});

function fakeDom(specs) {
  const elements = specs.map(([localName, text, classes = []]) => {
    const list = new Set(['tex', ...classes]), attributes = new Map();
    return {
      localName, textContent: text, html: null,
      classList: { contains: (c) => list.has(c), add: (c) => list.add(c), has: (c) => list.has(c) },
      hasAttribute: (n) => attributes.has(n), setAttribute: (n, v) => attributes.set(n, String(v)), getAttribute: (n) => attributes.get(n) ?? null,
      replaceChildren(...nodes) { this.html = nodes.map((n) => n.markup).join(''); this.textContent = ''; },
    };
  });
  const root = { queried: [], querySelectorAll(selector) { this.queried.push(selector); return elements; } };
  return { root, elements };
}

test('renderMath replaces TeX with MathML once and marks bad TeX', (t) => {
  const warnings = [];
  t.mock.method(console, 'warn', (message) => warnings.push(message));
  globalThis.document = { createElement: (tag) => { assert.equal(tag, 'template'); return { set innerHTML(markup) { this.content = { markup }; } }; } };
  t.after(() => { delete globalThis.document; });
  const { root, elements: [inline, block, flagged, broken] } = fakeDom([['span', ' x^2 '], ['div', String.raw`\frac{a}{b}`], ['span', 'y', ['display']], ['span', String.raw`\frac{a}`]]);
  renderMath(root);
  assert.deepEqual(root.queried, ['.tex']);
  assert.equal(inline.html, texToMathML('x^2'));
  assert.equal(block.html, texToMathML(String.raw`\frac{a}{b}`, { display: true }));
  assert.equal(flagged.html, texToMathML('y', { display: true }));
  assert.equal(broken.html, null);
  assert.equal(broken.textContent, String.raw`\frac{a}`);
  assert.ok(broken.classList.has('tex-error'));
  assert.equal(broken.getAttribute('data-math'), 'error');
  assert.equal(inline.getAttribute('data-math'), 'rendered');
  assert.equal(warnings.length, 1);
  assert.match(warnings[0], /\\frac\{a\}.*missing argument/);
  inline.html = 'kept';
  renderMath(root);
  assert.equal(inline.html, 'kept');
  assert.equal(warnings.length, 1);
});

test('a sign opening a colored term is binary after an operand and unary after a relation', () => {
  assert.equal(body(String.raw`a \class{term-b}{- b}`), '<mi>a</mi><mrow class="term-b"><mo form="infix">−</mo><mi>b</mi></mrow>');
  assert.equal(body(String.raw`= \class{term-b}{- b}`), `<mo>=</mo><mrow class="term-b">${unary('−')}<mi>b</mi></mrow>`);
});

test('every equation on the explainer pages is well formed', () => {
  let count = 0;
  for (const page of ['air', 'grid', 'sphere', 'primitive']) {
    const html = readFileSync(new URL(`../explain/${page}.html`, import.meta.url), 'utf8');
    for (const [, tag, flag, source] of html.matchAll(/<(span|div) class="tex( display)?">([^]*?)<\/\1>/g)) {
      wellFormed(unescape(source).trim(), tag === 'div' || flag ? { display: true } : undefined);
      count++;
    }
  }
  assert.ok(count > 0);
});
