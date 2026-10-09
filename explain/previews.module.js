const cache = new Map();
let card = null, current = null, showTimer = 0, hideTimer = 0, sheet = null, backdrop = null, sheetLink = null;

const titleOf = (href) => href.match(/wikipedia\.org\/wiki\/([^#?]+)/)?.[1] ?? null;

function summary(title) {
  if (!cache.has(title)) cache.set(title, fetch(`https://en.wikipedia.org/api/rest_v1/page/summary/${title}`, { headers: { Accept: 'application/json' } }).then((response) => response.ok ? response.json() : null).catch(() => null));
  return cache.get(title);
}

function content(data) {
  const fragment = document.createDocumentFragment();
  if (data.thumbnail) {
    const img = document.createElement('img');
    img.crossOrigin = 'anonymous';
    img.src = data.thumbnail.source;
    img.width = data.thumbnail.width;
    img.height = data.thumbnail.height;
    img.alt = '';
    fragment.append(img);
  }
  const body = document.createElement('div'), title = document.createElement('div'), text = document.createElement('p'), from = document.createElement('div');
  body.className = 'body';
  title.className = 'title';
  title.textContent = data.title;
  text.textContent = data.extract;
  from.className = 'from';
  from.textContent = data.from ?? 'From Wikipedia';
  body.append(title, text, from);
  fragment.append(body);
  return fragment;
}

function element() {
  if (card) return card;
  card = document.createElement('div');
  card.className = 'preview';
  card.setAttribute('role', 'tooltip');
  card.hidden = true;
  card.addEventListener('mouseenter', () => clearTimeout(hideTimer));
  card.addEventListener('mouseleave', scheduleHide);
  document.body.append(card);
  return card;
}

function place(box, link) {
  const rect = link.getBoundingClientRect(), gap = 8, width = Math.min(320, innerWidth - 32);
  box.style.width = `${width}px`;
  const height = box.offsetHeight;
  const below = rect.bottom + gap + height <= innerHeight || rect.top - gap - height < 0;
  box.style.top = `${scrollY + (below ? rect.bottom + gap : rect.top - gap - height)}px`;
  box.style.left = `${scrollX + Math.max(16, Math.min(rect.left, innerWidth - width - 16))}px`;
}

function render(data, link) {
  const box = element();
  box.replaceChildren(content(data));
  box.hidden = false;
  place(box, link);
}

function hide() {
  if (card) card.hidden = true;
  current = null;
}

function scheduleHide() {
  clearTimeout(hideTimer);
  hideTimer = setTimeout(hide, 250);
}

function show(link, delay) {
  const title = titleOf(link.href);
  if (!title) return;
  clearTimeout(hideTimer);
  clearTimeout(showTimer);
  current = link;
  showTimer = setTimeout(async () => { const data = await summary(title); if (data && current === link) render(data, link); }, delay);
}

function sheetElements() {
  if (sheet) return;
  backdrop = document.createElement('div');
  backdrop.className = 'preview-backdrop';
  backdrop.hidden = true;
  backdrop.addEventListener('click', closeSheet);
  sheet = document.createElement('div');
  sheet.className = 'sheet';
  sheet.setAttribute('role', 'dialog');
  sheet.setAttribute('aria-modal', 'true');
  sheet.hidden = true;
  let startY = 0, dy = 0, moved = 0;
  sheet.addEventListener('pointerdown', (e) => { startY = e.clientY; dy = 0; moved = 0; sheet.setPointerCapture(e.pointerId); sheet.style.transition = 'none'; });
  sheet.addEventListener('pointermove', (e) => { if (!sheet.hasPointerCapture(e.pointerId)) return; dy = Math.max(0, e.clientY - startY); moved = Math.max(moved, Math.abs(e.clientY - startY)); sheet.style.transform = `translateY(${dy}px)`; });
  sheet.addEventListener('pointerup', (e) => {
    if (!sheet.hasPointerCapture(e.pointerId)) return;
    sheet.releasePointerCapture(e.pointerId);
    sheet.style.transition = '';
    sheet.style.transform = '';
    if (dy > 80) closeSheet();
    else if (moved < 10 && sheetLink) location.assign(sheetLink.href);
  });
  sheet.addEventListener('pointercancel', () => { sheet.style.transition = ''; sheet.style.transform = ''; });
  sheet.addEventListener('transitionend', () => { if (!sheet.classList.contains('open')) { sheet.hidden = true; backdrop.hidden = true; } });
  document.body.append(backdrop, sheet);
}

async function openSheet(link) {
  sheetElements();
  sheetLink = link;
  const handle = document.createElement('div');
  handle.className = 'handle';
  sheet.replaceChildren(handle, content({ title: link.textContent, extract: 'Loading…', from: 'From Wikipedia · tap to open' }));
  backdrop.hidden = false;
  sheet.hidden = false;
  void sheet.offsetHeight;
  sheet.classList.add('open');
  const data = await summary(titleOf(link.href));
  if (data && sheetLink === link) sheet.replaceChildren(handle, content({ ...data, from: 'From Wikipedia · tap to open' }));
}

function closeSheet() {
  sheetLink = null;
  if (sheet) sheet.classList.remove('open');
}

export function previews(root) {
  const touch = matchMedia('(hover: none)').matches;
  for (const link of root.querySelectorAll('a[href*="wikipedia.org/wiki/"]')) {
    if (touch) { link.addEventListener('click', (e) => { e.preventDefault(); openSheet(link); }); continue; }
    link.addEventListener('mouseenter', () => show(link, 350));
    link.addEventListener('mouseleave', () => { clearTimeout(showTimer); scheduleHide(); });
    link.addEventListener('focus', () => show(link, 0));
    link.addEventListener('blur', scheduleHide);
  }
  addEventListener('keydown', (e) => { if (e.key === 'Escape') { hide(); closeSheet(); } });
  addEventListener('scroll', () => { clearTimeout(showTimer); hide(); }, { passive: true });
}
