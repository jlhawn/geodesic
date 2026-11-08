import { test } from 'node:test';
import assert from 'node:assert/strict';
import { mkdtempSync, mkdirSync, writeFileSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { terrainCache } from '../js/geography.module.js';

test('the fine grids are found beside the repository a git worktree belongs to, wherever the worktree is checked out', () => {
  const held = process.env.TERRAIN_CACHE, home = mkdtempSync(join(tmpdir(), 'terrain-home-')), away = mkdtempSync(join(tmpdir(), 'terrain-away-'));
  delete process.env.TERRAIN_CACHE;
  try {
    const repository = join(home, 'geodesic'), grids = join(home, 'geodesic-terrain-cache'), worktree = join(away, 'review'), gitdir = join(repository, '.git', 'worktrees', 'review');
    mkdirSync(gitdir, { recursive: true });
    mkdirSync(grids);
    mkdirSync(worktree);
    for (const name of ['h5', 'flt2']) writeFileSync(join(grids, `fine_${name}.f32`), '');
    assert.equal(terrainCache(repository), grids);
    assert.equal(terrainCache(worktree), null);
    writeFileSync(join(worktree, '.git'), `gitdir: ${gitdir}\n`);
    assert.equal(terrainCache(worktree), grids);
  } finally {
    if (held === undefined) delete process.env.TERRAIN_CACHE; else process.env.TERRAIN_CACHE = held;
    rmSync(home, { recursive: true, force: true });
    rmSync(away, { recursive: true, force: true });
  }
});
