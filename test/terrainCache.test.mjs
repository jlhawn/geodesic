import { test } from 'node:test';
import assert from 'node:assert/strict';
import { mkdtempSync, mkdirSync, writeFileSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { terrainCache } from '../js/geography.module.js';

test('the fine grids are found beside the repository a git worktree belongs to, wherever the worktree is checked out, and beside a bare repository', () => {
  const held = process.env.TERRAIN_CACHE, roots = ['home', 'store', 'away'].map((name) => mkdtempSync(join(tmpdir(), `terrain-${name}-`)));
  const [home, store, away] = roots;
  delete process.env.TERRAIN_CACHE;
  try {
    const repository = join(home, 'geodesic'), grids = join(home, 'geodesic-terrain-cache'), worktree = join(away, 'review'), gitdir = join(repository, '.git', 'worktrees', 'review');
    const bareGrids = join(store, 'geodesic-terrain-cache'), bareWorktree = join(away, 'bare'), bareGitdir = join(store, 'geodesic.git', 'worktrees', 'bare');
    for (const dir of [gitdir, grids, worktree, bareGitdir, bareGrids, bareWorktree]) mkdirSync(dir, { recursive: true });
    for (const dir of [grids, bareGrids]) for (const name of ['h5', 'flt2']) writeFileSync(join(dir, `fine_${name}.f32`), '');
    for (const dir of [gitdir, bareGitdir]) writeFileSync(join(dir, 'commondir'), '../..\n');
    assert.equal(terrainCache(repository), grids);
    assert.equal(terrainCache(worktree), null);
    writeFileSync(join(worktree, '.git'), `gitdir: ${gitdir}\n`);
    assert.equal(terrainCache(worktree), grids);
    assert.equal(terrainCache(bareWorktree), null);
    writeFileSync(join(bareWorktree, '.git'), `gitdir: ${bareGitdir}\n`);
    assert.equal(terrainCache(bareWorktree), bareGrids);
  } finally {
    if (held === undefined) delete process.env.TERRAIN_CACHE; else process.env.TERRAIN_CACHE = held;
    for (const root of roots) rmSync(root, { recursive: true, force: true });
  }
});
