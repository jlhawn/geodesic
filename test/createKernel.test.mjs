import { test } from 'node:test';
import assert from 'node:assert/strict';
import { createKernel } from '../js/gpu/device.module.js';

function fakeDevice() {
  let made = 0;
  return {
    made: () => made,
    createShaderModule: ({ code }) => ({ code }),
    createComputePipeline: () => ({ getBindGroupLayout: () => ({}) }),
    createBindGroup: ({ entries }) => { made++; return { buffers: entries.map((entry) => entry.resource.buffer) }; },
  };
}

function boundBy(kernel, buffers) {
  let group = null;
  kernel.run({ setPipeline() {}, setBindGroup: (_, g) => { group = g; }, dispatchWorkgroups() {} }, buffers, 1);
  return group.buffers;
}

test('a createKernel kernel binds each run\'s own buffers, unlabeled or sharing labels, and reuses a bind group for the same buffers', () => {
  const device = fakeDevice(), kernel = createKernel(device, '', 'probe');
  const sets = [
    [{ label: '' }, { label: '' }],
    [{ label: '' }, { label: '' }],
    [{ label: 'A' }, { label: 'B' }],
    [{ label: 'A' }, { label: 'B' }],
  ];
  for (const set of sets) {
    const bound = boundBy(kernel, set);
    assert.equal(bound.length, set.length);
    set.forEach((buffer, binding) => assert.equal(bound[binding], buffer, `binding ${binding} of a run`));
  }
  assert.equal(device.made(), sets.length);
  for (const set of sets) boundBy(kernel, set);
  assert.equal(device.made(), sets.length, 'a buffer set seen before keeps its bind group');
});
