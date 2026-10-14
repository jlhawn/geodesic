# Metal compiler on Apple M5: a down-counting loop over a thread-private array of more than 32 floats reads back zeros when its indices are clamped

## Apple Feedback (Developer Technologies → Metal)

**Title.** Apple M5 Metal compiler: a down-counting loop over a thread-private array of more than 32 floats reads back the array's zero initializer when the indices are clamped with `min(uint(j), N - 1u)` and the loop bounds fold to constants after inlining.

**Environment.** Apple M5 Max (MacBook Pro), macOS 26.7.1 (25G241), Metal 3. The same shader is correct on an Apple M1 Max running the same macOS. Reproduced in Google Chrome (WebGPU), in Node.js through Google's Dawn (the `webgpu` npm package 0.6.1), and with a plain Metal host that compiles the Metal source at run time (attached), so no WebGPU implementation is involved.

**Steps to reproduce.** Unzip the attachment and run, in Terminal:

    swift repro.swift min_robust_on.msl
    swift repro.swift min_robust_off.msl

`repro.swift` compiles the given Metal source with `MTLDevice.makeLibrary(source:options:)`, runs the kernel `dawn_entry_point` over 4096 threads in threadgroups of 64, and compares the output buffer with the same running product computed on the CPU in Float32.

**Expected.** Both files hold the same algorithm: thread `i` fills a thread-private `float gain[34]` with a running product over 15 rows from an input buffer, then copies the 15 values to the output buffer in a loop that counts down from row 14 to row 0. The output should equal the CPU's running product in every element. `min_robust_off.msl` does.

**Actual.** `min_robust_on.msl` returns 0 for rows 0 to 13 of every thread and the right value for row 14 only: 57,344 of 61,440 values wrong. The two files differ only in the index clamps that `min_robust_on.msl` carries (`gain[min(uint(j), 33u)]` and the equivalent clamps on the device buffers), which are no-ops here since every index is in range.

The failing loop, as compiled:

    {
      int j = as_type<int>((as_type<uint>(m) - 1u));      // m is 15 after inlining
      while(true) {
        if ((j >= 0)) { } else { break; }
        device float* const v_2 = (&(*tint_module_vars.OUT)[min(uint(...j * 4096 + i...), (tint_array_lengths.tint_array_length_0_1 - 1u))]);
        (*v_2) = gain[min(uint(j), 33u)];
        { j = as_type<int>((as_type<uint>(j) - 1u)); }
      }
    }

**What makes it happen.** All four together:

1. a thread-private array of more than 32 floats (32 is fine; 34 fails);
2. a forward loop that fills it with a recurrence (each element from the previous one);
3. a loop that then counts down to read it back;
4. loop bounds that are constants only after inlining (`solve(i, 15)`), with the clamped indices above.

Change any one of them and the result is right: an array of 32 elements, a read-back loop that counts up over the same indices, a bound only known at run time, or no index clamps. The source is valid Metal and every index is in bounds, so the output should not depend on the clamp.

**Why it matters.** The clamps are what WebGPU implementations insert for the specification's bounds-safety guarantee (Chrome's shader compiler, Tint, emits exactly this form), so any WebGPU compute shader with a tridiagonal back-substitution or a similar backward sweep over a private array is at risk on this GPU. It cost our model its boundary-layer mixing: a climate model's clouds disappeared on the M5 in Chrome, while Safari, which compiles WGSL differently, was unaffected.

**Attachments.**
- `repro.swift`: the Metal host (reads the `.msl` path from the command line).
- `min_robust_on.msl`: the failing shader, as Tint generated it with index clamps.
- `min_robust_off.msl`: the same shader without the clamps, which is correct.
- `repro_min.mjs`: the original WebGPU reproduction (Node.js, the `webgpu` package), with the WGSL source inline.

## Chromium / Dawn issue (issues.chromium.org, component Blink>WebGPU), for awareness

**Title.** Apple M5: Tint's robustness clamps trigger a Metal compiler bug in down-counting loops over private arrays longer than 32 elements.

**Summary.** On an Apple M5 Max (macOS 26.7.1), a WGSL compute shader that fills `var gain: array<f32, 34>` with a recurrence and reads it back with `for (var j = m - 1; j >= 0; j--)` returns zeros for all but the last element when robustness is on (the default). With `--enable-dawn-features=disable_robustness` the result is correct, as it is on an M1 Max, and in Safari. The generated MSL is valid and in bounds; the fault is Apple's compiler (Feedback filed, number to be added). Minimal repro attached (`repro_min.mjs`, runnable with the `webgpu` npm package; the WGSL is inline). A possible mitigation on Tint's side: skip the clamp when the index is provably in range, which here it is once the bound folds to a constant after inlining (the `enable_integer_range_analysis_in_robustness` toggle may already cover the plain case; this one needs the constant to propagate through the call). Workaround in our code: write the back-substitution as an up-counting loop over the same indices.
