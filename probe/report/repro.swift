// Pure-Metal reproduction, no WebGPU involved:
//   swift repro.swift min_robust_on.msl     (expected: many wrong values on an Apple M5)
//   swift repro.swift min_robust_off.msl    (the same shader without index clamps: 0 wrong)
// Compiles the given Metal source at run time, runs the kernel over 4096 threads,
// and compares the output with the running product computed on the CPU in Float32.
import Foundation
import Metal

let path = CommandLine.arguments.count > 1 ? CommandLine.arguments[1] : "min_robust_on.msl"
let source = try String(contentsOfFile: path, encoding: .utf8)
let device = MTLCreateSystemDefaultDevice()!
print("device:", device.name)
let library = try device.makeLibrary(source: source, options: MTLCompileOptions())
let function = library.makeFunction(name: "dawn_entry_point")!
let pipeline = try device.makeComputePipelineState(function: function)

let N = 4096, M = 15
var seed: UInt32 = 12345
func rand() -> Float { seed = seed &* 1664525 &+ 1013904223; return Float(Double(seed) / 4294967296.0) }
var L = [Float](repeating: 0, count: M * N)
for x in 0..<(M * N) { L[x] = 0.5 + rand() }
var expected = [Float](repeating: 0, count: M * N)
for i in 0..<N { var g = L[i]; expected[i] = g; for j in 1..<M { g = L[j * N + i] * g; expected[j * N + i] = g } }

let bytes = M * N * MemoryLayout<Float>.size
let bufL = device.makeBuffer(bytes: L, length: bytes, options: .storageModeShared)!
let bufOut = device.makeBuffer(length: bytes, options: .storageModeShared)!
var immediate: [UInt32] = [0, UInt32(bytes), UInt32(bytes)]   // Dawn's immediate data: a zero, then each storage buffer's size in bytes
let bufImm = device.makeBuffer(bytes: &immediate, length: 12, options: .storageModeShared)!

let queue = device.makeCommandQueue()!
let commands = queue.makeCommandBuffer()!
let encoder = commands.makeComputeCommandEncoder()!
encoder.setComputePipelineState(pipeline)
encoder.setBuffer(bufL, offset: 0, index: 0)
encoder.setBuffer(bufOut, offset: 0, index: 1)
encoder.setBuffer(bufImm, offset: 0, index: 30)
encoder.dispatchThreadgroups(MTLSize(width: N / 64, height: 1, depth: 1), threadsPerThreadgroup: MTLSize(width: 64, height: 1, depth: 1))
encoder.endEncoding()
commands.commit()
commands.waitUntilCompleted()

let out = bufOut.contents().bindMemory(to: Float.self, capacity: M * N)
var wrong = 0, first = -1
for x in 0..<(M * N) { if out[x] != expected[x] { wrong += 1; if first < 0 { first = x } } }
if first >= 0 {
  print("\(path): \(wrong) of \(M * N) wrong; first at thread \(first % N), row \(first / N): got \(out[first]), expected \(expected[first])")
} else {
  print("\(path): 0 of \(M * N) wrong")
}
