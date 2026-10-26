# The GPU queue of scripts/gpuQueue.mjs and the CPU slots of scripts/cpuQueue.mjs for zsh:
#   source scripts/gpuQueue.zsh
#   export GPU_WHO=<session name>                  once per session; labels the runs of both
#   GPU_JOB='<one line>' gpu_run LOG CMD [ARG...]  waits its turn, holds the GPU for all of CMD, exits with its status
#   GPU_PRIORITY=<1-3> gpu_run …                   the run's level, 3 by default: 1: a job under 5 minutes that tolerates
#                                                  contention, such as a GPU test file, up to 5 at once, each on CPU_SLOTS
#                                                  slots (1); 2: a measurement such as pace, at most 10 minutes, alone;
#                                                  3: a long job with heavy GPU use, alone
#   GPU_JOB='<one line>' gpu_share LOG CMD [ARG...]
#                                                  gpu_run at level 1; status 75 when a GPU abort in node fails it
#   gpu_lane [reset]                               prints how many level-1 jobs may share the GPU at once, and why;
#                                                  reset lets 5 share again after an abort dropped level 1 to one at a time
#   gpu_status                                     rewrites and prints <queue dir>/status.md: the holder and the queue in order
#   gpu_serve [PORT]                               serves a live page of the queue and its recent runs on http://127.0.0.1:PORT/ (8008)
#   CPU_SLOTS=<k> CPU_JOB='<one line>' cpu_run LOG CMD [ARG...]
#                                                  waits until k slots of 1 GB are free (1 by default), holds them for all of CMD,
#                                                  exits with its status; a node --test CMD passes --test-concurrency=<at most k>
#   cpu_status                                     prints the runs holding CPU slots and those waiting
# or, without sourcing, zsh scripts/gpuQueue.zsh run LOG CMD [ARG...] | share LOG CMD [ARG...] | lane [reset] | status | serve [PORT]
# and node scripts/cpuQueue.mjs run LOG CMD [ARG...] | status.
typeset -g GPU_QUEUE_SCRIPT=${${(%):-%x}:A:h}/gpuQueue.mjs CPU_QUEUE_SCRIPT=${${(%):-%x}:A:h}/cpuQueue.mjs

_queue_node() {
  local node=${GPU_QUEUE_NODE:-$HOME/.nvm/versions/node/v22.22.2/bin/node}
  [[ -x $node ]] || node=node
  $node "$@"
}
_gpu_queue() { _queue_node $GPU_QUEUE_SCRIPT "$@" }
gpu_run() { _gpu_queue run "$@" }
gpu_share() { _gpu_queue share "$@" }
gpu_lane() { _gpu_queue lane "$@" }
gpu_status() { _gpu_queue status }
gpu_serve() { _gpu_queue serve "$@" }
cpu_run() { _queue_node $CPU_QUEUE_SCRIPT run "$@" }
cpu_status() { _queue_node $CPU_QUEUE_SCRIPT status }

if [[ $ZSH_EVAL_CONTEXT == toplevel ]]; then
  _gpu_queue "$@"
fi
