# The GPU queue of scripts/gpuQueue.mjs and the CPU slots of scripts/cpuQueue.mjs for zsh:
#   source scripts/gpuQueue.zsh
#   export GPU_WHO=<session name>                  once per session; labels the runs of both
#   GPU_JOB='<one line>' gpu_run LOG CMD [ARG...]  waits its turn, holds the GPU for all of CMD, exits with its status
#   GPU_PRIORITY=<1-5> gpu_run …                   the run's level, 1 the highest and 3 the default: 1: a check Joshua is waiting on;
#                                                  2: test suites before a merge or push; 3: the default; 4: routine reports and figures;
#                                                  5: long background experiments, in short segments
#   gpu_status                                     rewrites and prints <queue dir>/status.md: the holder and the queue in order
#   gpu_serve [PORT]                               serves a live page of the queue and its recent runs on http://127.0.0.1:PORT/ (8008)
#   CPU_SLOTS=<k> CPU_JOB='<one line>' cpu_run LOG CMD [ARG...]
#                                                  waits until k slots of 1 GB are free (1 by default), holds them for all of CMD,
#                                                  exits with its status; a node --test CMD passes --test-concurrency=<at most k>
#   cpu_status                                     prints the runs holding CPU slots and those waiting
# or, without sourcing, zsh scripts/gpuQueue.zsh run LOG CMD [ARG...] | status | serve [PORT]
# and node scripts/cpuQueue.mjs run LOG CMD [ARG...] | status.
typeset -g GPU_QUEUE_SCRIPT=${${(%):-%x}:A:h}/gpuQueue.mjs CPU_QUEUE_SCRIPT=${${(%):-%x}:A:h}/cpuQueue.mjs

_queue_node() {
  local node=${GPU_QUEUE_NODE:-$HOME/.nvm/versions/node/v22.22.2/bin/node}
  [[ -x $node ]] || node=node
  $node "$@"
}
_gpu_queue() { _queue_node $GPU_QUEUE_SCRIPT "$@" }
gpu_run() { _gpu_queue run "$@" }
gpu_status() { _gpu_queue status }
gpu_serve() { _gpu_queue serve "$@" }
cpu_run() { _queue_node $CPU_QUEUE_SCRIPT run "$@" }
cpu_status() { _queue_node $CPU_QUEUE_SCRIPT status }

if [[ $ZSH_EVAL_CONTEXT == toplevel ]]; then
  _gpu_queue "$@"
fi
