# The GPU queue of scripts/gpuQueue.mjs for zsh:
#   source scripts/gpuQueue.zsh
#   export GPU_WHO=<session name>                  once per session
#   GPU_JOB='<one line>' gpu_run LOG CMD [ARG...]  waits its turn, holds the GPU for all of CMD, exits with its status
#   gpu_status                                     rewrites and prints <queue dir>/status.md: the holder and the queue in order
#   gpu_serve [PORT]                               serves a live page of the queue and its recent runs on http://127.0.0.1:PORT/ (8008)
# or, without sourcing, zsh scripts/gpuQueue.zsh run LOG CMD [ARG...] | status | serve [PORT].
typeset -g GPU_QUEUE_SCRIPT=${${(%):-%x}:A:h}/gpuQueue.mjs

_gpu_queue() {
  local node=${GPU_QUEUE_NODE:-$HOME/.nvm/versions/node/v22.22.2/bin/node}
  [[ -x $node ]] || node=node
  $node $GPU_QUEUE_SCRIPT "$@"
}
gpu_run() { _gpu_queue run "$@" }
gpu_status() { _gpu_queue status }
gpu_serve() { _gpu_queue serve "$@" }

if [[ $ZSH_EVAL_CONTEXT == toplevel ]]; then
  _gpu_queue "$@"
fi
