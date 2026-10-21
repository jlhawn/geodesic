# The GPU queue of scripts/gpuQueue.mjs for zsh:
#   source scripts/gpuQueue.zsh
#   gpu_run LOG CMD [ARG...]   waits its turn, holds the GPU for all of CMD, exits with its status
#   gpu_status                 the holder, its command and age, and the queue in order
# or, without sourcing, zsh scripts/gpuQueue.zsh run LOG CMD [ARG...] | status.
typeset -g GPU_QUEUE_SCRIPT=${${(%):-%x}:A:h}/gpuQueue.mjs

_gpu_queue() {
  local node=${GPU_QUEUE_NODE:-$HOME/.nvm/versions/node/v22.22.2/bin/node}
  [[ -x $node ]] || node=node
  $node $GPU_QUEUE_SCRIPT "$@"
}
gpu_run() { _gpu_queue run "$@" }
gpu_status() { _gpu_queue status }

if [[ $ZSH_EVAL_CONTEXT == toplevel ]]; then
  _gpu_queue "$@"
fi
