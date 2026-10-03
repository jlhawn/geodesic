# Sourced by scripts/verdaPush.sh and scripts/verdaPull.sh (with REPO set):
# finds the instance and reaches it. HOST, when set, is used as it is:
# user@address for ssh, or 'local' for this machine, the remote paths then
# being local ones. Otherwise each resolve asks `verda --agent vm list -o
# json` for the running instance called NAME and its address, reached as
# SSH_USER (root) with a known-hosts file per instance id
# ($HOME/.ssh/verda_known_hosts_<id>, its host key accepted on first
# contact), since an instance recreated on the kept OS volume gets new host
# keys and may get another address. VERDA (verda) and SSH (ssh) name the
# commands.
VERDA=${VERDA:-verda} SSH=${SSH:-ssh} SSH_USER=${SSH_USER:-root}
target= rsh= why=

resolve() {
  local list id status ip
  why=
  rsh="$SSH -o BatchMode=yes -o ConnectTimeout=20 -o ServerAliveInterval=30 -o ServerAliveCountMax=4"
  if [ -n "$HOST" ]; then target=$HOST; return 0; fi
  list=$($VERDA --agent vm list -o json 2>&1) || { why="verda vm list failed: $(echo "$list" | tr '\n' ' ' | cut -c1-200)"; return 1; }
  read -r id status ip <<< "$(echo "$list" | node "$REPO/scripts/verdaInstances.mjs" address "$NAME")"
  if [ "$status" != running ] || [ -z "$ip" ] || [ "$ip" = - ]; then why="$NAME is ${status:-absent}${id:+ ($id)}"; return 1; fi
  target=$SSH_USER@$ip
  mkdir -p "$HOME/.ssh"
  rsh="$rsh -o StrictHostKeyChecking=accept-new -o UserKnownHostsFile=$HOME/.ssh/verda_known_hosts_$id"
}

remote() { if [ "$target" = local ]; then bash -c "$1"; else $rsh "$target" "$1"; fi; }
at() { if [ "$target" = local ]; then echo "$1"; else echo "$target:$1"; fi; }
