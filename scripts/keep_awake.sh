#!/usr/bin/env bash
# Keep Mac awake while Geant4 campaigns run (safe to close the lid on AC power).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PIDFILE="$ROOT/data/.keep_awake.pid"
LOG="$ROOT/data/keep_awake.log"

stop() {
  if [[ -f "$PIDFILE" ]]; then
    local pid
    pid="$(cat "$PIDFILE")"
    if kill -0 "$pid" 2>/dev/null; then
      kill "$pid" 2>/dev/null || true
      echo "Stopped keep-awake (pid $pid)"
    fi
    rm -f "$PIDFILE"
  fi
}

start() {
  stop
  mkdir -p "$ROOT/data"

  # AC power: never sleep or spin down disks while plugged in.
  pmset -c sleep 0 disksleep 0 standby 0 2>/dev/null || true

  # Prevent idle/display/system sleep. -w waits on simulation PIDs if present.
  local wait_args=()
  for pat in './Main ' './MainHad '; do
    while IFS= read -r pid; do
      wait_args+=(-w "$pid")
    done < <(pgrep -f "$pat" 2>/dev/null || true)
  done

  {
    echo "=== keep_awake started $(date) ==="
    echo "pmset AC: $(pmset -g custom 2>/dev/null | sed -n '/AC Power:/,/Battery/p' | head -15)"
    echo "waiting on: ${wait_args[*]:-none (indefinite)}"
  } >>"$LOG"

  if ((${#wait_args[@]})); then
    nohup caffeinate -dims "${wait_args[@]}" >>"$LOG" 2>&1 &
  else
    nohup caffeinate -dims >>"$LOG" 2>&1 &
  fi
  echo $! >"$PIDFILE"
  echo "Keep-awake running (pid $(cat "$PIDFILE")). Log: $LOG"
  echo "Plug into AC power before closing the lid."
  echo "If the Mac still sleeps on lid close, run once: sudo pmset -a disablesleep 1"
  echo "Restore later with: sudo pmset -a disablesleep 0"
}

status() {
  if [[ -f "$PIDFILE" ]] && kill -0 "$(cat "$PIDFILE")" 2>/dev/null; then
    echo "keep_awake: running (pid $(cat "$PIDFILE"))"
  else
    echo "keep_awake: not running"
  fi
  pmset -g custom 2>/dev/null | sed -n '/AC Power:/,/Battery/p' | head -10
}

case "${1:-start}" in
  start) start ;;
  stop) stop ;;
  status) status ;;
  restart) stop; start ;;
  *) echo "Usage: $0 {start|stop|status|restart}"; exit 1 ;;
esac
