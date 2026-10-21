#!/bin/sh
# Runs on the host before the dev container starts (initializeCommand).
# Writes nothing into the repo: all outputs go to $dir, which is mounted at /ssh-agent-dir.
set -u
name=$(basename "$PWD")
dir="$XDG_RUNTIME_DIR/devcontainer-agent-$name"
sock="$dir/agent.sock"
mkdir -p "$dir" && chmod 700 "$dir"

# askpass program for signing confirmations
askpass=${SSH_ASKPASS:-}
if [ -z "$askpass" ]; then
    for c in lxqt-openssh-askpass ksshaskpass; do
        if p=$(command -v "$c"); then askpass=$p; break; fi
    done
fi
[ -n "$askpass" ] || echo "warning: no askpass program found; signing confirmation will fail" >&2

# Git identity and signing intent from the host config
f="$dir/gitidentity"
rm -f "$f"
for k in user.name user.email gpg.format commit.gpgsign; do
    v=$(git config --get "$k") && git config --file "$f" "$k" "$v"
done

# Signing key (ssh format, file path only)
key=""
if [ "$(git config --get gpg.format)" = ssh ] && key=$(git config --get user.signingkey); then
    case "$key" in
        key::*) echo "host signing key is a literal; not loading it" >&2; key="" ;;
        "~/"*)  key="$HOME/${key#\~/}" ;;
    esac
    key=${key%.pub}
fi
if [ -n "$key" ]; then
    if [ -f "$key.pub" ]; then
        git config --file "$f" user.signingkey "key::$(cat "$key.pub")"
    else
        echo "warning: $key.pub not found; not loading signing key" >&2
        key=""
    fi
fi

# Read-only deploy key for this repo, if one exists
deploy="$HOME/.ssh/deploy_$name"
if [ -f "$deploy" ]; then
    [ -f "$deploy.pub" ] && cp "$deploy.pub" "$dir/deploy.pub"
else
    echo "note: no deploy key at $deploy; fetching from inside the container won't work (pull on the host instead)" >&2
    deploy=""
fi

[ -n "$key$deploy" ] || exit 0

# Reuse the agent if alive; otherwise start one as a transient user unit
rc=0; SSH_AUTH_SOCK="$sock" ssh-add -l >/dev/null 2>&1 || rc=$?
if [ "$rc" -eq 2 ]; then
    rm -f "$sock"
    systemd-run --user --collect --unit="devcontainer-agent-$(systemd-escape "$name")" \
        --setenv=SSH_ASKPASS="$askpass" \
        --setenv=SSH_ASKPASS_REQUIRE=force \
        --setenv=WAYLAND_DISPLAY="${WAYLAND_DISPLAY:-}" \
        --setenv=DISPLAY="${DISPLAY:-}" \
        ssh-agent -D -a "$sock"
    i=0; while [ ! -S "$sock" ] && [ $i -lt 50 ]; do sleep 0.1; i=$((i+1)); done
fi

# Load keys into an empty or new agent: signing key with confirmation, deploy key without
if [ "$rc" -ne 0 ]; then
    [ -n "$key" ] && SSH_AUTH_SOCK="$sock" SSH_ASKPASS="$askpass" SSH_ASKPASS_REQUIRE=force \
        ssh-add -c "$key"
    [ -n "$deploy" ] && SSH_AUTH_SOCK="$sock" ssh-add "$deploy"
fi
exit 0