#!/bin/sh

# Sign commits via the mounted ssh-agent socket
k=$(ssh-add -L 2>/dev/null) || exit 1
printf 'key::%s\n' "$(printf '%s\n' "$k" | head -n1)"