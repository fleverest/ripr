#!/bin/sh
set -e

# Load git user.name and user.email from host
git config --global include.path "/ssh-agent-dir/gitidentity"

# Prevent git push in devcontainer
git config --global --unset-all url."no-push-from-container://".pushInsteadOf || true
git config --global --add url."no-push-from-container://".pushInsteadOf "git@github.com:"
git config --global --add url."no-push-from-container://".pushInsteadOf "https://github.com/"
