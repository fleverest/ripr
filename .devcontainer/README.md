# Getting set up with my Dev Container

## Requirements

1. A Linux host with rootless podman and a systemd user session.
2. Positron with dev containers enabled (`dev.containers.enable`). `${localEnv}` substitution and `initializeCommand` currently need some patches, until upstream fixes land.
3. An `askpass` program (`lxqt-openssh-askpass` or `ksshaskpass`)

## Git

`user.name`, `user.email`, `gpg.format` and `commit.gpgsign` come from your host git config (copied via `initializeCommand`).
Signing requires `gpg.format=ssh` and `user.signingkey` set to a key file, with its `.pub` next to it.
Each signature brings up a confirmation dialog on the host via `askpass`. Deny prompts you didn’t trigger (these could be commits by a coding agent, for example).
The signing key should be registered on GitHub as a signing key only, not an authentication key.

#### Optional: Deploy key

Create an ssh-key at `~/.ssh/deploy_<folder-name>` via
```
ssh-keygen -t ed25519 -f ~/.ssh/deploy_$(basename $PWD) -C "read-only deploy key"
```
Then add it as a deploy key on the github repo settings.
This enables fetching within the container.
Without it, you can only fetch/pull on the host.