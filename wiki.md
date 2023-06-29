# Local registry

You need to add the key to SSH agent before updating and messing with the local registry
* `eval "$(ssh-agent -s)"` and then `ssh-add ~/.ssh/id_ed25519`

After running these codes run the julia
