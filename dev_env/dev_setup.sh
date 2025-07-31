#!/usr/bin/env bash
set -euo pipefail

# howdy. this script is gonna initialize the dev environment so that everything works properly.

# first things first, we should set up the env and install our linters and whatnot.
uv venv
source .venv/bin/activate
uv pip install -r dev_deps.txt

# now set it up as a git repo.
git init

# now install the pre-commit hooks
pre-commit install

# echo the command to activate the venv and exit cleanly :)

echo "Dev environment has been setup! execute the following command to begin:"
echo "source .venv/bin/activate"
exit 0
