# TemplateTerraWorkflow

This repository contains a template workflow that serves as a comprehensive starting point for writing pipelines to run in Terra.
To begin using this, simply click "Use this template" at the top right of the page or alternatively use the following command if
you have [gh cli](https://cli.github.com/manual/):
```
gh repo create my-terra-workflow --template mjfos2r/TemplateTerraWorkflow
```

## Dependencies

- [uv](https://astral.sh/uv)
- [make](https://www.gnu.org/software/make/)
- [docker](https://www.docker.com)
- [shellcheck](https://www.shellcheck.net)

## Setup dev environment

To get this all set up and groovy, simply ensure you have the dependencies listed above and execute `./dev_env/dev_setup.sh`

Which does the following:

1. creates a virtual environment using uv
2. installs dev_deps.txt using uv.
   contents:
    - [pre-commit](https://pre-commit.com) to check files using the following linters at time before commit.
    - [miniwdl](https://github.com/chanzuckerberg/miniwdl) for linting and checking wdl files. (needs shellcheck)
    - [yamllint](https://github.com/adrienverge/yamllint) for linting and checking yaml files.
3. installs pre-commit to this repository.

Activate the virtual environment and you're good to go.
>protip: I have these handy aliases in my .bashrc:
>
> - `alias uva="source .venv/bin/activate"`
> - `alias uvd="deactivate"`
> - `alias uvi="uv pip install"`

## Develop inside of Docker container

I am also writing a docker container to serve as the development environment with all dependencies pre-installed.

To use it, you can execute the helper script: `./dev_env/launch_container.sh`

## Important placeholder locations
Since this was originally *my* template, I have *my* information as placeholder.
To use this template, you will need to change these to your own.
>In addition to everything else that needs changing. this is a template after all.

1. `Makefile`; Line 3; Line 4;
2. `.dockstore.yml`; Line 3; Line 5;
3. `wdls/tasks/Tasks.wdl`; Line 33;
4. `containers/**/Makefile`; Line 3; Line 4;

## Managing containers
I've included two example Dockerfiles and their makefiles.

This repo is set up to manage multiple containers required for the workflow you're writing.
If you don't need to write any containers, you can safely delete the containers directory and ignore the commands in the makefile.
If you *do* need to write containers, I've simplified the whole process.

1. within `containers/`, create separate directories for each that you write.
2. copy over the makefile from one of the template directories and change the metadata at the top.
3. write the Dockerfiles
4. from the project root (this dir), execute the following command: `make build-containers`

You can also do `make build` and `make push` in each of your container directories while testing so you don't have to wait for everything to build and push.

## Managing Tags and Versions

for versioning this workflow, be sure to change your version in `.VERSION` and use `make tag` to create and push the tag.

Manage versions for your containers in their respective makefiles.

## Running miniwdl check without committing changes
just use `make check` and both miniwdl and yamllint will run on every `*.wdl` and `*.yaml` file in the project.
