# ONT_Autocycler

Full workflow for complete autocycler denovo assembly of bacterial genomes.

## Workflow I/O

### Inputs

```wdl
input {
    File input_reads
    Int num_subsamples = 4
    Int min_read_depth_ss = 25
    Array[String] assemblers = ["raven", "miniasm", "flye", "metamdbg", "necat", "nextdenovo", "plassembler", "canu"]
    String read_type = "ont_r10"
    Int min_depth_abs = 5
    Float min_depth_rel = 0.1
    String? genomesize
    Int? max_contigs
    Int? kmer_size
    Float? cutoff
    Float? min_identity
    Int? max_unitigs
    Float? mad
}
```

### Outputs

```wdl
output {
    File atc_assemblies = FinalizeAssembly.assemblies
    File atc_autocycler_out = FinalizeAssembly.autocycler_out
    File atc_consensus_assembly_fa = FinalizeAssembly.consensus_assembly_fa
    File atc_consensus_assembly_gfa = FinalizeAssembly.consensus_assembly_gfa
    File atc_contig_headers = FinalizeAssembly.contig_headers
    Int atc_num_contigs = FinalizeAssembly.num_contigs
    File atc_final_log = ConsolidateLogs.final_log
}
```

***

## Development

### Dependencies

- [uv](https://astral.sh/uv)
- [make](https://www.gnu.org/software/make/)
- [docker](https://www.docker.com)
- [shellcheck](https://www.shellcheck.net)

### Setup dev environment

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

### Managing Tags and Versions

for versioning this workflow, be sure to change your version in `.VERSION` and use `make tag` to create and push the tag.

Manage versions for your containers in their respective makefiles.

### Running miniwdl check without committing changes
just use `make check` and both miniwdl and yamllint will run on every `*.wdl` and `*.yaml` file in the project.
