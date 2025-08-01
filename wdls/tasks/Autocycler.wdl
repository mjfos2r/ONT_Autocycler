version 1.0
import "../structs/Structs.wdl"

workflow Autocycler {

    meta {
        description: "Description of the workflow"
    }
    parameter_meta {
        input_reads: "reads input as fastq"
        num_subsamples: "how many subsamples to generate from our reads? Default: 4"
        min_read_depth_ss: "minimum read depth that will be allowed for subsampling. Default: 25"
        assemblers: "which assemblers will we be using for this sample? default: [raven miniasm flye metamdbg necat nextdenovo plassembler canu]"
        read_type: "available choices: ont_r9 ont_r10 pacbio_clr pacbio_hifi, default: ont_r10"
        min_depth_abs: "exclude contigs with read depth less than this absolute value. Default: 5"
        min_depth_rel: "exclude contigs with read depth less than this fraction of the longest contig's depth. Default: 0.1"
        prefix: "prefix to output files. Testing to see if dockstore is busted for me..."
    }

    input {
        File input_reads
        Int num_subsamples = 4
        Int min_read_depth_ss = 25
        Array[String] assemblers = ["raven", "miniasm", "flye", "metamdbg", "necat", "nextdenovo", "plassembler", "canu"]
        String read_type = "ont_r10"
        Int min_depth_abs = 5
        Float min_depth_rel = 0.1
        String prefix
    }

    call Subsample {
        input:
            input_reads = input_reads,
            subsample_count = num_subsamples,
            min_read_depth = min_read_depth_ss
    }

    # okay we now have an Array[File] containing our subsampled reads. Let's get to scattergatherin.

    scatter (subsample in Subsample.subsamples) {
        scatter (tool in assemblers) {
            call Assemble {
                input:
                    reads = subsample,
                    assembler = tool,
                    genome_size = Subsample.genome_size,
                    read_type = read_type,
                    min_depth_rel = min_depth_rel,
                    min_depth_abs = min_depth_abs,
                    #args = assembler_args[tool] # this isn't legal in WDL1.0 but it should be.
            }
        }
    }
    Array[File] gathered_assemblies = flatten(Assemble.assembler_output)

    call ConsolidateAssemblies {
        input:
            tarballs = gathered_assemblies,
    }


    output {
        File assemblies = ConsolidateAssemblies.assemblies
    }
}

task Subsample {
    input {
        File input_reads
        Int subsample_count
        Int min_read_depth


        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        input_reads: "input reads to be subsampled"
        subsample_count: "how many subsamples to create"
        min_read_depth: "minimum allowed read depth."
    }

    Int disk_size = 365 + 3 * ceil(size(input_reads, "GB"))

    command <<<
        set -euo pipefail
        shopt -s nullglob
        NPROCS=$(cat /proc/cpuinfo | awk '/^processor/{print}' | wc -l)

        # make our output directory
        mkdir -p subsamples

        # now let's estimate our genome size using autocycler's helper function
        autocycler helper genome_size \
            --threads "$NPROCS" \
            --reads ~{input_reads} > est_genome_size.txt
        # cool, we need this for everything else in this workflow. dump to a file.
        genome_size="$(cat est_genome_size.txt)"

        # ok now we can get to subsampling!
        autocycler subsample \
            --reads ~{input_reads} \
            --out_dir subsamples \
            --genome_size "$genome_size" \
            --min_read_depth ~{min_read_depth} \
            --count ~{subsample_count} 2>> autocycler_subsample.stderr
    >>>

    output {
        Int genome_size = read_int("est_genome_size.txt")
        Array[File] subsamples = glob("subsamples/sample_*.fastq")
        File log = "autocycler_subsample.stderr"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "mjfos2r/autocycler:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task Assemble {
    input {
        File reads
        String assembler
        Int genome_size
        String read_type
        Float min_depth_rel
        Int min_depth_abs

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 100 + 3 * ceil(size(reads, "GB"))

    command <<<
        set -euo pipefail
        shopt -s nullglob
        NPROCS=$(cat /proc/cpuinfo | awk '/^processor/{print}' | wc -l)

        # make our output directory
        mkdir -p assemblies

        sid="$(basename ~{reads} | sed -E 's/.*_([0-9]+).fastq/\1/')"

        if [[ "~{assembler}" == "plassembler" ]]; then
            ARGS="--args --chromosome 900000"
        else
            ARGS=""
        fi

        # now let's estimate our genome size using autocycler's helper function
        autocycler helper ~{assembler} \
            --reads ~{reads} \
            --out_prefix "assemblies/~{assembler}_${sid}" \
            --threads "$NPROCS" \
            --genome_size ~{genome_size} \
            --read_type ~{read_type} \
            --min_depth_rel ~{min_depth_rel} \
            --min_depth_abs ~{min_depth_abs} "$ARGS" # this will probably fail.

        mkdir tarball
        # now to pack everything up
        tar -czvf "tarball/~{assembler}_${sid}.tar.gz" assemblies
    >>>

    output {
        File assembler_output = glob("tarball/*.tar.gz")[0]
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "mjfos2r/autocycler:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
task ConsolidateAssemblies {
    input {
        Array[File] tarballs
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 365 + 3 * ceil(size(tarballs, "GB"))

    command <<<
        set -euo pipefail
        shopt -s nullglob
        NPROCS=$(cat /proc/cpuinfo | awk '/^processor/{print}' | wc -l)

        mkdir -p assemblies
        for ball in ~{sep=' 'tarballs}; do
            tar -xzvf "$ball" -C assemblies --strip-components=1
        done

        tar -zcvf assemblies.tar.gz assemblies/
    >>>

    output {
        File assemblies = "assemblies.tar.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "mjfos2r/autocycler:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}