version 1.0
import "../structs/Structs.wdl"

task Subsample {
    input {
        File input_reads
        Int subsample_count
        Int min_read_depth
        String? genomesize


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
        genomesize=~{genomesize}
        if [[ -n "$genomesize" ]]; then
            genome_size="$genomesize"
            echo "$genome_size" > est_genome_size.txt
        else
            # now let's estimate our genome size using autocycler's helper function
            autocycler helper genome_size \
                --threads "$NPROCS" \
                --reads ~{input_reads} > est_genome_size.txt
        fi
            # cool, we need this for everything else in this workflow. dump to a file.
            genome_size="$(cat est_genome_size.txt)"

        # ok now we can get to subsampling!
        autocycler subsample \
            --reads ~{input_reads} \
            --out_dir subsamples \
            --genome_size "$genome_size" \
            --min_read_depth ~{min_read_depth} \
            --count ~{subsample_count} 2>>autocycler.stderr
    >>>

    output {
        Int genome_size = read_string("est_genome_size.txt")
        Array[File] subsamples = glob("subsamples/sample_*.fastq")
        File log = "autocycler.stderr"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             96,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
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
        String genome_size
        String read_type
        Float min_depth_rel
        Int min_depth_abs

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 365 + 3 * ceil(size(reads, "GB"))

    command <<<
        set -euo pipefail
        shopt -s nullglob

        NPROCS=$(cat /proc/cpuinfo | awk '/^processor/{print}' | wc -l)

        # make our output directory
        mkdir -p assemblies

        sid="$(basename ~{reads} | sed -E 's/.*_([0-9]+)\.fastq/\1/')"
        echo "Using ~{assembler} on subsample ${sid}."
        echo "Using reads: $(basename ~{reads})"
        echo "Using genome_size: ~{genome_size}"

        if [[ "~{assembler}" == "plassembler" ]]; then
            echo "Plassembler execution detected. Setting plassembler specific arguments"
            ARGS="--args --chromosome 900000 "
        else
            ARGS=""
        fi
        # now let's assemble our reads using autocycler helper.
        echo "Beginning assembly using Autocycler helper."
        (
            autocycler helper ~{assembler} \
                --reads ~{reads} \
                --out_prefix "assemblies/~{assembler}_${sid}" \
                --threads "$NPROCS" \
                --genome_size ~{genome_size} \
                --read_type ~{read_type} \
                --min_depth_rel ~{min_depth_rel} \
                --min_depth_abs ~{min_depth_abs} $ARGS
        ) 2>>autocycler.stderr # this will probably fail.

        # now to pack everything up
        tar -czvf "assemblies.tar.gz" assemblies
    >>>

    output {
        File assembler_output = "assemblies.tar.gz"
        File log = "autocycler.stderr"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             64,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "mjfos2r/autocycler:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
        timeoutMinutes: 360 # timeout after 6 hours... Shouldn't need this long but just in case.
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
        cpu_cores:          16,
        mem_gb:             64,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
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

task FinalizeAssembly {
    input {
        File assemblies_in

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        assemblies_in: "tarball of assemblies to finalize into a consensus assembly."
    }

    Int disk_size = 365 + 2 * ceil(size(assemblies_in, "GB"))

    command <<<
        set -euo pipefail
        shopt -s nullglob

        # unpack the tarball
        tar -xzvf "~{assemblies_in}"

        # upweight plassembler circular contigs since it's adept at pulling them together
        for asm in assemblies/plassembler*.fasta; do
            echo "Upweighting plassembler contigs"
            sed -i 's/circular=True/circular=True Autocycler_cluster_weight=2/' "$asm"
        done

        # give contigs from canu and flye assemblies extra weight too since these are good assemblers
        for asm in assemblies/canu*.fasta assemblies/flye*.fasta; do
            echo "Upweighting contigs from Canu and Flye"
            sed -i 's/^>.*$/& Autocycler_consensus_weight=2/' "$asm"
        done
        echo "Compressing assemblies into a single graph."
        # Now we compress em into a single graph!
        (autocycler compress -i assemblies -a autocycler_out) 2>>autocycler.stderr

        # cluster em!
        echo "Generating Clusters using autocycler cluster."
        (autocycler cluster -a autocycler_out) 2>>autocycler.stderr

        # trim each good cluster
        for c in autocycler_out/clustering/qc_pass/cluster_*; do
            echo "Trimming cluster $(basename ${c})"
            (autocycler trim -c "$c") 2>>autocycler.stderr
            echo "Trimming resolving cluster $(basename ${c})"
            (autocycler resolve -c "$c") 2>>autocycler.stderr
        done

        # combine everything into a consensus assembly!
        echo "Combining all QC passed clusters into a single consensus."
        (autocycler combine -a autocycler_out -i autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa) 2>>autocycler.stderr

        # lets also count our contigs and dump their names
        echo "getting contig count"
        cat autocycler_out/consensus_assembly.fasta | grep '>' | tee contig_headers.txt | wc -l >contig_count.txt

        # Now that we've got everything finished. Let's pack up the autocycler_out directory into a tarball
        # and also provide the final assembly.fa and assembly.gfa as direct outputs.
        echo "compressing autocycler_out directory"
        tar -zcvf autocycler_out.tar.gz autocycler_out/
        # and lets pack up the assemblies directory too since we'll need them for debugging.
        echo "compressing assemblies directory"
        tar -zcvf assemblies.tar.gz assemblies/
    >>>

    output {
        File consensus_assembly_fa = "autocycler_out/consensus_assembly.fasta"
        File consensus_assembly_gfa = "autocycler_out/consensus_assembly.gfa"
        File autocycler_out = "autocycler_out.tar.gz"
        File assemblies = "assemblies.tar.gz"
        File log = "autocycler.stderr"
        File contig_headers = "contig_headers.txt"
        Int num_contigs = read_int("contig_count.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             64,
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
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task ConsolidateLogs {
    input {
        Array[File] logs
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 365 + 2 * ceil(size(logs, "GB"))

    command <<<
        set -euo pipefail
        shopt -s nullglob

        cat ~{sep=' ' logs} > autocycler_full.stderr
    >>>

    output {
        File final_log = "autocycler_full.stderr"
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