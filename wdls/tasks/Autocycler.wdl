version 1.0
import "../structs/Structs.wdl"

task Subsample {
    meta {
        description: "Task to generate sets of subsampled reads using 'autocycler subsample'."
        author: "Michael J. Foster"
    }
    parameter_meta {
        input_reads: "input reads to be subsampled"
        subsample_count: "how many subsamples to create"
        min_read_depth: "minimum allowed read depth."
        genomesize: "Optional string detailing the genomesize used for assembly and subsampling. [Example:'1500000'] [Default: calculated using autocycler helper]"
    }
    input {
        File input_reads
        Int subsample_count
        Int min_read_depth
        String? genomesize
        Int max_fs_gb = 4
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 365 + 3 * ceil(size(input_reads, "GB"))

    command <<<
        set -euo pipefail
        shopt -s nullglob
        NPROCS=$(cat /proc/cpuinfo | awk '/^processor/{print}' | wc -l)
        echo "Using ${NPROCS} for subsampling."
        # set up handling to prevent OOM for large read sets. For some reason the
        # documented estimation for memory usage in raven is incorrect.
        # it is using far more than 1.2 * size + 16GB.
        # So we're just gonna estimate genome_size based on a subsample of half of our reads
        # and go from there.

        reads_fs_bytes=$(stat -c%s ~{input_reads})
        max_fs_gb=~{max_fs_gb} # we're gonna try with 8. It requires roughly 9x input_fs GB in ram.
        max_fs_bytes=$((max_fs_gb * 1024 * 1024 * 1024))

        # make our output directory
        mkdir -p subsamples
        genomesize=~{genomesize}
        if [[ -n "$genomesize" ]]; then
            echo "$genomesize" > est_genome_size.txt
            echo "Using provided genome size for subsampling."
        else
            if (( reads_fs_bytes > max_fs_bytes )); then
                echo "input is too large. determining proportion to subsample by."
                prop=$(awk -v s="$reads_fs_bytes" -v t="$max_fs_bytes" 'BEGIN { print t/s }')
                echo "$prop"
                # clamp between 0.1 and 0.9
                clamped=$(awk -v p="$prop" 'BEGIN { if (p < 0.1) p = 0.1; if (p > 0.9) p = 0.9; print p}')
                echo "$clamped"
                echo "subsampling input reads by ${clamped} and using to determine genome_size, please stand by..."
                seqkit sample --threads 14 -p "$clamped" -s 318 ~{input_reads} -o subsampled_input.fastq
                # now let's estimate our genome size using autocycler's helper function
                echo "Determining genome_size of subsampled input using autocycler helper genome_size."
                (
                    autocycler helper genome_size \
                        --threads "$((NPROCS - 2))" \
                        --reads subsampled_input.fastq > est_genome_size.txt
                ) 2>>autocycler.stderr
            else
                # now let's estimate our genome size using autocycler's helper function
                echo "Determining genome_size using autocycler helper genome_size."
                (
                    autocycler helper genome_size \
                        --threads "$((NPROCS - 2))" \
                        --reads ~{input_reads} > est_genome_size.txt
                ) 2>>autocycler.stderr
            fi
        fi

        # cool, we need this for everything else in this workflow. dump to a file.
        genome_size="$(cat est_genome_size.txt)"
        echo "Now generating subsamples using autocycler subsample."
        echo "Number of Subsamples: ~{subsample_count}"
        echo "Minimum read depth:   ~{min_read_depth}"
        echo "genome_size:          $genome_size"
        # ok now we can get to subsampling!
        (
            autocycler subsample \
                --reads ~{input_reads} \
                --out_dir subsamples \
                --genome_size "$genome_size" \
                --min_read_depth ~{min_read_depth} \
                --count ~{subsample_count}
        ) 2>>autocycler.stderr
        echo "Autocycler subsample finished!"
    >>>

    output {
        String genome_size = read_string("est_genome_size.txt")
        Array[File] subsamples = glob("subsamples/sample_*.fastq")
        File log = "autocycler.stderr"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             96,
        disk_gb:            disk_size,
        boot_disk_gb:       50,
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
    }
}

task Assemble {
    meta {
        description: "Task that executes 'autocycler helper {assembler}' on an input set of subsampled reads using the provided assembler."
        author: "Michael J. Foster"
    }
    parameter_meta {
        reads: "File containing the set of subsampled reads to be assembled."
        assembler: "String containing the name of the assembler to use."
        genome_size: "String containing the size of the genome."
        read_type: "Type of read to be assembled.  [Options: ont_r9, ont_r10, pacbio_clr, pacbio_hifi]"
        min_depth_abs: "exclude contigs with read depth less than this absolute value. [Default: 5]"
        min_depth_rel: "exclude contigs with read depth less than this fraction of the longest contig's depth. [Default: 0.1]"
    }
    input {
        File reads
        String assembler
        Int genome_size
        String read_type
        Int min_depth_abs
        Float min_depth_rel
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
        echo "Assembly complete! Compressing assemblies directory."
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
    meta {
        description: "Gather task to consolidate scattered assemblies generated for each assembler, for each subsampled reads set."
        author: "Michael J. Foster"
    }
    parameter_meta {
        tarballs: "Array[File] containing the compressed assemblies.tar.gz file for each assembly generated."
    }
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
    meta {
        description: "Task to finalize this assembly. Performs contig biasing for Flye, Canu, and Plassembler. Compresses all assemblies into a single graph, then generates clusters. Clusters passing QC are trimmed and resolved. Finally, all clusters are combined into a single consensus fasta and graph."
        author: "Michael J. Foster"
    }
    parameter_meta {
        assemblies_in: "tarball of assemblies to finalize into a consensus assembly."
        max_contigs: "integer specifying the maximum number of contigs allowed per assembly. [Default: 25]"
        kmer_size: "integer specifying the kmer size for De Bruijn graph construction. [Default: 51]"
        cutoff: "float specifying the cutoff distance threshold for hiearchical clustering. [Default: 0.2]"
        min_identity: "float specifying the minimum alignment identity for trimming alignments [Default: 0.75]"
        max_unitigs: "integer specifying the maximum number of unitigs used for overlap alignment, set to 0 to disable trimming. [Default: 5.0]"
        mad: "float specifying the allowed variability in cluster length, measured in Median Absolute Deviations. Set to 0 to disable exclusion of length outliers. [Default: 5.0]"
    }
    input {
        File assemblies_in
        Int max_contigs = 25
        Int kmer_size = 51
        Float cutoff = 0.2
        Float min_identity = 0.75
        Int max_unitigs = 5000
        Float mad = 5.0
        #Array[Int]? manual_nodes # list of tree node numbers to manually define clusters.
        ##  usage: ~{sep=',' manual_nodes}
        ##    Could be useful in Bb assembly if I'm going to resolve plasmids after multiassembly before finalizing?
        #Float? min_assemblies # do I need to implement this using some sort of arg construction based on input presence? it defaults to automatic determination.
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 365 + 2 * ceil(size(assemblies_in, "GB"))

    command <<<
        set -euo pipefail
        shopt -s nullglob

        NPROCS=$(cat /proc/cpuinfo | awk '/^processor/{print}' | wc -l)

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
        (autocycler compress -i assemblies -a autocycler_out --max_contigs ~{max_contigs} --kmer ~{kmer_size} --threads "$NPROCS") 2>>autocycler.stderr

        # cluster em!
        echo "Generating Clusters using autocycler cluster."
        (autocycler cluster -a autocycler_out --cutoff ~{cutoff} --max_contigs ~{max_contigs}) 2>>autocycler.stderr

        # trim each good cluster
        for c in autocycler_out/clustering/qc_pass/cluster_*; do
            echo "Trimming cluster $(basename ${c})"
            (autocycler trim -c "$c" --min_identity ~{min_identity} --max_unitigs ~{max_unitigs} --mad ~{mad} --threads "$NPROCS") 2>>autocycler.stderr
            echo "Trimming resolving cluster $(basename ${c})"
            (autocycler resolve -c "$c") 2>>autocycler.stderr
        done

        # combine everything into a consensus assembly!
        echo "Combining all QC passed clusters into a single consensus."
        (autocycler combine -a autocycler_out -i autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa) 2>>autocycler.stderr

        # lets also count our contigs and dump their names
        echo "getting contig count and final assembly length"
        read -r num_contigs asm_length < <(seqkit stats -T autocycler_out/consensus_assembly.fasta | tail -n1 | cut -f4,5)
        echo "$num_contigs" > contig_count.txt
        echo "$asm_length" > asm_length.txt
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
        Int num_contigs = read_int("contig_count.txt")
        Int asm_length = read_int("asm_length.txt")
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
    meta {
        description: "Task to consolidate all autocycler.stderr logs from each task into a single file."
        author: "Michael J. Foster"
    }
    parameter_meta {
        logs: "Array[File] of logfiles gathered from each executed autocycler task."
    }
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