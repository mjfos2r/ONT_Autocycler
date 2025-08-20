version 1.0
import "../tasks/Autocycler.wdl" as ATC
import "../tasks/Dnaapler.wdl" as DNAP

workflow Autocycler {
    meta {
        description: "Pipeline for the execution of autocycler to generate robust denovo assemblies of bacterial genomes. Adapted from Ryan Wick's script at: https://github.com/rrwick/Autocycler/blob/main/pipelines/Automated_Autocycler_Bash_script_by_Ryan_Wick/autocycler_full.sh"
        author: "Michael J. Foster"
    }

    parameter_meta {
        sample_id: "sample id for this assembly."
        input_reads: "reads input as fastq."
        # subsample
        num_subsamples: "how many subsamples to generate from our reads? Default: 4"
        min_read_depth_ss: "minimum read depth that will be allowed for subsampling. Default: 25"
        genome_size: "Optional string specifying the estimated length of the genome in bp. Pass this variable if subsample crashes due to OOM."
        # assemble
        assemblers: "which assemblers will we be using for this sample? default: [raven miniasm flye metamdbg necat nextdenovo plassembler canu]"
        min_depth_abs: "exclude contigs with read depth less than this absolute value. [Default: 5]"
        min_depth_rel: "exclude contigs with read depth less than this fraction of the longest contig's depth. [Default: 0.1]"
        read_type: "Type of read to be assembled. [Default: 'ont_r10'] [Options: ont_r9, ont_r10, pacbio_clr, pacbio_hifi]"
        # Finalization: compress, cluster, trim, resolve, combine
        max_contigs: "integer specifying the maximum number of contigs allowed per assembly. [Default: 25]"
        kmer_size: "integer specifying the kmer size for De Bruijn graph construction. [Default: 51]"
        cutoff: "float specifying the cutoff distance threshold for hiearchical clustering. [Default: 0.2]"
        min_identity: "float specifying the minimum alignment identity for trimming alignments [Default: 0.75]"
        max_unitigs: "integer specifying the maximum number of unitigs used for overlap alignment, set to 0 to disable trimming. [Default: 5.0]"
        mad: "float specifying the allowed variability in cluster length, measured in Median Absolute Deviations. Set to 0 to disable exclusion of length outliers. [Default: 5.0]"
    }

    input {
        String sample_id
        File input_reads
        Int num_subsamples = 4
        Int min_read_depth_ss = 25
        Array[String] assemblers = ["raven", "miniasm", "flye", "metamdbg", "necat", "nextdenovo", "plassembler", "canu"]
        String read_type = "ont_r10"
        Int min_depth_abs = 5
        Float min_depth_rel = 0.1
        String? genomesize
        Int max_contigs
        Int kmer_size
        Float cutoff
        Float min_identity
        Int max_unitigs
        Float mad
    }

    call ATC.Subsample {
        input:
            input_reads = input_reads,
            subsample_count = num_subsamples,
            min_read_depth = min_read_depth_ss,
            genomesize = genomesize
    }
    # okay we now have an Array[File] containing our subsampled reads. Let's get to scattergatherin.
    scatter (subsample in Subsample.subsamples) {
        scatter (tool in assemblers) {
            call ATC.Assemble {
                input:
                    reads = subsample,
                    assembler = tool,
                    genome_size = Subsample.genome_size,
                    read_type = read_type,
                    min_depth_abs = min_depth_abs,
                    min_depth_rel = min_depth_rel,
            }
        }
    }
    # gather up everything.
    Array[File] gathered_assemblies = flatten(Assemble.assembler_output)
    Array[File] assembly_logs = flatten(Assemble.log)

    # consolidate our assemblies into a single tarball
    call ATC.ConsolidateAssemblies { input: tarballs = gathered_assemblies }

    # upweight contigs, compress everything into a single graph, then cluster, trim the good clusters, and generate our consensus
    call ATC.FinalizeAssembly {
        input:
            assemblies_in = ConsolidateAssemblies.assemblies,
            max_contigs = max_contigs,
            kmer_size = kmer_size,
            cutoff = cutoff,
            min_identity = min_identity,
            max_unitigs = max_unitigs,
            mad = mad,
    }
    call DNAP.Dnaapler {
        input:
            draft_assembly = FinalizeAssembly.consensus_assembly_fa,
            sample_id = sample_id
    }

    # first we make an Array[File] of every log we've generated so far. by adding the newest logs to our gathered assembly logs.
    Array[File] all_logs = flatten([ [Subsample.log], assembly_logs, [FinalizeAssembly.log] ])
    # now we consolidate our logs into a single final log
    call ATC.ConsolidateLogs { input: logs = all_logs }

    output {
        File atc_assemblies = FinalizeAssembly.assemblies
        File atc_autocycler_out = FinalizeAssembly.autocycler_out
        File atc_consensus_assembly_fa = FinalizeAssembly.consensus_assembly_fa
        File atc_consensus_assembly_fa_reoriented = Dnaapler.assembly_reoriented
        File atc_consensus_assembly_gfa = FinalizeAssembly.consensus_assembly_gfa
        File atc_contig_headers = FinalizeAssembly.contig_headers
        Int atc_num_contigs = FinalizeAssembly.num_contigs
        File atc_final_log = ConsolidateLogs.final_log
    }
}