version 1.0
import "../structs/Structs.wdl"
import "../tasks/Autocycler.wdl" as ATC

workflow Autocycler {
    meta {
        description: "Pipeline for the execution of autocycler to generate robust denovo assemblies of bacterial genomes. Adapted from Ryan Wick's script at: https://github.com/rrwick/Autocycler/blob/main/pipelines/Automated_Autocycler_Bash_script_by_Ryan_Wick/autocycler_full.sh"
        author: "Michael J. Foster"
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

    call ATC.Subsample {
        input:
            input_reads = input_reads,
            subsample_count = num_subsamples,
            min_read_depth = min_read_depth_ss
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
                    min_depth_rel = min_depth_rel,
                    min_depth_abs = min_depth_abs,
            }
        }
    }
    # gather up everything.
    Array[File] gathered_assemblies = flatten(Assemble.assembler_output)
    Array[File] assembly_logs = flatten(Assemble.log)

    # consolidate our assemblies into a single tarball
    call ATC.ConsolidateAssemblies { input: tarballs = gathered_assemblies }

    # upweight contigs, compress everything into a single graph, then cluster, trim the good clusters, and generate our consensus
    call ATC.FinalizeAssembly { input: assemblies_in = ConsolidateAssemblies.assemblies }

    # first we make an Array[File] of every log we've generated so far. by adding the newest logs to our gathered assembly logs.
    Array[File] all_logs = flatten([ [Subsample.log], assembly_logs, [FinalizeAssembly.log] ])
    # now we consolidate our logs into a single final log
    call ATC.ConsolidateLogs { input: logs = all_logs }

    output {
        File assemblies = FinalizeAssembly.assemblies
        File autocycler_out = FinalizeAssembly.autocycler_out
        File consensus_assembly_fa = FinalizeAssembly.consensus_assembly_fa
        File consensus_assembly_gfa = FinalizeAssembly.consensus_assembly_gfa
        File contig_headers = FinalizeAssembly.contig_headers
        Int num_contigs = FinalizeAssembly.num_contigs
        File final_log = ConsolidateLogs.final_log
    }
}