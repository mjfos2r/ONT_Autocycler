version 1.0
import "../tasks/Autocycler.wdl" as ATC
import "../tasks/Dnaapler.wdl" as DNAP
import "../tasks/BandageNG.wdl" as BNG

workflow AutocyclerClean {
    meta {
        description: "Manual curation workflow for cleaning Autocycler assemblies. Use after inspecting consensus_assembly.gfa in Bandage to identify tigs for removal or duplication."
        author: "Michael J. Foster"
    }

    parameter_meta {
        sample_id: "Sample identifier."
        consensus_gfa: "GFA file from Autocycler combine to be cleaned."
        min_depth: "Automatically remove tigs below this depth. Set to 0 to disable. [Default: 0]"
        max_depth: "Automatically remove tigs above this depth. Set to 0 to disable. [Default: 100]"
        min_length: "Minimum length for filtering. All contigs shorter than this are dropped."
        remove_tigs: "Comma-separated list of tig numbers to remove (e.g. '5,7,8,11')."
        duplicate_tigs: "Comma-separated list of tig numbers to duplicate for TIRs (e.g. '4')."
        telos_to_split: "Comma-separated list of tig numbers representing telomere nodes that need to be split (e.g. '4')."
    }

    input {
        String sample_id
        File consensus_gfa
        Int min_depth = 0
        Int max_depth = 100
        Int min_length = 2500
        String? remove_tigs
        String? duplicate_tigs
        String? telos_to_split
    }

    call ATC.CleanGFA {
        input:
            sample_id = sample_id,
            consensus_gfa = consensus_gfa,
            min_depth = min_depth,
            max_depth = max_depth,
            min_length = min_length,
            remove_tigs = remove_tigs,
            duplicate_tigs = duplicate_tigs,
            telos_to_split = telos_to_split
    }

    call DNAP.Dnaapler {
        input:
            draft_assembly = CleanGFA.final_fasta,
            sample_id = sample_id,
    }

    call BNG.DrawGFA {
        input:
            gfa = CleanGFA.final_gfa,
            assembly_id = sample_id + "_manual_clean"
    }

    output {
        File atc_final_gfa = CleanGFA.final_gfa
        File atc_final_gfa_svg = DrawGFA.image
        File atc_final_fasta = CleanGFA.final_fasta
        File atc_final_fasta_reoriented = Dnaapler.assembly_reoriented
        File atc_final_log = CleanGFA.final_log
        Int atc_final_num_contigs = CleanGFA.final_num_contigs
        Int atc_final_asm_length = CleanGFA.final_asm_length
    }
}