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
        min_depth: "Automatically remove tigs at or below this depth. Set to 0 to disable. [Default: 0]"
        remove_tigs: "Comma-separated list of tig numbers to remove (e.g. '5,7,8,11')."
        duplicate_tigs: "Comma-separated list of tig numbers to duplicate for TIRs (e.g. '4')."
    }

    input {
        String sample_id
        File consensus_gfa
        Int min_depth = 0
        String? remove_tigs
        String? duplicate_tigs
    }

    call ATC.CleanAssembly {
        input:
            consensus_gfa = consensus_gfa,
            sample_id = sample_id,
            min_depth = min_depth,
            remove_tigs = remove_tigs,
            duplicate_tigs = duplicate_tigs,
    }

    call DNAP.Dnaapler {
        input:
            draft_assembly = CleanAssembly.cleaned_fasta,
            sample_id = sample_id,
    }

    call BNG.DrawGFA {
        input:
            gfa = CleanAssembly.cleaned_gfa,
            assembly_id = sample_id + "_manual_clean"
    }

    output {
        File atc_cleaned_gfa = CleanAssembly.cleaned_gfa
        File atc_cleaned_fasta = CleanAssembly.cleaned_fasta
        File atc_cleaned_fasta_reoriented = Dnaapler.assembly_reoriented
        File atc_cleaned_gfa_svg = DrawGFA.image
        Int atc_num_contigs_clean = CleanAssembly.num_contigs
        Int atc_asm_length_clean = CleanAssembly.asm_length
        File atc_log_clean = CleanAssembly.log
    }
}