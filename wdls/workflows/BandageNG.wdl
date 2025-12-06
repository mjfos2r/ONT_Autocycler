version 1.0
import "../tasks/BandageNG.wdl" as BNG

workflow BandageNG {
    meta {
        description: "Generate static image from an assembly graph in GFA format."
    }

    parameter_meta {
        gfa: "assembly graph."
        assembly_id: "assembly id for output file."
    }

    input {
        File gfa
        String assembly_id
    }

    call BNG.DrawGFA {
        input:
            gfa = gfa,
            assembly_id = assembly_id
    }

    output {
        File gfa_plot = DrawGFA.image
    }
}