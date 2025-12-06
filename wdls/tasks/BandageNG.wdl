version 1.0
import "../structs/Structs.wdl"

task DrawGFA {
    input {
        File gfa
        String assembly_id
        Int width = 1920
        Int height = 1080
        String textcol = "black"
        String toutcol = "white"
        Int toutline = 5
        Int fontsize = 11
        Int iter = 4
        String ext = ".svg"
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        gfa: "assembly graph."
        assembly_id: "assembly id for output file."
        width: "width = 1920 "
        height: "height = 1080 "
        textcol: "textcol = black "
        toutcol: "toutcol = white "
        toutline: "toutline = 5 "
        fontsize: "fontsize = 11 "
        iter: "iter = 4 "
        ext: "ext = .svg #(CAN ONLY BE .svg, .png, .jpg)"
    }

    Int disk_size = 50 + ceil(size(gfa, "GB"))

    command <<<
        set -euo pipefail
        shopt -s nullglob

        BandageNG image \
            --width ~{width} \
            --height ~{height} \
            --textcol ~{textcol} \
            --toutcol ~{toutcol} \
            --toutline ~{toutline} \
            --fontsize ~{fontsize} \
            --iter ~{iter} \
            --centre \
            --lengths \
            --depth \
            ~{gfa} \
            "~{assembly_id}_assembly_graph~{ext}"
    >>>

    output {
        File image = "~{assembly_id}_assembly_graph_~{ext}"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "mjfos2r/bandageng:latest"
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