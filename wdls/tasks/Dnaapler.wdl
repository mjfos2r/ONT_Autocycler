version 1.0
import "../structs/Structs.wdl"

workflow ReorientContigs {
    meta { description: "Simple workflow to reorient bacterial assemblies using Dnaapler." }
    input {
        File draft_assembly
        String? dnaapler_mode
        String sample_id
    }
    call Dnaapler { input: draft_assembly = draft_assembly, sample_id = sample_id, dnaapler_mode = dnaapler_mode}
    output { File assembly_reoriented = Dnaapler.assembly_reoriented }
}

task Dnaapler {
    input {
        File draft_assembly
        String sample_id
        String dnaapler_mode = "all"
        RuntimeAttr? runtime_attr_override
    }
    parameter_meta {
        draft_assembly: "fasta file containing the draft assembly to reorient"
        sample_id: "sample_id for the assembly being reordered."
        mode: "mode to run dnaapler in [ Default: 'all' ]"
    }
    Int disk_size = 365 + 2 * ceil(size(draft_assembly, "GB"))
    command <<<
        set -euo pipefail
        shopt -s nullglob
        NPROCS=$(cat /proc/cpuinfo | awk '/^processor/{print}' | wc -l)

        #mkdir -p dnaapler_out # dnaapler makes this directory. no need to create it here.

        dnaapler ~{dnaapler_mode} \
        --threads "${NPROCS}" \
        --prefix ~{sample_id} \
        --force \
        --output dnaapler_out \
        --input ~{draft_assembly}
    >>>
    output {
        File assembly_reoriented = "dnaapler_out/~{sample_id}_reoriented.fasta"
    }
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             64,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "mjfos2r/dnaapler:latest"
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