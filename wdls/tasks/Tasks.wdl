version 1.0
import "../structs/Structs.wdl"

task FirstTask {
    meta {
        # task-level metadata can go here
        description: "first task"
    }

    parameter_meta {
        # metadata about each input/output parameter can go here
        input_file: "input file for analysis"
    }

    input {
        # task inputs are declared here
        File input_file
        RuntimeAttr? runtime_attr_override
    }

    # other "private" declarations can be made here
    #Int disk_size = 365 + 2 * ceil(size(fastq_files, "GB")) # change this to whatever you need.

    command <<<
        # the command template - this section is required
        # put your shell commands here.
        cp ~{input_file} first_task_output.txt
    >>>

    output {
        # task outputs are declared here
        File output_file = "first_task_output.txt"
    }
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            100, #disk_size,
        boot_disk_gb:       50,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "mjfos2r/mycontainer:latest"
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

task SecondTask {
    meta {
        # task-level metadata can go here
        description: "second task"
    }

    parameter_meta {
        # metadata about each input/output parameter can go here
        input_file: "input file for analysis"
    }

    input {
        # task inputs are declared here
        File input_file
        RuntimeAttr? runtime_attr_override
    }

    # other "private" declarations can be made here
    #Int disk_size = 365 + 2 * ceil(size(fastq_files, "GB")) # change this to whatever you need.

    command <<<
        # the command template - this section is required
        # put your shell commands here.
        cp ~{input_file} second_task_output.txt
    >>>

    output {
        # task outputs are declared here
        File results = "second_task_output.txt"
    }
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            100, #disk_size,
        boot_disk_gb:       50,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "mjfos2r/mycontainer:latest"
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
