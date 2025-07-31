version 1.0
import "../structs/Structs.wdl"
import "../tasks/Tasks.wdl" as EXAMPLE

workflow WorkflowName {

    meta {
        description: "Description of the workflow"
    }
    parameter_meta {
        input_file: "description of input"
    }

    input {
        File input_file
    }

    # call our first task/workflow
    call EXAMPLE.FirstTask {
        input:
            input_file = input_file,
    }
    # pass the output of first task into our second task/workflow.
    call EXAMPLE.SecondTask {
        input:
            input_file = FirstTask.output_file,
    }

    output {
        File output_file = SecondTask.results
    }
}