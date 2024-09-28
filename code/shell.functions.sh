#!/bin/bash
# Create a function for waiting job to finish (avoid to use all memory of the group in the server)
# made by Giovanni Ladeira, July 12, 2024
# Thank you Giovanni
function wait_for_job_completion() {
    # Crete a local variable (job_id) to store the first argument passed to the function
    local job_id=$1
    local time_min=$2
    echo "          Waiting for job $job_id to complete..."
    # While loop only stops when break or exit commands are executed
    while true; do
        # Check if the job is still in the queue
        job_status=$(squeue -j $job_id -h -o %T 2>/dev/null)
        if [[ -z $job_status ]]; then
            echo ""
            break
        elif [[ $job_status == "FAILED" ]]; then
            echo "          Job $job_id failed."
            exit 1
        else
            # Job is still running or pending
            echo "          Job $job_id is $job_status. Checking again in ${time_min} minutes..."
            sleep ${time_min}m
        fi
    done
}

function wait_for_job_completion_less_output() {
    # Crete a local variable (job_id) to store the first argument passed to the function
    local job_id=$1
    local time_min=$2
    # While loop only stops when break or exit commands are executed
    while true; do
        # Check if the job is still in the queue
        job_status=$(squeue -j $job_id -h -o %T 2>/dev/null)
        if [[ -z $job_status ]]; then
            break
        elif [[ $job_status == "FAILED" ]]; then
            echo "          Job $job_id failed."
            exit 1
        else
            # Job is still running or pending
            echo "          Job $job_id is $job_status. Checking again in ${time_min} minutes..."
            sleep ${time_min}
        fi
    done
}