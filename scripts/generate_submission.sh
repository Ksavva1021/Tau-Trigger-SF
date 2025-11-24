#!/bin/bash

# Useful Commands:
# ------------------
# If you want to create multiple sub files at the same time for different samples you can do the following:

# year=2023BPix  # Set the year variable
# 
# for data_type in $(ls Run3_"$year"/*.txt | xargs -n 1 basename -s .txt); do
#     echo "Running command for: $data_type"
#     bash scripts/generate_submission.sh --year "$year" --data_type "$data_type"
# done

# ------------------
# Generate sub files: bash scripts/generate_submission.sh --year 2023BPix --data_type dy*
# Submit all sub files: find submissions -name "*.sub" | xargs -I {} condor_submit {} -spool

# Default values
year=""
data_type=""

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --year) year="$2"; shift ;;   # Assign year from argument
        --data_type) data_type="$2"; shift ;; # Assign data_type from argument
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Check if required arguments are provided
if [[ -z "$year" || -z "$data_type" ]]; then
    echo "Usage: $0 --year <YEAR> --data_type <DATA_TYPE>"
    exit 1
fi

# Handle wildcard for data_type
if [[ "$data_type" == *"*"* ]]; then
    for expanded_data_type in $(ls Run3_${year}/${data_type}.txt | xargs -n 1 basename -s .txt); do
        echo "Processing data_type: $expanded_data_type"
        bash "$0" --year "$year" --data_type "$expanded_data_type"
    done
    exit 0
fi

# Determine isMC based on data_type
if [[ "$data_type" =~ ^dy ]]; then
    isMC=1
else
    isMC=0
fi

# Create the submissions directory if it doesn't exist
mkdir -p submissions/Run3_${year}
mkdir -p logs/Run3_${year}/${data_type}
mkdir -p /eos/cms/store/group/phys_tau/irandreo/TauTrgSF/Run3_${year}/${data_type}

# Generate the submission file
cat <<EOF > submissions/Run3_${year}/job_submission_${data_type}.sub
executable              = scripts/run_job.sh
arguments               = --input \$(inputFile) --output /eos/cms/store/group/phys_tau/irandreo/TauTrgSF/Run3_${year}/${data_type}/ --isMC ${isMC} --era ${year}
log                     = logs/Run3_${year}/${data_type}/\$(ProcId).log
error                   = logs/Run3_${year}/${data_type}/\$(ProcId).err
output                  = logs/Run3_${year}/${data_type}/\$(ProcId).out

should_transfer_files = YES
when_to_transfer_output = ON_EXIT

# Specify files to transfer
transfer_output_files   = _condor_stdout, _condor_stderr
transfer_output_remaps  = "_condor_stdout=/eos/user/i/irandreo/TauTrgSF/PhysicsTools/Tau-Trigger-SF/logs/Run3_${year}/${data_type}/\$(ProcId).out, _condor_stderr=/eos/user/i/irandreo/TauTrgSF/PhysicsTools/Tau-Trigger-SF/logs/Run3_${year}/${data_type}/\$(ProcId).err"

request_memory = 8GB

# Job runtime flavor
+JobFlavour             = "longlunch"

# Queue jobs from input file list
queue inputFile from Run3_${year}/${data_type}.txt
EOF

echo "Submission file 'submissions/Run3_${year}/job_submission_${data_type}.sub' has been generated."

