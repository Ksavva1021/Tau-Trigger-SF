import os
import subprocess
import shutil

# ---------------
# Example Command to Run: python3 logs/check_logs.py logs/Run3_2022/* --era Run3_2022 --output /eos/user/i/irandreo/TauTrgSF/PhysicsTools/Tau-Trigger-SF
# Resubmission needs editing given new types of MC wildcards
# --------------

def create_sh_file(resubmit_path, out_filename, txt_file_path, isMC, era, output_dir):
    """
    Creates a .sh file for local submission based on the .out file name.

    Args:
        out_filename: The name of the .out file (e.g., "nano_0.out").
        txt_file_path: Path to the .txt file containing input paths.
        isMC: Integer (0 for data, 1 for dy).
        era: The era for the command (e.g., "Run3_2022").
        output_dir: Directory for the command output.
    """
    # Extract the number in front of the .out filename
    try:
        file_number = int(out_filename.split("_")[-1].replace(".out", ""))
    except ValueError:
        print(f"Error: Unable to extract file number from {out_filename}")
        return

    # Read the corresponding line from the txt file
    try:
        with open(txt_file_path, "r") as f:
            lines = f.readlines()
            if file_number >= len(lines):
                print(f"Error: File number {file_number} exceeds the number of lines in {txt_file_path}")
                return
            input_path = lines[file_number].strip()
    except Exception as e:
        print(f"Error reading {txt_file_path}: {e}")
        return

    if era in ["Run3_2022","Run3_2022EE"]:
        era = "2022"
    elif era in ["Run3_2023", "Run3_2023BPix"]:
        era = "2023"

    cwd = os.getcwd()

    if "NanoAODTools" in cwd:
        dir_NanoAODTools = cwd[:cwd.index("NanoAODTools") + len("NanoAODTools")]
    else:
        raise ValueError("NanoAODTools not found in the path")

    # Construct the command
    command = (
        f"python3 scripts/nano_postproc.py --input {input_path} "
        f"--isMC {isMC} --era {era} --output {output_dir}/"
    )

    # Create the .sh file
    sh_filename = out_filename.replace(".out", ".sh")
    sh_filename = resubmit_path + "/" + sh_filename
    try:
        with open(sh_filename, "w") as sh_file:
            sh_file.write("#!/bin/bash\n")
            sh_file.write("set -e # Exit on errors\n")
            sh_file.write("set -x # Debug mode: print commands as they are executed\n")
            sh_file.write(f"cd /eos/user/i/irandreo/TauTrgSF/PhysicsTools/NanoAODTools\n")
            sh_file.write("source standalone/env_standalone.sh\n")
            sh_file.write(f"cd /eos/user/i/irandreo/TauTrgSF/PhysicsTools/Tau-Trigger-SF\n")
            sh_file.write(command + "\n")
        os.chmod(sh_filename, 0o755)  # Make the .sh file executable
    except Exception as e:
        print(f"Error creating {sh_filename}: {e}")


def create_sub_file(resubmit_path, out_filename, year, data_type):
    """
    Creates an HTCondor .sub file for job submission, ensuring log filenames match the job script.

    Args:
        resubmit_path (str): Directory where the .sub file will be created.
        sh_filename (str): The name of the .sh script (e.g., "18.sh").
        year (str): The dataset year (e.g., "2024").
        data_type (str): The type of data (e.g., "MC" or "Data").
    """

    # Extract the number in front of the .out filename
    try:
        file_number = int(out_filename.split("_")[-1].replace(".out", ""))
    except ValueError:
        print(f"Error: Unable to extract file number from {out_filename}")
        return

    # Define the .sub file name and its full path
    sub_filename = f"{file_number}.sub"
    sub_filepath = os.path.join(resubmit_path, sub_filename)
    sh_filepath = sub_filepath.replace(".sub",".sh")

    # Define log, output, and error file paths
    log_dir = f"logs/{year}/{data_type}"
    log_file = f"{log_dir}/{file_number}.log"
    error_file = f"{log_dir}/{file_number}.err"
    output_file = f"{log_dir}/{file_number}.out"

    # Create the HTCondor submission file content (split into sections for readability)
    sub_content = [
        f"executable              = {sh_filepath}",
        f"log                     = {log_file}",
        f"error                   = {error_file}",
        f"output                  = {output_file}",
        "",
        "should_transfer_files = YES",
        "when_to_transfer_output = ON_EXIT",
        "# Specify files to transfer",
        "transfer_output_files   = _condor_stdout, _condor_stderr",
        f"transfer_output_remaps  = \"_condor_stdout=/eos/user/i/irandreo/TauTrgSF/PhysicsTools/Tau-Trigger-SF/{output_file}, _condor_stderr=/eos/user/i/irandreo/TauTrgSF/PhysicsTools/Tau-Trigger-SF/{error_file}\"",

        "request_memory = 8GB",

        "# Job runtime flavor",
        '+JobFlavour             = "longlunch"',
        "",
        "queue"
    ]

    # Write the .sub file
    try:
        with open(sub_filepath, "w") as sub_file:
            sub_file.write("\n".join(sub_content) + "\n")
        print(f"HTCondor submission file created: {sub_filepath}")
    except Exception as e:
        print(f"Error creating {sub_filepath}: {e}")

    return sub_filepath


def process_directory(directory_path, era, output_dir):
    """
    Processes a directory to find and handle .sub files and create .sh files.

    Args:
        directory_path: The path to the directory to process.
        era: The era for the command (e.g., "2022").
        output_dir: Directory for the command output.

    Returns:
        A tuple containing the number of successes and the total number of files processed.
    """
    success_count = 0
    total_count = 0

    # Create the resubmit directory if it doesn't exist
    resubmit_dir = os.path.join(directory_path, "resubmit")
    os.makedirs(resubmit_dir, exist_ok=True)

    sample_type = directory_path.split("/")[-1]
    if sample_type == "":
        sample_type = directory_path.split("/")[-2]

    output_dir = output_dir + "/" + era + "/" + sample_type

    for filename in os.listdir(directory_path):
        if filename.endswith(".log"):
            total_count += 1
            try:
                # Replace .log with .out
                out_filename = filename.replace(".log", ".out")
                out_filepath = os.path.join(directory_path, out_filename)

                cwd = os.getcwd() 
                # Determine the txt file and isMC based on the log directory
                if sample_type == "dy":
                    txt_file = os.path.join(cwd, era, "dy.txt")
                    isMC = 1
                else:
                    txt_file = os.path.join(cwd, era, "{}.txt".format(sample_type))
                    isMC = 0

                if os.path.exists(out_filepath):
                    # Check for "Done !" in the .out file
                    with open(out_filepath, "r") as f:
                        if "Done !" in f.read():
                            success_count += 1
                        else:
                            # Copy the .out file to the resubmit directory
                            shutil.copy(out_filepath, os.path.join(resubmit_dir, out_filename))
                            # Create a .sh file
                            create_sh_file(resubmit_dir,out_filename, txt_file, isMC, era, output_dir)
                            sub_filepath = create_sub_file(resubmit_dir, out_filename,era,sample_type)
                            subprocess.run(["condor_submit", "-spool", sub_filepath], check=True)
                else:
                    # Create an empty .out file in the resubmit directory
                    empty_out_filepath = os.path.join(resubmit_dir, out_filename)
                    open(empty_out_filepath, "w").close()
                    # Create a .sh file
                    create_sh_file(resubmit_dir,out_filename, txt_file, isMC, era, output_dir)
                    sub_filepath = create_sub_file(resubmit_dir, out_filename,era,sample_type)
                    subprocess.run(["condor_submit", "-spool", sub_filepath], check=True)


            except Exception as e:
                print(f"Error processing {filename}: {e}")

    return success_count, total_count

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Process a directory of .log files and handle corresponding .out files.")
    parser.add_argument("directory_path", nargs="+", help="Path to the directory to process.")
    parser.add_argument("--era", required=True, help="Era for the command (e.g., 2022).")
    parser.add_argument("--output_dir", required=True, help="Directory for the command output.")
    args = parser.parse_args()

    total_successes = 0
    total_files = 0
    for directory_path in args.directory_path:
        successes, total = process_directory(directory_path, args.era, args.output_dir)
        total_successes += successes
        total_files += total
        
        print(f"Processed directory: {directory_path}")
        print(f"Successes: {successes}/{total}")
    
    print(f"Total Successes: {total_successes}/{total_files}")
