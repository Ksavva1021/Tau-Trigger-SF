import subprocess

def condor_submit(commands, condor_dir, job_name):
    # Create executable shell script for the command
    with open(f"{condor_dir}/{job_name}.sh", "w") as f:
        f.write("#!/bin/bash\n")
        for command in commands:
            f.write(f"{command}\n")
    subprocess.run(["chmod", "+x", f"{condor_dir}/{job_name}.sh"])

    # Create condor submission file
    with open(f"{condor_dir}/{job_name}.sub", "w") as f:
        f.write("getenv = true\n")
        f.write(f"executable = {condor_dir}/{job_name}.sh\n")
        f.write(f"output = {condor_dir}/{job_name}.$(Cluster).out\n")
        f.write(f"error = {condor_dir}/{job_name}.$(Cluster).err\n")
        f.write(f"log = {condor_dir}/{job_name}.$(Cluster).log\n")
        f.write("request_cpus = 1\n")
        f.write("request_memory = 2GB\n")
        f.write("+MaxRuntime = 86400\n")  # Set max runtime to 24 hours
        f.write("queue\n")

    subprocess.run(["condor_submit", f"{condor_dir}/{job_name}.sub"])