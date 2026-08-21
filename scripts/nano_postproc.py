#!/usr/bin/env python3
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from importlib import import_module
import os
import sys; sys.path.append('python')
import ROOT
import argparse
import pickle
import json

from summaryProducer import *
from selectionFilter import *
from tupleProducer import *
from submitJob import condor_submit

ROOT.PyConfig.IgnoreCommandLineOptions = True

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Post Processing.')
    parser.add_argument('--input', required=False, type=str, help="NANO input")
    parser.add_argument('--inputJson', required=False, type=str, help="NANO input file list (HiggsDNA style sample json)")
    parser.add_argument('--output', required=True, type=str, help="Output directory")
    parser.add_argument('--isMC', action='store_true', help="is it MC?")
    parser.add_argument('--era', required=True, type=str, help="2022, 2023, 2024")
    parser.add_argument('--batch', action='store_true', help="run in batch mode")
    args = parser.parse_args()
    print("args = ",args)

    if (args.input is None) and (args.inputJson is None):
        raise RuntimeError("Please check the input!")
    if (args.input is not None) and (args.inputJson is not None):
        raise RuntimeError("Please check the input!")

    isMC = args.isMC
    output = args.output
    era = args.era
    if args.input:
        files = [ args.input ]
    if args.inputJson:
        with open(args.inputJson, 'r') as f:
            samples = json.load(f)
        if isMC:
            files = samples['DYto2Tau_MLL_50_amcatnloFXFX']
        else:
            files = []
            for key in samples.keys():
                if 'Muon' in key:  # Muon datasets
                    print(f"Adding files from key: {key}")
                    files.extend(samples[key])

    if era == "2022EE": era = "2022"
    if era == "2023BPix": era = "2023"

    if not args.batch:
        if isMC:
            jsoninput = None,
            if era == '2016':
                Modules = [summary2016MC(), selection2016MC(), tuple2016MC()]
            elif era == '2017':
                Modules = [summary2017MC(), selection2017MC(), tuple2017MC()]
            elif era == '2018':
                Modules = [summary2018MC(), selection2018MC(), tuple2018MC()]
            elif era == '2022':
                Modules = [summary2022MC(), selection2022MC(), tuple2022MC()]
            elif era == '2023':
                Modules = [summary2023MC(), selection2023MC(), tuple2023MC()]
            elif era == '2024':
                Modules = [summary2024MC(), selection2024MC(), tuple2024MC()]
            else:
                raise RuntimeError("Please check the right Year!")
            p = PostProcessor(output, files, "1", 
                              branchsel = "keep_and_drop.txt", 
                              modules= Modules, 
                              provenance=True,
                              outputbranchsel = "output_branch.txt"
            )

        else:
            if era == '2022':
                Modules = [summary2022data(), selection2022data(), tuple2022data()]
                lumi_json_path = './lumi_jsons/2022.pkl' 
            elif era == '2023':
                Modules = [summary2023data(), selection2023data(), tuple2023data()]
                lumi_json_path = './lumi_jsons/2023.pkl'
            elif era == '2024':
                Modules = [summary2024data(), selection2024data(), tuple2024data()]
                lumi_json_path = './lumi_jsons/Cert_Collisions2024_378981_386951_Golden.json'
            else:
                raise RuntimeError("Please check the right Year!")

            if lumi_json_path.endswith('.pkl'):
                with open(lumi_json_path, 'rb') as file:
                    jsoninput = pickle.load(file)
            elif lumi_json_path.endswith('.json'):
                with open(lumi_json_path, 'r') as file:
                    jsoninput = json.load(file)

            p = PostProcessor(output, files, "1", 
                              branchsel = "keep_and_drop.txt", 
                              modules= Modules, 
                              jsonInput=jsoninput,
                              provenance=True,
                              outputbranchsel = "output_branch.txt"
            )

        p.run()
        print("Done !")

    else:
        BATCH_SIZE = 5  # Number of files to process in each batch
        commands = []
        if args.isMC:
            for file in files:
                commands.append(
                    f"python3 scripts/nano_postproc.py --input {file} --isMC --era {args.era} --output {args.output}"
                )
        else:
            for file in files:
                commands.append(
                    f"python3 scripts/nano_postproc.py --input {file} --era {args.era} --output {args.output}"
                )
        
        for i in range(0, len(commands), BATCH_SIZE):
            batch_commands = commands[i:min(i + BATCH_SIZE, len(commands))]
            if args.isMC:
                job_name = f"nano_postproc_mc_{args.era}_{i//BATCH_SIZE}"
            else:
                job_name = f"nano_postproc_data_{args.era}_{i//BATCH_SIZE}"
            print(
                "\033[92m"
                + f"Submitting job: {job_name}"
                + "\033[0m"
            )
            condor_dir = f"{args.output}/condor"
            os.makedirs(condor_dir, exist_ok=True)
            condor_submit(batch_commands, condor_dir, job_name)

