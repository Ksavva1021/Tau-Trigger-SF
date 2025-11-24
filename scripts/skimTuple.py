#!/usr/bin/env python

import argparse
from array import array
import math
import numpy as np
import os
import re
import sys
import ROOT
import glob

# -------------------------
# Example commands:
# python3 scripts/skimTuple.py --input_dir "dir" --selection DeepTau --type mc --pudata pudata.root --pumc pumc.root --output dir/skim_mc.root
# python3 scripts/skimTuple.py --input_dir "dir" --selection DeepTau --type data --output dir/skim_data.root
# ------------------------


parser = argparse.ArgumentParser(description='Skim full tuple.')
parser.add_argument('--input', required=False, type=str, nargs='+', help="input files")
parser.add_argument('--input_dir', required=False, type=str, help="input directory containing sub-directories with ROOT files (works with a wildcard)")
parser.add_argument('--selection', required=True, type=str, help="tau selection")
parser.add_argument('--output', required=True, type=str, help="output file")
parser.add_argument('--type', required=True, type=str, help="data or mc")
parser.add_argument('--pudata', required=False, type=str, default=None,
                    help="file with the pileup profile for the data taking period")
parser.add_argument('--pumc', required=False, type=str, default=None,
                    help="file with the pileup profile for the data taking period")
args = parser.parse_args()

sys.path.insert(0, 'Common/python')
from AnalysisTypes import *
from AnalysisTools import *
import TriggerConfig
ROOT.ROOT.EnableImplicitMT(4)
ROOT.gROOT.SetBatch(True)
ROOT.gInterpreter.Declare('#include "interface/PyInterface.h"')
ROOT.gInterpreter.Declare('#include "interface/picoNtupler.h"')

LUMI_PB = 7980.4
DY_CFG = {
    "dy_NLO_0J": {"xs": 1664.6841730883225, "eff": 31651962.0, "filter_eff": 1.0},
    "dy_NLO_1J": {"xs":  316.2403378778898, "eff": 26576552.0, "filter_eff": 1.0},
    "dy_NLO_2J": {"xs":  116.47203023117962, "eff": 32763798.0, "filter_eff": 1.0},
}

input_files = []

if args.input:
    input_files.extend(args.input)

if args.input_dir and ',' in args.input_dir:
    input_dirs = args.input_dir.split(',')
else:
    input_dirs = [args.input_dir] if args.input_dir else []

for input_dir in input_dirs:
    if '*' in input_dir or '?' in input_dir or '[' in input_dir:
        matched_dirs = glob.glob(input_dir)
        for matched_dir in matched_dirs:
            if not os.path.isdir(matched_dir):
                raise RuntimeError(f"Invalid directory: {matched_dir}")
            for root, _, files in os.walk(matched_dir):
                dir_files = [os.path.join(root, f) for f in files if f.endswith('.root')]
                input_files.extend(dir_files)
    else:
        if not os.path.isdir(input_dir):
            raise RuntimeError(f"Invalid directory: {input_dir}")
        for root, _, files in os.walk(input_dir):
            dir_files = [os.path.join(root, f) for f in files if f.endswith('.root')]
            input_files.extend(dir_files)

if not input_files:
    raise RuntimeError("No input files provided. Use --input or --input-dir to specify ROOT files.")

print(input_files)
files_0J = [f for f in input_files if "/dy_NLO_0J/" in f]
files_1J = [f for f in input_files if "/dy_NLO_1J/" in f]
files_2J = [f for f in input_files if "/dy_NLO_2J/" in f]
other    = [f for f in input_files if f not in files_0J+files_1J+files_2J]

dfs = []
if args.type == 'mc':
    input_0L = ListToStdVector(files_0J)
    df_0J = ROOT.RDataFrame('Events', input_0L)
    df_0J = df_0J.Define('wt_sf', f"{LUMI_PB*DY_CFG['dy_NLO_0J']['xs']/DY_CFG['dy_NLO_0J']['eff']}")
    df_0J = df_0J.Redefine('genWeight', 'genWeight > 0 ? 1 : -1')
    df_0J = df_0J.Redefine('wt_sf', f'wt_sf * genWeight')
    input_1L = ListToStdVector(files_1J)
    df_1J = ROOT.RDataFrame('Events', input_1L)
    df_1J = df_1J.Define('wt_sf', f"{LUMI_PB*DY_CFG['dy_NLO_1J']['xs']/DY_CFG['dy_NLO_1J']['eff']}")
    df_1J = df_1J.Redefine('genWeight', 'genWeight > 0 ? 1 : -1')
    df_1J = df_1J.Redefine('wt_sf', f'wt_sf * genWeight')
    input_2L = ListToStdVector(files_2J)
    df_2J = ROOT.RDataFrame('Events', input_2L)
    df_2J = df_2J.Define('wt_sf', f"{LUMI_PB*DY_CFG['dy_NLO_2J']['xs']/DY_CFG['dy_NLO_2J']['eff']}")
    df_2J = df_2J.Redefine('genWeight', 'genWeight > 0 ? 1 : -1')
    df_2J = df_2J.Redefine('wt_sf', f'wt_sf * genWeight')
    dfs.extend([df_0J, df_1J, df_2J])
else:
    input_other = ListToStdVector(other)
    df_other = ROOT.RDataFrame('Events', input_other)
    df_other = df_other.Define('genWeight', '1.0')
    df_other = df_other.Define('wt_sf', "1.0")
    dfs.append(df_other)

if args.type not in ['data', 'mc']:
    raise RuntimeError("Invalid sample type")

input_vec = ListToStdVector(input_files)
if args.type == 'mc':
    if args.pudata is None or args.pumc is None:
        raise RuntimeError("Pileup file should be provided for mc.")
    data_pu_file = ROOT.TFile(args.pudata, 'READ')
    data_pu = data_pu_file.Get('pileup')
    #df_all = ROOT.RDataFrame('Events', input_vec)
    mc_pu_file = ROOT.TFile(args.pumc, 'READ')
    mc_pu = mc_pu_file.Get('pileup')
    ROOT.PileUpWeightProvider.Initialize(data_pu, mc_pu)


selection_id = ParseEnum(TauSelection, args.selection)
df = ROOT.RDataFrame('Events', input_vec)
# instead concatenate all dataframes
for i,df in enumerate(dfs):
    df = df.Filter('''
                (tau_sel & {}) != 0 && muon_pt > 24 && muon_iso < 0.1 && muon_mt < 30
                && tau_pt > 20 && abs(tau_eta) < 2.1 && tau_decayMode != 5 && tau_decayMode != 6
                && vis_mass > 40 && vis_mass < 80
                '''.format(selection_id))
    if selection_id == TauSelection.DeepTau:
        df = df.Filter('( tau_idDeepTau2018v2p5VSmu  & 4) != 0')
    if args.type == 'mc':
        df = df.Filter('tau_charge + muon_charge == 0 && tau_gen_match == 5')
        df = df.Define('weight', "PileUpWeightProvider::GetDefault().GetWeight(npu) * 1.0")
    else:
        df = df.Define('weight', "muon_charge != tau_charge ? 1. : -1.")

    skimmed_branches = [
        'muon_pt', 'tau_pt', 'tau_eta', 'tau_phi', 'tau_mass', 'tau_charge', 'tau_decayMode','tau_decayModePNet', 'weight', 'tau_idDeepTau2017v2p1VSjet', 'tau_idDeepTau2018v2p5VSjet',"tau_ipLengthSig","tau_hasRefitSV",'TrigObj_l1pt', 'TrigObj_l1iso', 'nTrigObj'
        # use monitoring path, as TnP won't work in HLT path
    ]
    # Adding branches for studies
    skimmed_branches.extend(['vis_mass', 'muon_mt', 'wt_sf', 'genWeight'])

    df = df.Define("pass_mutau", "PassMuTauTrig2022(nTrigObj, TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, tau_pt, tau_eta, tau_phi)")
    df = df.Define("pass_etau", "PassEleTauTrig2022(nTrigObj, TrigObj_l1pt, TrigObj_l1iso, TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, tau_pt, tau_eta, tau_phi)")

    # ditau -> drop !bit18 cut, change l1pt>32 with l1pt>=32
    df = df.Define("pass_ditau", "PassDiTauTrig2022(nTrigObj, TrigObj_l1pt, TrigObj_l1iso, TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, tau_pt, tau_eta, tau_phi)")
    # 
    df = df.Define("pass_ditaujet", "PassDiTauJetTrig2022(nTrigObj, TrigObj_l1pt, TrigObj_l1iso, TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, tau_pt, tau_eta, tau_phi)")

    # VBFSingleTau and VBFDiTau
    df = df.Define("pass_vbftau", "PassVBFTauTrig2022(nTrigObj, TrigObj_l1pt, TrigObj_l1iso, TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, tau_pt, tau_eta, tau_phi)")
    df = df.Define("pass_vbfditau", "PassVBFDiTauTrig2022(nTrigObj, TrigObj_l1pt, TrigObj_l1iso,TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, tau_pt, tau_eta, tau_phi)")

    # Studies with only L1
    df = df.Define("pass_mutau_l1", "PassMuTauTrig2022_L1Only(nTrigObj, TrigObj_l1pt, TrigObj_l1iso, TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, tau_pt, tau_eta, tau_phi)")
    df = df.Define("pass_etau_l1", "PassETauTrig2022_L1Only(nTrigObj, TrigObj_l1pt, TrigObj_l1iso, TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, tau_pt, tau_eta, tau_phi)")
    df = df.Define("pass_ditau_l1", "PassDiTauTrig2022_L1Only(nTrigObj, TrigObj_l1pt, TrigObj_l1iso, TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, tau_pt, tau_eta, tau_phi)")
    df = df.Define("pass_ditaujet_l1", "PassDiTauJetTrig2022_L1Only(nTrigObj, TrigObj_l1pt, TrigObj_l1iso, TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, tau_pt, tau_eta, tau_phi)")
    df = df.Define("pass_vbftau_l1", "PassVBFTauTrig2022_L1Only(nTrigObj, TrigObj_l1pt, TrigObj_l1iso, TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, tau_pt, tau_eta, tau_phi)")
    df = df.Define("pass_vbfditau_l1", "PassVBFDiTauTrig2022_L1Only(nTrigObj, TrigObj_l1pt, TrigObj_l1iso, TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, tau_pt, tau_eta, tau_phi)")


    skimmed_branches.append("pass_ditau")
    skimmed_branches.append("pass_etau")
    skimmed_branches.append("pass_mutau")
    skimmed_branches.append("pass_ditaujet")
    skimmed_branches.append("pass_vbftau")
    skimmed_branches.append("pass_vbfditau")
    skimmed_branches.append("pass_mutau_l1")
    skimmed_branches.append("pass_etau_l1")
    skimmed_branches.append("pass_ditau_l1")
    skimmed_branches.append("pass_ditaujet_l1")
    skimmed_branches.append("pass_vbftau_l1")
    skimmed_branches.append("pass_vbfditau_l1")

    #Create output directory if it does not exist
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # snapshot_options = ROOT.RDF.RSnapshotOptions()
    # snapshot_options.fMode = "UPDATE" if i > 0 else "RECREATE"
    file_name = args.output.replace('.root', f'_{i}.root')
    df.Snapshot('Events', file_name, ListToStdVector(skimmed_branches))
# hadd the output files
if len(dfs) > 1:
    os.system(f"hadd -f {args.output} " + ' '.join([args.output.replace('.root', f'_{i}.root') for i in range(len(dfs))]))
    for i in range(len(dfs)):
        os.remove(args.output.replace('.root', f'_{i}.root'))
print("Check Point")
os._exit(0)
