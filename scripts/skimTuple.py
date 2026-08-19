#!/usr/bin/env python

import argparse
from array import array
import math
import numpy as np
import os
import re
import sys
import ROOT

# -------------------------
# Example commands:
# python3 skimTuple.py --input_dir dir1,dir2,dir3 --selection DeepTau --type mc --pudata pudata.root --pumc pumc.root --output dir/skim_mc.root
# python3 skimTuple.py --input_dir dir1_data,dir2_data --selection DeepTau --type data --output dir/skim_data.root
# ------------------------


parser = argparse.ArgumentParser(description='Skim full tuple.')
parser.add_argument('--input', required=False, type=str, nargs='+', help="input files")
parser.add_argument('--input_dir', required=False, type=str, help="input directory/ies containing ROOT files (comma-separated)")
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

input_files = []

if args.input:
    input_files.extend(args.input)

if args.input_dir:
    for input_dir in args.input_dir.split(','):
        if not os.path.isdir(input_dir):
            raise RuntimeError(f"Invalid directory: {input_dir}")
        dir_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith('.root')]
        input_files.extend(dir_files)

if not input_files:
    raise RuntimeError("No input files provided. Use --input or --input-dir to specify ROOT files.")

print(input_files)

if args.type not in ['data', 'mc']:
    raise RuntimeError("Invalid sample type")

input_vec = ListToStdVector(input_files)
if args.type == 'mc':
    if args.pudata is None or args.pumc is None:
        raise RuntimeError("Pileup file should be provided for mc.")
    data_pu_file = ROOT.TFile(args.pudata, 'READ')
    data_pu = data_pu_file.Get('pileup')
    df_all = ROOT.RDataFrame('Events', input_vec)
    mc_pu_file = ROOT.TFile(args.pumc, 'READ')
    mc_pu = mc_pu_file.Get('pileup')
    ROOT.PileUpWeightProvider.Initialize(data_pu, mc_pu)


selection_id = ParseEnum(TauSelection, args.selection)
df = ROOT.RDataFrame('Events', input_vec)
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
    'tau_pt', 'tau_eta', 'tau_phi', 'tau_mass', 'tau_charge', 'tau_decayMode',
    'tau_decayModePNet', 'weight', 'tau_idDeepTau2017v2p1VSjet',
    'tau_idDeepTau2018v2p5VSjet', "tau_ipLengthSig", "tau_hasRefitSV",
    'TrigObj_l1pt', 'TrigObj_l1iso', 'nTrigObj'
    # use monitoring path, as TnP won't work in HLT path
]

df = df.Define(
    "pass_mutau",
    "PassMuTauTrig2022(nTrigObj, TrigObj_id, TrigObj_filterBits,TrigObj_pt, TrigObj_eta, TrigObj_phi, tau_pt, tau_eta, tau_phi)"
)
df = df.Define(
    "pass_etau",
    "PassEleTauTrig2022(nTrigObj, TrigObj_l1pt, TrigObj_l1iso, TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, tau_pt, tau_eta, tau_phi)"
)
# ditau -> drop !bit18 cut, change l1pt>32 with l1pt>=32
df = df.Define(
    "pass_ditau",
    "PassDiTauTrig2022(nTrigObj, TrigObj_l1pt, TrigObj_l1iso, TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, tau_pt, tau_eta, tau_phi)"
) 
df = df.Define(
    "pass_ditaujet",
    "PassDiTauJetTrig2022(nTrigObj, TrigObj_l1pt, TrigObj_l1iso, TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, tau_pt, tau_eta, tau_phi)"
)

skimmed_branches.append("pass_ditau")
skimmed_branches.append("pass_etau")
skimmed_branches.append("pass_mutau")
skimmed_branches.append("pass_ditaujet")

df.Snapshot('Events', args.output, ListToStdVector(skimmed_branches))
print("Check Point")
os._exit(0)
