#! /usr/bin/env python3

import sys; sys.path.append('python')
from utils import *
import uproot
import os, ROOT
import numpy

from correctionlib.schemav2 import (
    VERSION,
    Binning,
    Category,
    Correction,
    CorrectionSet,
    Formula,
    Transform
)

# ------------------
# Example Command:
# python3 scripts/convert_to_json.py --years 2022preEE 2022postEE 2023preBPix 2023postBPix --outdir jsons
# ------------------

#
# Some default values
#
default_workingpoints     = [
  'VVVLoose', 'VVLoose', 'VLoose',
  'Loose', 'Medium', 'Tight',
  'VTight', 'VVTight'
]
corrtypes = ['sf', 'eff_mc', 'eff_data']
trigtypes = {
    '2022preEE' : ['ditau', 'etau', 'mutau', 'ditaujet'],
    '2022postEE' : ['ditau', 'etau', 'mutau', 'ditaujet'],
    '2023preBPix' : ['ditau', 'etau', 'mutau', 'ditaujet'],
    '2023postBPix' : ['ditau', 'etau', 'mutau', 'ditaujet'],
    '2024' : ['ditau', 'etau', 'mutau', 'ditaujet', 'ditauANDditaujet']
}

dms_nonmerged = [-1, 0, 1, 2, 10, 11]
dms_merged = [-1, 0, 1, 2, 10]

dm_dict_nonmerged = {
  -1: 'dmall',
  0 : 'dm0',
  1 : 'dm1',
  2 : 'dm2',
  10 : 'dm10',
  11 : 'dm11'
}
dm_dict_merged = {
  -1: 'dmall',
  0 : 'dm0',
  1 : 'dm1',
  2 : 'dm2',
  10 : 'dm1011',
  11 : 'dm1011'
}

corrtype_dict = {
  'sf' : 'sf',
  'eff_mc' : 'mc',
  'eff_data' : 'data',
}

idalgo_dict = {
  '2022preEE' : 'DeepTau2018v2p5',
  '2022postEE' : 'DeepTau2018v2p5',
  '2023preBPix' : 'DeepTau2018v2p5',
  '2023postBPix' : 'DeepTau2018v2p5',
  '2024' : 'PNet'
}


def in_file_name(year):
  return f'jsons/fitTurnOn_{year}.root'


def in_hist_name(corrtype, typ, wp, dm_str): 
  return '_'.join([corrtype,typ,wp,dm_str,'fitted'])


# Helper function to merge bins that have similar values
def merge_pt_bins(edges, values, errors, pt_threshold = 20.):

  tmp_values = []
  tmp_errors = []

  bin_nb = 0
  n_merges = 0
  original_val = 0
  for val, err in zip(values, errors):
    if edges[bin_nb] < pt_threshold:
      edges.pop(bin_nb)
      continue

    if bin_nb > 0 and val == tmp_values[-1] and err == tmp_errors[-1]:
      edges.pop(bin_nb)
    elif bin_nb > 0 and abs((val-original_val)/original_val) < 0.005 and n_merges < 4:
      tmp_values[-1] = (tmp_values[-1] + val)/2
      tmp_errors[-1] = max(err, tmp_errors[-1])
      edges.pop(bin_nb)
      n_merges += 1
    else:
      tmp_values.append(val)
      tmp_errors.append(err)
      bin_nb += 1
      n_merges = 0
      original_val = val
  
  return edges, tmp_values, tmp_errors


# Helper function to unpack kwargs
def kwargs_get(kwargs, kw, default):
  val = kwargs.get(kw)
  if val is None: val = default
  return val


def getPtThreshold(triggertype):
  if triggertype == 'mutau':
    return 32.
  else:
    return 35.


# Functions to build the correction objects
def build_pts(in_file, hist_name, pt_threshold = 20.):

  f = uproot.open(in_file)
  hist = f[hist_name]
  edges = [round(x,5) for x in hist.to_numpy()[1]]

  edges, tmp_values, tmp_errors = merge_pt_bins(edges, hist.values(), hist.errors(), pt_threshold = pt_threshold)

  content = {
    'nom' : [round(val, 8) for val in tmp_values],
    'up' : [round(val+err, 8) for val, err in zip(tmp_values, tmp_errors)],
    'down' : [round(val-err, 8) for val, err in zip(tmp_values, tmp_errors)]
  }

  edges[-1] = float('inf')
  return Category.parse_obj(
    {
      "nodetype": "category",
      "input": "syst",
      "content": [
        {"key": syst, 
        "value": {
            "nodetype": "binning",
            "input": "pt",
            "edges": edges,
            "content": content[syst],
            "flow": "error",
          } 
        } for syst in ['nom', 'up', 'down']
      ],
    }
  )


def build_dms(trigtype, wp, year, corrtype):
  print('Filling {0} {1} trigger for {2} WP'.format(year, trigtype, wp))

  return Category.parse_obj(
    {
      'nodetype': 'category',
      'input': "dm",
      'default': -1, # default DM if unrecognized category
      'content': [
        { 'key': dm,
          'value': 
            build_pts(in_file_name(year), in_hist_name(corrtype_dict[corrtype],trigtype,wp,dm_dict_nonmerged[dm]), pt_threshold=getPtThreshold(trigtype))
        } for dm in dms_nonmerged
      ]
    }
  )


def convert_trigger(corrs, year, **kwargs):

  workingpoints = kwargs_get(kwargs, 'workingpoints', default_workingpoints)
  trigger_types = kwargs_get(kwargs, 'triggertypes', trigtypes[year])
  correction_types = kwargs_get(kwargs, 'correctiontypes', corrtypes)
  outdir = kwargs_get(kwargs, 'outdir', 'jsons/')

  """Tau trigger SF, pT- and dm dependent."""
  header("Tau trigger SF, pT- and dm dependent")
  fname = os.path.join(outdir, f'tau_trigger_{idalgo_dict[year]}_{year}.json')
  corr    = Correction.parse_obj({
    'version': 1,
    'name':    "tauTriggerSF",
    'description' : "Tau Trigger SFs and efficiencies for ditau, ditaujet, etau or mutau triggers. " +\
                    "To get the usual DM-specific SF's, specify the DM, otherwise set DM to -1 to get the inclusive SFs. " +\
                    "Default corrections are set to SF's, if you require the input efficiencies, you can specify so in " +\
                    "the corrtype input variable" +\
                    "Note: These SFs are specific for the Htautau CP Analysis (IP significance cuts for DM0 and requiring a refitted SV for DMs 10 & 11)",
    'inputs': [
      {'name': "pt",       'type': "real",   'description': "tau pt"},
      {'name': "dm",       'type': "int",    'description': "tau PNet decay mode (0, 1, 2, 10, or 11, -1)"},
      {'name': "trigtype",       'type': "string",    'description': f"Type of trigger: {trigger_types}"},
      {'name': "wp",       'type': "string", 'description': f"{idalgo_dict[year]}VSjet WP: {workingpoints}"},
      {'name': "corrtype",       'type': "string",    'description': "Type of information: 'eff_data', 'eff_mc', 'sf'"},
      {'name': "syst",     'type': "string", 'description': "systematic 'nom', 'up', 'down'"},
    ],
    'output': {'name': "weight", 'type': "real"},
    'data': { # category:trigtype -> category:corrtype -> category:wp -> category dm -> binning:pt -> category:syst
      'nodetype': 'category', # category:trigtype
      'input': "trigtype",
      'content': [
        { 'key': trigtype,
          'value': {
            'nodetype': 'category', # category:corrtype
            'input': "corrtype",
            'content' : [
              { 'key' : corrtype,
                'value' : {
                  'nodetype': 'category', # category:dm
                  'input': "wp",
                  'content': [
                    { 'key': wp,
                      'value' : build_dms(trigtype, wp, year, corrtype)
                    } for wp in workingpoints
                  ]
                }
              } for corrtype in correction_types
            ]  
          } # category:wp
        } for trigtype in trigger_types
      ]
    } #category:trigtype
  })
  print(f">>> Writing {fname}...")
  with open(fname,'w') as fout:
    JSONEncoder.write(corr,fname,maxlistlen=20)
  corrs.append(corr)


def evaluate(corrs, year, **kwargs):

  workingpoints = kwargs_get(kwargs, 'workingpoints', default_workingpoints)
  trigger_types = kwargs_get(kwargs, 'triggertypes', trigtypes[year])

  header("Evaluate")
  cset = wrap(corrs) # wrap to create C++ object that can be evaluated
  ptbins = [10.,21.,26.,31.,36.,41.,501.,750.,999.,2000.]
  for name in list(cset):
    corr = cset[name]
    print(f">>>\n>>> {name}: {corr.description}")
    for wp in workingpoints:
      print(f">>>\n>>> WP={wp}")
      for dm in dms_nonmerged:
        print(f">>>\n>>> DM={dm}")
        print(">>> %8s"%("trigger type")+" ".join("  %-15.1f"%(p) for p in ptbins))
        for tt in trigger_types:
          row = ">>> %s"%(tt)
          for pt in ptbins:
            sfnom = 0.0
            for syst in ['nom','up','down']:
              #print(">>>   gm=%d, eta=%4.1f, syst=%r sf=%s"%(gm,eta,syst,eval(corr,eta,gm,wp,syst)))
              try:
                sf = corr.evaluate(pt, dm, tt,wp,'sf',syst)
                if 'nom' in syst:
                  row += "%6.2f"%(sf)
                  sfnom = sf
                elif 'up' in syst:
                  row += "%+6.2f"%(sf-sfnom)
                else:
                  row += "%+6.2f"%(sf-sfnom)
              except Exception as err:
                row += "\033[1m\033[91m"+"  ERR".ljust(6)+"\033[0m"
          print(row)
  print(">>>")


def makeRootFiles(corrs, year, **kwargs):

  workingpoints = kwargs_get(kwargs, 'workingpoints', default_workingpoints)
  trigger_types = kwargs_get(kwargs, 'triggertypes', trigtypes[year])
  correction_types = kwargs_get(kwargs, 'correctiontypes', corrtypes)

  cset = wrap(corrs) # wrap to create C++ object that can be evaluated
  out_file = ROOT.TFile('data/tau/correctionHistogramsRebinned_'+str(year)+'.root', 'recreate')
  for name in list(cset):
    corr = cset[name]
    print(f">>>\n>>> {name}: {corr.description}")
    for corrtype in correction_types:
      print(f">>>\n>>> {corrtype}")
      for wp in workingpoints:
        print(f">>>\n>>> WP={wp}")
        for dm in dms_nonmerged:
          print(f">>>\n>>> DM={dm}")
          print(">>> %8s"%("trigger type"))
          for tt in trigger_types:
            print('Building histogram')
            f = uproot.open(in_file_name(year))
            hist = f[in_hist_name(corrtype_dict[corrtype], tt, wp, dm_dict_nonmerged[dm])]
            edges = [round(x,5) for x in hist.to_numpy()[1]]
            ptbin_edges, tmp_values, tmp_errors = merge_pt_bins(edges, hist.values(), hist.errors(), pt_threshold = getPtThreshold(tt))

            hist = ROOT.TH1D('-'.join([corrtype_dict[corrtype], wp, str(dm), tt]), '-'.join([corrtype_dict[corrtype], wp, str(dm), tt]), len(ptbin_edges)-1, numpy.array(ptbin_edges))
            for pt_bin in range(1, len(ptbin_edges)):
              pt_bin_center = ptbin_edges[pt_bin-1]+ (ptbin_edges[pt_bin] - ptbin_edges[pt_bin-1])/2.
              try:
                hist.SetBinContent(pt_bin, corr.evaluate(pt_bin_center, dm, tt, wp, corrtype,'nom'))
                hist.SetBinError(pt_bin, corr.evaluate(pt_bin_center, dm, tt, wp, corrtype,'up')-corr.evaluate(pt_bin_center, dm, tt, wp, corrtype,'nom') )
              except:
                print("Errors for {0} trigger with {1} wp and dm={2} for pt={3}GeV".format(tt, wp, str(dm), str(pt_bin_center)))
            hist.Write()
  out_file.Close()
  print(">>>")


def compareSFs(corrs, year, **kwargs):

  workingpoints = kwargs_get(kwargs, 'workingpoints', default_workingpoints)
  trigger_types = kwargs_get(kwargs, 'triggertypes', trigtypes[year])
  correction_types = kwargs_get(kwargs, 'correctiontypes', corrtypes)

  from TauAnalysisTools.TauTriggerSFs.SFProvider import SFProvider
  cset = wrap(corrs) # wrap to create C++ object that can be evaluated
  for name in list(cset):
    corr = cset[name]
    print(f">>>\n>>> {name}: {corr.description}")
    for wp in workingpoints:
      print(f">>>\n>>> WP={wp}")
      print(">>> %8s"%("trigger type"))
      for tt in trigger_types:
        ptbins = numpy.arange(40., 10000., 0.1) if 'ditau' in tt else numpy.arange(25., 10000., 0.1)
        dms = [0, 1, 2, 10, 11]
        for dm in dms:
          print(f">>>\n>>> DM={dm}")
          old_sfs = SFProvider(in_file_name(year), tt, wp)
          for pt in ptbins:
            if 'sf' in correction_types and abs((old_sfs.getSF(pt, dm, 0) - corr.evaluate(pt, dm, tt,wp,'sf','nom'))/old_sfs.getSF(pt, dm, 0)) > 0.01:
              print("Large difference in SF ({0}) for {1}, {2}, {3}, {4}, {5}".format(str((old_sfs.getSF(pt, dm, 0) - corr.evaluate(pt, dm, tt,wp,'sf','nom'))/old_sfs.getSF(pt, dm, 0)), year, tt, wp, str(dm), str(pt)))
              print("Old: {0} New: {1}".format(old_sfs.getSF(pt, dm, 0), corr.evaluate(pt, dm,tt,wp,'sf','nom')))
            if 'eff_mc' in correction_types and abs((old_sfs.getEfficiencyMC(pt, dm, 0) - corr.evaluate(pt, dm, tt,wp,'eff_mc','nom'))/old_sfs.getEfficiencyMC(pt, dm, 0)) > 0.01:
              print("Large difference in MC EFF ({0}) for {1}, {2}, {3}, {4}, {5}".format(str((old_sfs.getEfficiencyMC(pt, dm, 0) - corr.evaluate(pt, dm, tt,wp,'eff_mc','nom'))/old_sfs.getEfficiencyMC(pt, dm, 0)), year, tt, wp, str(dm), str(pt)))
              print("Old: {0} New: {1}".format(old_sfs.getEfficiencyMC(pt, dm, 0), corr.evaluate(pt, dm,tt,wp,'eff_mc','nom')))
            if 'eff_data' in correction_types and abs((old_sfs.getEfficiencyData(pt, dm, 0) - corr.evaluate(pt, dm, tt,wp,'eff_data','nom'))/old_sfs.getEfficiencyData(pt, dm, 0)) > 0.01:
              print("Large difference in Data EFF ({0}) for {1}, {2}, {3}, {4}, {5}".format(str((old_sfs.getEfficiencyData(pt, dm, 0) - corr.evaluate(pt, dm, tt,wp,'eff_data','nom'))/old_sfs.getEfficiencyData(pt, dm, 0)), year, tt, wp, str(dm), str(pt)))
              print("Old: {0} New: {1}".format(old_sfs.getEfficiencyData(pt, dm, 0), corr.evaluate(pt, dm,tt,wp,'eff_data','nom')))

  print(">>>")


if __name__ == '__main__':
  import argparse
  argParser = argparse.ArgumentParser(description = "Argument parser")
  argParser.add_argument('--years',   action='store', nargs='*', required = True,
                            help='Select years/eras to convert', choices=['2022preEE','2022postEE','2023preBPix','2023postBPix', '2024'])
  argParser.add_argument('--workingpoints',   action='store', nargs='*', default = None, help='Select offline working points to convert', 
                            choices=['VVVLoose', 'VVLoose', 'VLoose', 'Loose', 'Medium', 'Tight', 'VTight', 'VVTight'])
  argParser.add_argument('--triggertypes',   action='store', nargs='*', default = None, help='Select trigger types to convert', 
                            choices=['ditau', 'etau', 'mutau','ditaujet', 'ditauANDditaujet'])
  argParser.add_argument('--correctiontypes',   action='store', nargs='*', default = None, help='Select correction types to convert', 
                            choices=['sf', 'eff_mc', 'eff_data'])
  argParser.add_argument('--outdir',   action='store', default = None, help='Select directory to store output')
  args = argParser.parse_args()

  os.makedirs(args.outdir, exist_ok=True)
  
  for year in args.years:
    idalgo = idalgo_dict.get(year, None)
    corrs = [ ] # list of corrections
    convert_trigger(corrs, year, workingpoints=args.workingpoints, triggertypes=args.triggertypes, correctiontypes=args.correctiontypes, outdir=args.outdir)
    if not idalgo:
      makeRootFiles(corrs, year, workingpoints=args.workingpoints, triggertypes=args.triggertypes, correctiontypes=args.correctiontypes, outdir=args.outdir)
      compareSFs(corrs, year, workingpoints=args.workingpoints, triggertypes=args.triggertypes, correctiontypes=args.correctiontypes, outdir=args.outdir)
  print()
  
