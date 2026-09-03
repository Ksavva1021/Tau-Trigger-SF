
# Tau-Trigger-SF
A tool for deriving trigger scale factors for the hadronic tau leg, binned not only in $p_\mathrm T$, but also in the PNet decay mode. The steps for setting up and using this tool are outlined below.
___
## 0) Setup

This repository is designed for use inside of CMSSW. I am not sure what version of CMSSW the original author was using, but I am using `CMSSW_14_1_0_pre4` and it seems to work fine.

Once you have set up your installation of CMSSW, clone this repo inside of the `src` directory. Remember to run `cmsenv` before running any of the below steps.
___
## 1) Producing n-tuples from the central NanoAOD files

Triggers have efficiencies which are dependent on the $p_\mathrm T$ of the object being triggered upon. Plotting a trigger's efficiency as a function of $p_\mathrm T$ produces a curve. This curve is referred to as a "turn-on" curve. Producing accurate Monte Carlo simulations of how particles interact with the CMS detector is challenging. As such, the "turn-on" curve in data looks slightly different to how it looks in MC. To make MC look more like real data, we multiply the weight of MC events by a correction factor we call a Scale Factor (which in this case are a function of $p_\mathrm T$ AND the PNet DM of the $\tau_\mathrm h$), which we derive by comparing the "turn-on" curves in data and MC. As detector conditions are era dependent, so too are Scale Factors. Thus, the first step is to procure samples of both data and MC events corresponding to the era for which we are measuring our Scale Factors.

Like HiggsDNA, this tool takes as its input NanoAOD files. The script `scripts/nano_postproc.py` takes these NanoAOD files, drops information about the events we don't need, and restructures them into a columnar format (called n-tuples). Depending on where the NanoAOD files you are trying to access are stored, you made need a valid grid certificate to access them. This involves running something like:

```
voms-proxy-init --rfc --voms cms --valid 192:00 --out ${HOME}/cms.proxy
export X509_USER_PROXY=${HOME}/cms.proxy
```

The script `scripts/nano_postproc.py` is designed to be run twice. Once for the MC events, and another time for the data events.

As we are measuring scale factors specifically for taus, we typically use a DYToTauTau MC sample for deriving the MC turn-on curve. The `--inputJson` flag allows you to pass a sample list identical to those found in HiggsDNA. However, if the `--isMC` flag is also present, it looks specifically for the `DYto2Tau_MLL_50_amcatnloFXFX` NanoAOD files and just uses those. As such, running this step for producing the MC n-tuples looks something like:

```
python3 scripts/nano_postproc.py --inputJson input/2024.json --output output/ntuples_MC_2024 --isMC --era 2024 --batch
```

The `--batch` parallelises the task by submitting jobs with HTCondor. This is recommended, as there are many events, so otherwise will take ages. If you run without the `--isMC` flag (and with `--inputJson` specified) the script will instead search the JSON for the Muon datasets (tau scale factors are most easily measured in the $\tau_\mu \tau_\mathrm h$–channel). This is why I have added all the `HiggsDNA/scripts/ditau/config/Run3_<year>/samples/samples_mt.json` files (which point to copies of the central NanoAOD files stored on the Imperial dCache) to `input/`. Thus, running this step for producing the data n-tuples looks something like:

```
python3 scripts/nano_postproc.py --inputJson input/2024.json --output output/ntuples_data_2024 --era 2024 --batch
```
___
## 2) Skimming the n-tuples
The script `scripts/skimTuple.py` takes n-tuples produced in the previous step, drops the columns we don't need, and adds 4 additional columns: `pass_mutau, pass_etau, pass_ditau, pass_ditaujet`. These correspond to the trigger decision for the 4 monitoring triggers of interest for a given event. All the n-tuples are then merged into a single ROOT file.

Like the previous step, this script is designed to be ran twice, once for data and again for MC. E.g. for 2024 data:

```
python3 scripts/skimTuple.py --input output/ntuples_data_2024/*/*.root --era 2024 --type data --output output/skim_data_2024.root
```

For running the MC skimming, additional ROOT files quantifying the levels of pileup in the data and MC samples need to be provided (someone in the TauPOG is responsible for making these I think). For example:

```
python3 scripts/nano_postproc.py --input_dir output/ntuples_MC_2024 --era 2024 --type mc --pudata pileup/Data_PileUp_2024_69p2.root --pumc pileup/MC_PileUp_2024.root --output output/skim_mc_2024.root
```
___
## 3) Creating the turn-on curves
The script `scripts/createTurnOn_multi.py` computes the trigger efficiencies in different $p_\mathrm T$ and PNet DM bins for data and MC, and outputs the result in a single ROOT file, as well as producing some plots in PDF format for a quick sanity check. It uses the skims produced in the previous step. For example:

```
python3 scripts/createTurnOn_multi.py --input-data output/skim_data_2024.root --input-dy-mc output/skim_mc_2024.root --channels etau,mutau,ditau,ditaujet,ditauANDditaujet --id-algo PNet --working-points VTight --output output/turn_on_2024/TurnOn
```

Unlike in early Run 3 (2022/23), in 2024 the set of events passing the `ditaujet` trigger is not a simple subset of those passing the `ditau` trigger. Hence, the command above includes an additional scale factor calculation corresponding to events which pass both the `ditau` and `ditaujet` triggers.

The choice of `--id-algo` and `--working-points` reflects the current preference of the $\mathrm{H}\rightarrow\tau\tau \ \mathcal{CP}$–working group for analysing late Run 3 data, where we have decided to use the offline selection criterion that hadronic tau candidates must pass the "VTight" working point for "PNetVSjet". For early Run 3 however, "DeepTau2018v2p5VSjet" is used instead. 
___
## 4) Fitting the turn-on curves
The script `scripts/fitTurnOn_multi.py` simply fits a curve to the histograms output by the previous step (and actually calculates the SFs). For example:

```
python3 scripts/fitTurnOn_multi.py --input output/turn_on_2024/TurnOn.root --channels etau,mutau,ditau,ditaujet,ditauANDditaujet --working-points VTight --output jsons/fitTurnOn_2024.root
```
___
## 5) Converting to JSON
The script `scripts/convert_to_json.py` takes the SFs contained in the ROOT file produced by the previous step and puts them into a JSON file in the required `correctionlib` format. The script doesn't allow you to specify an input. Instead it looks for the input ROOT file at `jsons/fitTurnOn_<year>.root`, so ensure your output from the previous step has been moved here and follows the correct naming convention. An example of running this script:

```
python3 scripts/convert_to_json.py --years 2024 --workingpoints VTight --outdir jsons/
```