
# Tau-Trigger-SF
A tool for deriving trigger scale factors for the hadronic tau leg, binned not only in $p_\mathrm T$, but also in the PNet decay mode. The steps for setting up and using this tool are outlined below.
___
## 0) Setup

This repository is designed for use inside of CMSSW. I am not sure what version of CMSSW the original author was using, but I am using `CMSSW_14_1_0_pre4` and it seems to work fine.

Once you have set up your installation of CMSSW, clone this repo inside of the `src` directory. Remember to run `cmsenv` before running any of the below steps.
___
## 1) Skimming the central NanoAOD files

Triggers have efficiencies which are dependent on the $p_\mathrm T$ of the object being triggered upon. Plotting a trigger's efficiency as a function of $p_\mathrm T$ produces a curve. This curve is referred to as a "turn-on" curve. Producing accurate Monte Carlo simulations of how particles interact with the CMS detector is challenging. As such, the "turn-on" curve in data looks slightly different to how it looks in MC. To make MC look more like real data, we multiply the weight of MC events by a correction factor we call a Scale Factor (which is usually a function of $p_\mathrm T$), which we derive by comparing the "turn-on" curves in data and MC. As detector conditions are era dependent, so too are Scale Factors. Thus, the first step is to procure samples of both data and MC events corresponding to the era for which we are measuring our Scale Factors.

Like HiggsDNA, this tool takes as its input NanoAOD files. The script `scripts/nano_postproc.py` takes these NanoAOD files, drops information about the events we don't need, and restructures them into a columnar format (called `ntuples`). Depending on where the NanoAOD files you are trying to access are stored, you made need a valid grid certificate to access them. This involves running something like:

```
voms-proxy-init --rfc --voms cms --valid 192:00 --out ${HOME}/cms.proxy
export X509_USER_PROXY=${HOME}/cms.proxy
```

The script `scripts/nano_postproc.py` is designed to be run twice. Once for the MC events, and another time for the data events.

As we are measuring scale factors specifically for taus, we typically use a DYToTauTau MC sample for deriving the MC turn-on curve. The `--inputJson` flag allows you to pass a sample list identical to those found in HiggsDNA. However, if the `--isMC` flag is also present, it looks specifically for just the `DYto2Tau_MLL_50_amcatnloFXFX` NanoAOD files and just uses those. As such, running this step for producing the MC ROOT file looks something like:

```
python3 scripts/nano_postproc.py --inputJson input/2024.json --output output/ntuples_MC_2024 --isMC --era 2024 --batch
```

The `--batch` parallelises the task by submitting jobs with HTCondor. This is recommended, as there are many events, so otherwise will take ages. If you run without the `--isMC` flag (and with `--inputJson` specified) the script will instead search the JSON for the Muon datasets (tau scale factors are most easily measured in the $\tau_\mu \tau_\mathrm h$–channel). This is why I have added all the `HiggsDNA/scripts/ditau/config/*/samples/samples_mt.json` files to `input/`. Thus, running this step for producing the data ROOT file looks something like:

```
python3 scripts/nano_postproc.py --inputJson input/2024.json --output output/ntuples_data_2024 --era 2024 --batch
```