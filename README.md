
# Tau-Trigger-SF
A tool for deriving trigger scale factors for the hadronic tau leg, binned not only in $p_\mathrm T$, but also in the PNet decay mode. The steps for setting up and using this tool are outlined below.
___
## 0) Setup

This repository is designed for use inside of CMSSW. I am not sure what version of CMSSW the original author was using, but I am using `CMSSW_14_1_0_pre4` and it seems to work fine.

Once you have set up your installation of CMSSW, clone this repo inside of the `src` directory. Remember to run `cmsenv` before running any of the below steps.
___
## 1) Skimming the central NanoAOD files

Triggers have efficiencies which are dependent on the $p_\mathrm T$ of the object being triggered upon. Plotting a trigger's efficiency as a function of $p_\mathrm T$ produces a curve. This curve is referred to as a "turn-on" curve. Producing accurate Monte Carlo simulations of how particles interact with the CMS detector is challenging. As such, the "turn-on" curve in data looks slightly different to how it looks in MC. To make MC look more like real data, we multiply the weight of MC events by a correction factor we call a Scale Factor (which is usually a function of $p_\mathrm T$), which we derive by comparing the "turn-on" curves in data and MC. As detector conditions are era dependent, so too are Scale Factors. Thus, the first step is to procure samples of both data and MC events corresponding to the era for which we are measuring our Scale Factors.

Like HiggsDNA, 



