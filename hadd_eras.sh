#!/bin/bash

set -e  # stop on error

OUTDIR="skimTuples/All_combined"
mkdir -p $OUTDIR

echo "Hadding data files..."
hadd -f $OUTDIR/skim_data.root skimTuples/Run3_2022/skim_data.root skimTuples/Run3_2022EE/skim_data.root skimTuples/Run3_2023/skim_data.root skimTuples/Run3_2023BPix/skim_data.root

echo "Hadding mcNLO files..."
hadd -f $OUTDIR/skim_mcNLO.root skimTuples/Run3_2022/skim_mcNLO.root skimTuples/Run3_2022EE/skim_mcNLO.root skimTuples/Run3_2023/skim_mcNLO.root skimTuples/Run3_2023BPix/skim_mcNLO.root

echo "All done. Output in $OUTDIR/"

