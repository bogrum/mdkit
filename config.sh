#!/usr/bin/env bash
# mdkit — projeye ozel TUM varsayimlar burada. Baska bir veri setiyle calistirmak
# icin bu dosyanin bir kopyasini duzenleyip --config ile verin; kod degismez.

# Veri kokü ve cikti
DATA_ROOT="/mnt/data/scratch-simulations-TUSEB"
RESULTS_DIR="$DATA_ROOT/results"

# Yorumlayicilar. GMX bos birakilirsa PATH'ten bulunur.
GMX="/usr/local/gromacs-2025.4-cuda/bin/gmx"
PYTHON="/home/emre/anaconda3/bin/python"

# Dizin ve dosya adlandirma
COMPLEX_GLOB="*_pandora"
REPS=(rep1 rep2 rep3)
TRAJ_NAME="traj_compact_center_dry.xtc"
TPR_NAME="md_0_10.tpr"
REF_NAME="check_ref.pdb"

# Zincir eslemesi. Analiz scriptleri bu harfleri gormez; index kurulumu
# bunlari kanonik adlara cevirir: RECEPTOR / AUX / LIGAND / RECEPTOR_BB / LIGAND_BB
CHAIN_RECEPTOR="A"   # MHC agir zincir
CHAIN_AUX="B"        # beta-2 mikroglobulin
CHAIN_LIGAND="C"     # peptid
