#!/usr/bin/env bash
 
#get port from arg, or use default 8002
port=${1:-8002}
 
#source ~/mambaforge/etc/profile.d/conda.sh
source $(conda info --base)/etc/profile.d/conda.sh
conda activate depdash_env
 
cmd='python squlite_visualization.py -p '$port
echo run: $cmd
echo
eval $cmd
#python dash/app.py -p 8003
