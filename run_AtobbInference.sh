#!/bin/bash

# Background samples
#SKNano.py -a AtobbInference -i TTLJ_powheg -n 200 -e 2024
#SKNano.py -a AtobbInference -i TTLL_powheg -n 200 -e 2024
#sleep 2h
#SKNano.py -a AtobbInference -i DYto2Mu_MLL50to120 -n 100 -e 2024
#SKNano.py -a AtobbInference -i DYto2E_MLL50to120 -n 100 -e 2024
#
#SKNano.py -a AtobbInference -i ST_tW_top_Semilep -n 10 -e 2024
#SKNano.py -a AtobbInference -i ST_tW_antitop_Semilep -n 10 -e 2024
#SKNano.py -a AtobbInference -i ST_tW_top_Lep -n 10 -e 2024
#SKNano.py -a AtobbInference -i ST_tW_antitop_Lep -n 10 -e 2024
#SKNano.py -a AtobbInference -i ST_tch_top_lep -n 5 -e 2024
#SKNano.py -a AtobbInference -i ST_tch_antitop_lep -n 5 -e 2024
#sleep 1h
#
## Data samples
#SKNano.py -a AtobbInference -i Muon0 -n 20 -e 2024
#SKNano.py -a AtobbInference -i Muon1 -n 20 -e 2024
#SKNano.py -a AtobbInference -i EGamma0 -n 20 -e 2024
#SKNano.py -a AtobbInference -i EGamma1 -n 20 -e 2024

# Signal samples (MHc130, MA18-125)
SKNano.py -a AtobbInference -i TTToHcToWAToBB-MHc130_MA18_SingleLepFilter -n 5 -e 2024
SKNano.py -a AtobbInference -i TTToHcToWAToBB-MHc130_MA25_SingleLepFilter -n 5 -e 2024
SKNano.py -a AtobbInference -i TTToHcToWAToBB-MHc130_MA36_SingleLepFilter -n 5 -e 2024
SKNano.py -a AtobbInference -i TTToHcToWAToBB-MHc130_MA50_SingleLepFilter -n 5 -e 2024
SKNano.py -a AtobbInference -i TTToHcToWAToBB-MHc130_MA70_SingleLepFilter -n 5 -e 2024
SKNano.py -a AtobbInference -i TTToHcToWAToBB-MHc130_MA90_SingleLepFilter -n 5 -e 2024
SKNano.py -a AtobbInference -i TTToHcToWAToBB-MHc130_MA110_SingleLepFilter -n 5 -e 2024
SKNano.py -a AtobbInference -i TTToHcToWAToBB-MHc130_MA125_SingleLepFilter -n 5 -e 2024
