SKNano.py -a Tutorial_reco_tt -i TTLJ_powheg -n 40 --reduction 10 -e 2024
SKNano.py -a Tutorial_reco_tt -i TTLL_powheg -n 40 --reduction 10 -e 2024
SKNano.py -a Tutorial_reco_tt -i DYto2Mu_MLL50to120 -n 40 --reduction 10 -e 2024
SKNano.py -a Tutorial_reco_tt -i DYto2Mu_MLL10to50  -n 40 --reduction 10 -e 2024
#SKNano.py -a Tutorial_reco_tt -i WtoLNu_1J -n 1 --reduction 1 -e 2024
#SKNano.py -a Tutorial_reco_tt -i WtoLNu_2J -n 1 --reduction 1 -e 2024
#SKNano.py -a Tutorial_reco_tt -i WtoLNu_3J -n 1 --reduction 1 -e 2024
#SKNano.py -a Tutorial_reco_tt -i WtoLNu_4J -n 1 --reduction 1 -e 2024
SKNano.py -a Tutorial_reco_tt -i ST_tW_top_Semilep -n 4 --reduction 10 -e 2024
SKNano.py -a Tutorial_reco_tt -i ST_tW_antitop_Semilep -n 4 --reduction 10 -e 2024
SKNano.py -a Tutorial_reco_tt -i ST_tW_top_Lep -n 4 --reduction 10 -e 2024
SKNano.py -a Tutorial_reco_tt -i ST_tW_antitop_Lep -n 4 --reduction 10 -e 2024
SKNano.py -a Tutorial_reco_tt -i ST_tch_top_lep -n 2 --reduction 10 -e 2024
SKNano.py -a Tutorial_reco_tt -i ST_tch_antitop_lep -n 2 --reduction 10 -e 2024
#SKNano.py -a Tutorial_reco_tt -i WWto2L2Nu -n 10 --reduction 10 -e 2024
#SKNano.py -a Tutorial_reco_tt -i WZto2L2Q -n 10 --reduction 10 -e 2024
#SKNano.py -a Tutorial_reco_tt -i ZZto2L2Q -n 10 --reduction 10 -e 2024
sleep 8m
SKNano.py -a Tutorial_reco_tt -i Muon0 -n 60 --reduction 1 -e 2024
sleep 20m
SKNano.py -a Tutorial_reco_tt -i Muon1 -n 60 --reduction 1 -e 2024
