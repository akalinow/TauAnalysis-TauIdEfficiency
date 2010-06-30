#!/bin/bash

# Run all of the relevant TaNC eff plots and distributions

##### PLOTS FOR THE PAS #####
python tanc_fakerates_for_pas.py &

# Wait for completion
wait

############################
##### PLOTS FOR THE AN #####
############################

# Control plots
# (assume someone else will make these)
#python shrinkingConePlots.py -p jetPt 
#python shrinkingConePlots.py -p jetEta 
#python shrinkingConePlots.py -p jetPhi 

# Fake rate plots
python shrinkingConePlots.py -p TaNCQuarter -v Pt -v Eta -v Phi &
python shrinkingConePlots.py -p TaNCHalf -v Pt -v Eta -v Phi &
python shrinkingConePlots.py -p TaNCOne -v Pt -v Eta -v Phi &
wait

##### Resolution plots #####
python make_tanc_resolution_plot.py &

##### ZTT efficiency plots #####
# NB this makes plots for ZTT for all algos
# I think Abdollah is running these
#python plot_Ztautau_v3.py &


# finish up
wait



