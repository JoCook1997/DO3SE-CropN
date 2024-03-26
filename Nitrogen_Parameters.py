#Parameters required by the DO3SE model Nitrogen module
#Jo Cook
#date created 08/07/2022
#last edited 01/03/2023

#---------------------------------------------------------------------------------------------
import json
#---------------------------------------------------------------------------------------------
# get config with N parameters in it
with open(r'configN.json') as file:
    configN = json.load(file)

#file path and ID variables
file_path=configN["runInfo"]["file_path"]
calib_eval=configN["runInfo"]["calib_eval"]
file_ids=configN["runInfo"]["file_ids"]

#---------------------------------------------------------------------------------------------
#Constant variables

#carbon to dry matter conversions
kg_to_g=1000 #conversion factor
frac_carbon=0.5 #carbon fraction of dry biomass is 0.5

#target N
leafN_targ=configN["concN"]["leafN_targ"] #target leaf N in green leaves 
stemN_targ=configN["concN"]["stemN_targ"] #target stem N conc.

#O3 remobilisation values
leaf_O3gradient=configN["O3Remob"]["leaf_O3gradient"]
leaf_O3intercept=configN["O3Remob"]["leaf_O3intercept"]
stem_O3gradient=configN["O3Remob"]["stem_O3gradient"]
stem_O3intercept=configN["O3Remob"]["stem_O3intercept"]

#parameters describing uptake
Nup_premax=configN["uptakeN"]["Nup_premax"] #preanthesis maximum uptake of nitrogen due to leaf area expansion g N m-2 day-1
Nup_postmax=configN["uptakeN"]["Nup_postmax"] #postanthesis maximum uptake of nitrogen g N m-2 day-1

#parameter to account for the fact that N pool is for chaff and grain
ear_grain_N=configN["transfertograinN"]["ear_grain_N"] #could calibrate this, look in EarGrainN.ipynb for plot of ranges, about 0.75-0.98

#parameters controlling sigmoid grain filling function
alpha=configN["transfertograinN"]["alpha"]
beta=configN["transfertograinN"]["beta"]

#---------------------------------------------------------------------------------------------
#initialise variables

NUP=0 #initial N requirement from uptake is 0
Nup_preanth=0 #no N uptake
Nup_postanth=0
CumNup_pre_anth=0
Nharv=0 #current accumulated N in the harvest pool is 0 g m-2
weightStem=0
stemN=weightStem*stemN_targ #N content of the stem in g m-2
leafN=0
CumNup=stemN+leafN
weightStem_yesterday=0
LAI_yesterday=0
p_harv=0 #no partitioning of N to grain until grain starts to fill
hdd_sum=0


