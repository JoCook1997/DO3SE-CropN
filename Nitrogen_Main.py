#Main code for the DO3SE model Nitrogen module
#Jo Cook
#date created 08/07/2022
#last edited 26/03/2024

#import the functions and parameters for the nitrogen module
import pandas as pd
import json
import warnings
from datetime import datetime

from Nitrogen_Functions import *

start_time=datetime.now()

for file_id in file_ids:

    from Nitrogen_Parameters import * #this is so it re-initialises for every file_id

    #read in the DO3SE output file
    DO3SE_Output=read_DO3SE_output(file_path+"DO3SE_Outputs"+calib_eval, file_id)

    #get partition fractions used in model run
    with open (file_path+"DO3SE_Outputs"+calib_eval+file_id+"/processed_config.json","r") as config:
        DO3SE_config=json.load(config)

    #counters
    row_counter=0   
    counter=0
    
    #collect the data to calculate the new min stem
    #and leaf N conc.'s as affected by the average O3 conc.
    M12O3=avgM12O3(DO3SE_Output)
    LeafGamma=DO3SE_config["carbon_allocation"]["gamma"]
    LeafDelta=DO3SE_config["carbon_allocation"]["delta"]
    #calculates minimum leaf and stem N
    #this is dependent on average M12 value
    leafN_min=leafNconc_O3(DO3SE_Output["canopy_lai"].iloc[-1], DO3SE_Output["canopy_lai_brown"].iloc[-1], M12O3)
    stemN_min=stemNconc_O3(M12O3)

    #get the critical ozone flux for antioxidant function
    cL3=DO3SE_config["Land_Cover"]["parameters"][0]["pn_gsto"]["cL3"]

    #get the ratio of mass in grain versus ear
    grain_to_ear=DO3SE_config["carbon_allocation"]["grain_to_ear"]

    while row_counter < len(DO3SE_Output):
        
        #set dvi, lai and current fst_acc
        dvi=DO3SE_Output.loc[row_counter,'dvi']
        LAI=DO3SE_Output.loc[row_counter,'canopy_lai']
        fstAcc=DO3SE_Output.loc[row_counter,'fst_acc']

        if dvi>0 and dvi<=2: #thermal time is after emergence

            #get growth/decrease of stem, and increase/ decrease of LAI
            weightStem, growthStem, decreaseStem, growthLAI, decreaseLAI=stem_leaf_growth(DO3SE_Output.loc[counter,'c_stem'], 
                                                DO3SE_Output.loc[counter,'c_resv'], weightStem_yesterday, LAI_yesterday, LAI)

            if dvi>0 and dvi <1: #before grain starts to fill
                
                #calculate required N
                Nup_preanth=req_N_uptake(growthStem, growthLAI,weightStem,weightStem_yesterday,stemN,Nup_premax)
                Nup_postanth=0

                #distribute N to stem and leaf
                into_Nharv, intoStem, leavingStem, intoLeaf, leavingLeaf=distribute_N_uptake_pre_anth(weightStem,Nup_preanth,stemN,\
                                                                                                      growthLAI,decreaseLAI,stemN_min,leafN_min)
                
                #update N pools
                stemN=stemN+intoStem-leavingStem
                leafN=leafN+intoLeaf-leavingLeaf
                
            elif dvi >=1: #grain is filling
                #calculate thermal times for post-anth Nuptk eqn
                td_since_anth=DO3SE_Output.loc[row_counter,"td"] \
                                - DO3SE_config["Land_Cover"]["parameters"][0]["phenology"]["key_dates_td"]["Astart"]
                grainfill_td=DO3SE_config["Land_Cover"]["parameters"][0]["phenology"]["key_dates_td"]["harvest"] \
                                - DO3SE_config["Land_Cover"]["parameters"][0]["phenology"]["key_dates_td"]["Astart"]

                #N uptake when grain is filling
                Nup_postanth=post_anth_N_uptk(stemN,weightStem,grainfill_td,td_since_anth)
                Nup_preanth=0
                
                #distribute N from leaf and stem and uptake
                into_Nharv, intoLeaf, intoStem, leavingLeaf, leavingStem, p_harv, increaseLeafAntioxidant, increaseStemAntioxidant \
                    = distribute_N_post_anth(dvi,leafN, stemN, decreaseLAI, weightStem,  \
                                             Nup_postanth,leafN_min,stemN_min, \
                                             cL3, fstAcc)
               
                #update N pools
                stemN=stemN+intoStem-leavingStem
                leafN=leafN+intoLeaf-leavingLeaf
                leafAntioxidants=leafAntioxidants+increaseLeafAntioxidant
                stemAntioxidants=stemAntioxidants+increaseStemAntioxidant
                
                #get reduction in LAI
                growthLAI, decreaseLAI=change_LAI(LAI,LAI_yesterday)
                
        elif dvi<=0 or dvi>2: #if thermal time greater than end of seed growth or before emergence...
            Nup_preanth,Nup_postanth,leavingLeaf,leavingStem,intoLeaf,intoStem,into_Nharv,decreaseLAI, \
                growthLAI,p_harv, increaseLeafAntioxidant, increaseStemAntioxidant=set_not_growing_params()
            
            #update N pools
            stemN=stemN+intoStem-leavingStem
            leafN=leafN+intoLeaf-leavingLeaf
            leafAntioxidants=leafAntioxidants+increaseLeafAntioxidant
            stemAntioxidants=stemAntioxidants+increaseStemAntioxidant
        
        #update the grain N, and cumulative N uptake
        Nharv=Nharv+into_Nharv
        CumNup=CumNup+Nup_preanth+Nup_postanth
        if dvi<=1: CumNup_pre_anth=CumNup
        #save some of the parameters
        for counter in range(counter, counter+24):
            DO3SE_Output.loc[counter,'stemN']=stemN+stemAntioxidants #I seperate usable N from antioxidant N for ease of modelling
            DO3SE_Output.loc[counter,'leafN']=leafN+leafAntioxidants  #so I need to add it back to get correct measure of N
            DO3SE_Output.loc[counter,'decreaseLAI']=decreaseLAI
            DO3SE_Output.loc[counter,'Nharv']=Nharv
            DO3SE_Output.loc[counter,'N_grain']=Nharv*ear_grain_N
            DO3SE_Output.loc[counter,'N_chaff']=Nharv*(1-ear_grain_N)
            DO3SE_Output.loc[counter,'CumNup']=CumNup
            DO3SE_Output.loc[counter,'p_nharv']=p_harv
            DO3SE_Output.loc[counter,'NleavingLeaf']=leavingLeaf
            DO3SE_Output.loc[counter,'NleavingStem']=leavingStem
            DO3SE_Output.loc[counter,'LAI']=LAI
            DO3SE_Output.loc[counter,'growthLAI']=growthLAI
            DO3SE_Output.loc[counter,'NintoStem']=intoStem
            DO3SE_Output.loc[counter,'NintoLeaf']=intoLeaf
            DO3SE_Output.loc[counter,'LeafAntioxidantN']=leafAntioxidants
            DO3SE_Output.loc[counter,'StemAntioxidantN']=stemAntioxidants
            DO3SE_Output.loc[counter,'LeafNonAntioxidantN']=leafN
            DO3SE_Output.loc[counter,'StemNonAntioxidantN']=stemN

            with warnings.catch_warnings(): 
                warnings.filterwarnings('error')
                #N concentrations as a percentage
                try: DO3SE_Output.loc[counter,'grain_N_conc']=100*Nharv*ear_grain_N/(DO3SE_Output.loc[counter,'c_harv']*kg_to_g*grain_to_ear/frac_carbon) #gNm-2/gDMm-2 converted to a %
                except: DO3SE_Output.loc[counter,'grain_N_conc'] = 0
                try: DO3SE_Output.loc[counter,'chaff_N_conc']=100*Nharv*(1-ear_grain_N)/(DO3SE_Output.loc[counter,'c_harv']*kg_to_g*(1-grain_to_ear)/frac_carbon) #gNm-2/gDMm-2 converted to a %
                except: DO3SE_Output.loc[counter,'chaff_N_conc'] = 0
                try: DO3SE_Output.loc[counter,'stem_N_conc']=100*(stemN+stemAntioxidants)/((DO3SE_Output.loc[counter,'c_stem']+DO3SE_Output.loc[counter,'c_resv'])*kg_to_g/frac_carbon) #gNm-2/gDMm-2 converted to a %
                except: DO3SE_Output.loc[counter,'stem_N_conc']=0 
                DO3SE_Output.loc[counter,'leaf_N_conc']=100*get_leafN_conc(DO3SE_Output.loc[counter,"leaf_dm"],DO3SE_Output.loc[counter,"lbrn_dm"],
                                                                            DO3SE_Output.loc[counter,"leafN"]) #gNm-2/gDMm-2 converted to a %
            counter+=1
        
        LAI_yesterday=LAI
        weightStem_yesterday=weightStem
        
        #go to next day
        row_counter+=24

    DO3SE_Output.to_csv(file_path+'/N_Outputs'+calib_eval+'/'+file_id+'/'+file_id+'_N.csv', index=False)
    anthind=DO3SE_Output[DO3SE_Output["dvi"] >= 1].index[0]
    stemanthconc=DO3SE_Output.at[DO3SE_Output.index[anthind],'stem_N_conc']
    leafanthconc=DO3SE_Output.at[DO3SE_Output.index[anthind],'leaf_N_conc']

    #make_plots(DO3SE_Output, DO3SE_config,file_id)
end_time=datetime.now()
print(end_time-start_time)