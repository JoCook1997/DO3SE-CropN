#Functions used by the DO3SE model Nitrogen module
#Jo Cook
#date created 08/07/2022
#last edited 21/08/2023
#--------------------------------------------------------------------------------------------
import pandas as pd
import math
import numpy as np
import matplotlib.pyplot as plt

from Nitrogen_Parameters import *

#when the crop is not growing these parameters should be zero
def set_not_growing_params():
    N_uptake_pre_anth=0
    N_uptake_post_anth=0
    leaving_leaf=0
    leaving_stem=0
    into_leaf=0
    into_stem=0
    into_grain=0
    decrease_LAI=0
    growth_LAI=0
    p_nharv=0
    leafAntioxidantIncrease=0
    stemAntioxidantIncrease=0

    return N_uptake_pre_anth,N_uptake_post_anth, \
            leaving_leaf, leaving_stem, into_leaf, into_stem, \
                into_grain, decrease_LAI, growth_LAI, p_nharv, \
                    leafAntioxidantIncrease, stemAntioxidantIncrease

#-------------------------------------------------------------------------------------------

#get the N uptake from the soil post anthesis using Sirius Quality
#Waterlogging and high temperature stress can be modelled using
#Pan et al (2006) with the modifications from Osman et al. (2021)
#but are not included here

def post_anth_N_uptk(N_stem,stem_weight,T_grainfill,T_current):
    
    #potential uptake post anthesis from Sirius Quality is designed to
    #go to 0 as grain filling approaches it's end
    #in SiriusQuality they define a max stem N to get the capacity
    #I don't do that, I just allow it to be variable and usually use
    #the stem as like a storage bin for extra N

    N_pot=Nup_postmax*(T_grainfill-T_current)/T_grainfill 
    
    #added this to Sirius equations
    #calculates how different current stem N is from the target
    diff_target_capacity=(stem_weight*stemN_targ)-N_stem
    
    if diff_target_capacity<=0: #stem is holding more N than target
        post_anth_uptk=0
    else: #stem is not holding more than the target N, use Sirius equation
        post_anth_uptk=np.minimum(N_pot,diff_target_capacity)

    if post_anth_uptk<0: #limit post anth uptk so it can't be less than zero
        post_anth_uptk=0
    
    return post_anth_uptk

#-------------------------------------------------------------------------------------------

#determine the decrease or increase in leaf area index
def change_LAI(leaf_area_ind,leaf_area_ind_yest):
    growth_leaf_area_ind=leaf_area_ind-leaf_area_ind_yest
    if growth_leaf_area_ind<0:
        decrease_leaf_area_ind=growth_leaf_area_ind*-1
        growth_leaf_area_ind=0
    elif growth_leaf_area_ind==0:
        decrease_leaf_area_ind=0
    else:
        decrease_leaf_area_ind=0
    return growth_leaf_area_ind, decrease_leaf_area_ind

#-------------------------------------------------------------------------------------------
#get the growth in the stem and leaf area
def stem_leaf_growth(stem_carbon, resv_carbon, stem_weight_yest, LAI_yest, LAI_current):
    stem_weight=(stem_carbon+resv_carbon)*kg_to_g/frac_carbon
    stem_growth=stem_weight-stem_weight_yest
    stem_decrease=0
    if stem_growth<0:
        stem_growth=0
        stem_decrease=stem_growth*-1

    #get increase/ decrease in leaf area index compared with prev. day
    growth_LAI, decrease_LAI=change_LAI(LAI_current, LAI_yest)
    
    return stem_weight, stem_growth, stem_decrease, growth_LAI, decrease_LAI

#-------------------------------------------------------------------------------------------

#function to get the N demand
def req_N_uptake(stem_growth, leaf_area_growth, stem_weight,stem_weight_yest, N_stem, max_uptk):
    #get the N deficit of the stem  from yesterday
    N_deficit=(stem_weight_yest*stemN_targ)-N_stem
    if N_deficit<0: 
        N_deficit=0
        
    #get uptake incl. N deficit
    req_uptake=(stem_growth*stemN_targ)+(leaf_area_growth*leafN_targ)+N_deficit
    if req_uptake>max_uptk: 
        req_uptake=max_uptk
    
    if req_uptake<0: 
        req_uptake=0 #Would this ever happen?
    
    return req_uptake

#-------------------------------------------------------------------------------------------

def distribute_N_uptake_pre_anth(stem_weight,N_uptk,N_stem,leaf_area_growth,leaf_area_decrease,\
                                 stemN_senes,leafN_senes):
    
    #pre-anthesis no N goes to grain
    into_grain=0
    
    #if there is less N in the stem than minimum allowed stem concentration
    if N_stem<=stem_weight*stemN_senes:
        #get N to bring stem N conc. to it's minimum
        into_stem=(stem_weight*stemN_senes)-N_stem
        leaving_stem=0
        #if there's not enough N to meet stem growth at minimum N
        if into_stem>N_uptk:
            #code to senesce LAI to make up for N deficit should go here, however
            #in the current model leaf N and LAI development are not linked
            #so I will not allow release of leafN unless senescence is occuring
            #or the LAI in DO3SE is decreasing 14/06/2023
            into_leaf=0
            leaving_leaf=0
            into_stem=N_uptk
            
        #there is enough N to continue stem growth at it's minimum N conc.
        else:
            into_leaf=leaf_area_growth*leafN_targ
            #if there isn't enough N leftover to meet leaf demand reduce it
            if into_leaf>N_uptk-into_stem: 
                into_leaf=N_uptk-into_stem
            into_stem=N_uptk-into_leaf
            leaving_leaf=0
    #stem N is not below it's minimum
    else:
        into_leaf=leaf_area_growth*leafN_targ
        leaving_leaf=0
        #if leaf requires more N than taken up, fulfill excess demand using stem N
        if into_leaf>=N_uptk:
            into_stem=0
            leaving_stem=into_leaf-N_uptk
            #only reduce stem N to it's minimum N conc. to fulfill leaf demand
            if leaving_stem>N_stem-(stem_weight*stemN_senes):
                leaving_stem=N_stem-(stem_weight*stemN_senes)
            into_leaf=N_uptk+leaving_stem 
        #if there is enough uptake to fulfill leaf demand transfer excess to stem for storage
        else:
            into_stem=N_uptk-into_leaf
            leaving_stem=0

    #lai can senesce before anthesis so release the leaf N to the stem
    if leaf_area_decrease>0:
        leaving_leaf=leaving_leaf+(leaf_area_decrease*(leafN_targ-leafN_senes))
        into_stem=into_stem+(leaf_area_decrease*(leafN_targ-leafN_senes))
    return into_grain, into_stem, leaving_stem, into_leaf, leaving_leaf

#-------------------------------------------------------------------------------------------

def distribute_N_post_anth(development_index, N_leaf, N_stem, decrease_LAI, \
                           stem_weight,post_anth_uptk,leafN_senes,stemN_senes, \
                            threshold_O3flux, current_O3flux):
    #fraction N to grain
    p_nharv=N_to_grain(development_index)
    
    #N released from leaf area senescing
    released_leafN=decrease_LAI*(leafN_targ-leafN_senes) #because of this structure, 
                                                                           #N can leave leaf when it didn't accumulate N
                                                                            #as decrease in LAI is associated with carbon mass 
                                                                            #and I calculate equivalent released N for that decrease 
                                                                            #in leaf LAI
    
    #here I subtract all released N, regardless of whether it is later used for antioxidants and stays in leaf
    #Nleaf therefore becomes a measure of usable leaf N, same with stemN
    N_leaf=N_leaf-released_leafN

    #calculate available nitrogen for the grain
    available_stemN=N_stem-(stem_weight*stemN_senes)

    if Antioxidants=="True":
        #redefine released_leafN since antioxidants mean some of that will stay in leaf
        #not move to stem. But don't add it back to leaf pool. Keep leaf and stem
        #antioxidants seperate so that they don't accidentally get mobilised to grain
        #total leaf N can be calculated before writing to output file
        leafAntioxidantN, stemAntioxidantN=antioxidant_effect(released_leafN,available_stemN,threshold_O3flux,current_O3flux)
    else:
        leafAntioxidantN=0
        stemAntioxidantN=0
    
    #amount of N that grain should take
    into_grain=(available_stemN-stemAntioxidantN)*p_nharv

    #N released from senescing LAI and taken up post anthesis is moved to stem as it would have to
    #pass through the stem in order to get to the grain
    #update the N pool for the stem
    N_stem=N_stem+(released_leafN-leafAntioxidantN)-stemAntioxidantN+post_anth_uptk-into_grain

    #record the movement of Nitrogen to be used to update the leaf and stem pools in the main code
    #leaf and stem N in this function are local not global variables
    into_leaf=0
    leaving_leaf=released_leafN
    into_stem=(released_leafN-leafAntioxidantN)+post_anth_uptk
    leaving_stem=into_grain+stemAntioxidantN
    
    return into_grain, into_leaf, into_stem, leaving_leaf, leaving_stem,p_nharv, leafAntioxidantN, stemAntioxidantN

#-------------------------------------------------------------------------------------------

#define the fraction of labile N that goes to the grain as a function of DVI
def N_to_grain(development_index):
    frac_to_grain=1/(1+math.exp(-alpha*(development_index-beta))) #sigmoid
    return frac_to_grain

#-------------------------------------------------------------------------------------------
#calculate the overall leaf N concentration
#at some point could also do a green and brown leaf N conc.
def get_leafN_conc(greenleafDM,brownleafDM,total_leafN):
    if (brownleafDM+greenleafDM)>0:
        leafN_conc=total_leafN/(brownleafDM+greenleafDM)
    else:
        leafN_conc=0
    return leafN_conc
#-------------------------------------------------------------------------------------------
#calculate the minimum leaf N post-anthesis under O3 exposure
#O3 inhibits N remobilisation so minimum leaf N concentration
#increases with O3 exposure
def leafNconc_O3(final_greenLAI,final_brownLAI,avgO3):
    #see %RemainingLeafStem .xlsx file for info 
    pcnt_to_frac=1/100
    
    leafN_minconc=pcnt_to_frac*((leaf_O3gradient*avgO3)+leaf_O3intercept)
                    
    #because not all LAI senesces, achieving this leaf N conc. as the final value means
    #that the N concentration in brown leaf area needs to be lower
    #this equation uses the average leaf N conc at harvest and compares it with leafN_minconc to calculate
    #what the brown leaf min N conc needs to be 
    O3_senes_leafN_minconc=((leafN_minconc*(final_greenLAI+final_brownLAI))-(leafN_targ*final_greenLAI))/final_brownLAI #problem here for HD cv
    
    if O3_senes_leafN_minconc<0:
        O3_senes_leafN_minconc=0
    return O3_senes_leafN_minconc

#-------------------------------------------------------------------------------------------
#calculate the minimum stem N post-anthesis under O3 exposure
#O3 inhibits N remobilisation so minimum stem N concentration
#increases with O3 exposure
def stemNconc_O3(avgO3):
    #as above, calculating a new minimum stem N concentration
    pcnt_to_frac=1/100
    O3_senes_stemN_minconc=pcnt_to_frac*((stem_O3gradient*avgO3)+stem_O3intercept)
    return O3_senes_stemN_minconc
#-------------------------------------------------------------------------------------------
#calculate antioxidants either directly, or via a modified N:protein conversion factor
#takes inputs of N that has left the leaf and stem that day
def antioxidant_effect(releasedLeafN,releasedStemN,cL3,fst_acc):

    #unit convert. fst_acc is in nmol but cL3 is in micro
    fst_acc=fst_acc/1000

    #initialise variables
    antioxidantleafN=0
    antioxidantstemN=0
    if Antioxidants=="False":
        pass
    elif Antioxidants=="True": #if user has chosen to deal with antioxidant effect then proceed
            if fst_acc>cL3: #if fst has exceeded damage threshold
                O3stressfactorLeaf=((fst_acc-cL3)/((gradient_modifier_leaf*endPoint)-cL3))
                O3stressfactorStem=((fst_acc-cL3)/((gradient_modifier_stem*endPoint)-cL3))
                if O3stressfactorLeaf>1: 
                    O3stressfactorLeaf=1
                if O3stressfactorStem>1: 
                    O3stressfactorStem=1
                if O3stressfactorLeaf<0: 
                    O3stressfactorLeaf=0
                if O3stressfactorStem<0: 
                    O3stressfactorStem=0
                antioxidantleafN=O3stressfactorLeaf*releasedLeafN
                antioxidantstemN=O3stressfactorStem*releasedStemN
    else:
        print("Invalid choice for dealing with antioxidant response. Type 'False' in config for no inclusion of this feature and 'True' to include it")
    
    return antioxidantleafN, antioxidantstemN
    
#-------------------------------------------------------------------------------------------
def avgM12O3(df):
    subset1=df[df["hr"].between(9,20)] #selects 9am-9pm (look at an output file to understand 20:00 is 8pm goes up to 9pm)
    subset2=subset1[subset1['o3_ppb'] > 20] #20ppb is our gap filling pre O3 so 
                                            #select >20ppb so we hopefully capture
                                            #avg O3 after it has been turned on
    M12=subset2.loc[:,'o3_ppb'].mean()
    return M12
#-------------------------------------------------------------------------------------------
#can make this handle a list of input files eventually if necessary
def read_DO3SE_output(path, i_d):
    read_output_data=pd.DataFrame(pd.read_csv(str(path+'/'+i_d+'/'+i_d+'_out.csv')))
    return read_output_data
#-------------------------------------------------------------------------------------------