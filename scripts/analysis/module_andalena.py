#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 11:21:47 2023

@author: magdalena
"""

import csv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import string
from scipy import stats
from scipy.odr import *
from itertools import product
from sklearn.metrics import r2_score


plt.rcParams["font.family"] = "Arial"

def plate_OD_time_t(file_name):
    #returns OD data from excel exported from plate reader
    df = pd.read_excel(file_name)  
    OD_data = df.iloc[9:17,1:13]
    return(OD_data)

def add_raw_data_to_dic(d_replicates_filename):
    #create dictionary with the raw data of each replicate
    
    #create empty dictionary
    d_raw_data = {}
    for r in d_replicates_filename:
        d_raw_data[r] = {}

    #get data of replicate r
    plate_numb = ['1','2','3']
    for r in d_replicates_filename:
        for plate in plate_numb:
            d_raw_data[r][plate] = {}
            file_name = d_replicates_filename[r][int(plate)-1]

            OD_data = plate_OD_time_t(file_name) 

            d_raw_data[r][plate] = OD_data
                
    #
    return(d_raw_data)

def get_processed_data(d_replicates_filename, d_blank):
    #get raw data
    d_raw_data = add_raw_data_to_dic(d_replicates_filename)
    
    #now process it
    d_data = {}
    for k in d_raw_data:
        d_data[k] = {}

    row_plate = string.ascii_uppercase[:8]
    column_plate = range(1,13)
    
    for rep in d_raw_data:
        for plate in d_raw_data[rep]:
            for ic,c in enumerate(d_raw_data[rep][plate]):
                str_c = str(column_plate[ic])
                if len(str_c)==1:
                    str_c = '0'+str_c
                for ir,r in enumerate(d_raw_data[rep][plate].loc[:,c]):
                    if plate=='3' and column_plate[ic]>8:
                        continue
                    else:
                        make_key = 'plate'+plate+'-'+row_plate[ir]+str_c
                        d_data[rep][make_key] = r - d_blank[rep]
                        if d_data[rep][make_key]<0:
                            d_data[rep][make_key]=0
                        
    return(d_data)

def variance_between_replicates(d_data):
    variance = []
    for well in d_data['replicate1']:
        values = []
        for rep in d_data.keys():
            values.append(d_data[rep][well])
        
        variance.append(np.std(values))
    return(variance)
    
def high_variance_wells(variance_higher_than, variance, d_data):
    well_high_variance = []
    for iw,well in enumerate(d_data['replicate1']):
        if variance[iw]>variance_higher_than:
            print(well, [d_data[rep][well] for rep in d_data.keys()])
            well_high_variance.append(well)
            
    return(well_high_variance)
            
def compare_replicates(d_data,path_save_results):
    replicates_names = list(d_data.keys())
    for ir1 in range(len(replicates_names)-1):
        for ir2 in range(ir1+1,len(replicates_names)):
            fig, ax = plt.subplots()
            
            r1 = replicates_names[ir1]
            r2 = replicates_names[ir2]
            
            replicate_r1 = [d_data[r1][k] for k in d_data[r1]]
            replicate_r2 = [d_data[r2][k] for k in d_data[r2]]
            plt.plot(replicate_r1,replicate_r2,'ko')
                
            min_min = min([min(replicate_r1),min(replicate_r2)])
            max_max = max([max(replicate_r1),max(replicate_r2)])
            plt.plot([min_min,max_max],[min_min,max_max],'k:')
            
            plt.xlabel('OD600 '+r1)
            plt.ylabel('OD600 '+r2)
            
            plt.xlim(min_min,max_max)
            plt.ylim(min_min,max_max)
            
            # Remove axes splines
            for s in ['top', 'right']:#,'top', 'bottom', 'left', 'right']:
                ax.spines[s].set_visible(False)


            plt.show()
            fig.savefig(path_save_results+'/compare_'+r1+'_'+r2+'.pdf', bbox_inches='tight')

            
            
                        
def discard_well_from_data(d_data,discard_well_list):
    for d in discard_well_list:
        for rep in d_data.keys():
            del d_data[rep][d]
        
    return(d_data)

def write_clean_data_to_csv(save_folder, file_name, species, d_data, d_plateW_RCB):
    fil = open(save_folder+file_name,'w')
    fil.write('Well,Resources,')
    species_local = species.copy()
    species_local.reverse()

    # Escribir nombres de los recursos
    for s in species_local:
        s_name_only = s.split(' ')[0]
        fil.write(s_name_only + ',')
    
    
    
    
    
    
    #species.reverse()
    #for s in species:
        #s_name_only = s.split(' ')[0]
        #fil.write(s_name_only+',')
    
    num_reps = len(d_data.keys())
    for i in range(num_reps):
        fil.write(f'OD600 24h rep{i+1},')
    fil.write('OD600 24h mean,OD600 24h std\n')

    for well in d_data['replicate1']:
        fil.write(well+',')
        well_in_binary = d_plateW_RCB[well]
        fil.write(well_in_binary+',')
        #write species combination
        for b in well_in_binary:
            fil.write(b+',')
            
        #average and standard deviation
        rep_list = [d_data[f'replicate{i+1}'][well] for i in range(num_reps)]
        mean_replicates = np.mean(rep_list)
        std_replicates = np.std(rep_list)
        
        for val in rep_list:
            fil.write(str(val) + ',')
        fil.write(str(mean_replicates) + ',' + str(std_replicates) + '\n')

        
    fil.close()
    
def species_combinations_in_binary(n_species):
    n_combinations = 2**n_species
    
    #write combinations in binary numbers
    combinations = []
    for n in range(n_combinations):
        n_in_bin = bin(n)[2:]
        len_bin = len(n_in_bin)
        total_bin = '0'*(n_species-len_bin)+n_in_bin
        combinations.append(total_bin)
        
    return(combinations)
        

def plate_well_species_combination(n_species):
    #link plate well to species combination
    
    row_plate = string.ascii_uppercase[:8]
    column_plate = range(1,13)
    
    plate_numb = ['1','2','3']
    
    combinations = species_combinations_in_binary(n_species)
    
    d_plate = {}
    
    #d_plate has the combination of species as keys and the plate well as value
    count = 0
    for p in plate_numb:
        
        for c in column_plate:
            str_c = str(c)
            if len(str_c)==1:
                str_c = '0'+str_c
            
            for r in row_plate:
                if count<len(combinations):
                    d_plate[combinations[count]] = 'plate'+p+'-'+r+str_c
                    count = count + 1
                    
    #d_plate_reverse has the plate well as key and the combination of species as value
    d_plate_reverse = {}
    for k in d_plate:
        d_plate_reverse[d_plate[k]]=k

    return(d_plate,d_plate_reverse)

def get_fitness_in_environment(df, environment_binary):
    fitness_in_environment = df.loc[df['Resources']==environment_binary]['OD600 24h mean'].iloc[0]
    return(fitness_in_environment)

def get_fitness_std_in_environment(df, environment_binary):
    fitness_std_in_environment = df.loc[df['Resources']==environment_binary]['OD600 24h std'].iloc[0]
    return(fitness_std_in_environment)
    
def deltaF(df,env1,env2):
    f_env1 = get_fitness_in_environment(df, env1)
    f_env2 = get_fitness_in_environment(df, env2)
    Fenv1_minus_Fenv2 = f_env1 - f_env2
    return(Fenv1_minus_Fenv2)

def deltaF_std(df,env1,env2):
    f_env1_std = get_fitness_std_in_environment(df, env1)
    f_env2_std = get_fitness_std_in_environment(df, env2)
    Fenv1_minus_Fenv2_std = np.sqrt(f_env1_std**2+f_env2_std**2)
    return(Fenv1_minus_Fenv2_std)

def interaction(df,background_env,env_with_r1,env_with_r2,env_with_r1_r2):
    #calculate expected fitness
    g_background = get_fitness_in_environment(df, background_env)
    g_background_std = get_fitness_std_in_environment(df, background_env)
    
    delta_r1 = deltaF(df,env_with_r1,background_env)
    delta_r2 = deltaF(df,env_with_r2,background_env)
    
    g_expected = g_background + delta_r1 + delta_r2
    g_expected_std = np.sqrt(g_background_std**2 + deltaF_std(df,env_with_r1,background_env)**2 + deltaF_std(df,env_with_r2,background_env)**2)
    
    #observed
    g_observed = get_fitness_in_environment(df, env_with_r1_r2)
    g_observed_std = get_fitness_std_in_environment(df, env_with_r1_r2)
    
    #interaction
    epsilon = g_observed - g_expected
    epsilon_std = np.sqrt(g_observed_std**2+g_expected_std**2)
    
    return(g_observed,g_observed_std,g_expected,g_expected_std,epsilon,epsilon_std)

def high_order_interaction(df, background_env, env_with_resources, envs_ri, all_pairwise_interactions, all_pairwise_interactions_std):
    #calculate expected fitness
    g_background = get_fitness_in_environment(df, background_env)
    g_background_std = get_fitness_std_in_environment(df, background_env)
    
    all_deltas = []
    all_deltas_std_sqrt = []
    for env_ri in envs_ri:
        #deltas
        delta_ri = deltaF(df,env_ri,background_env)
        all_deltas.append(delta_ri)
        #deltas std
        delta_ri_std = deltaF_std(df,env_ri,background_env)
        all_deltas_std_sqrt.append(delta_ri_std**2)
    
    g_expected = g_background + sum(all_deltas) + sum(all_pairwise_interactions)
    g_expected_std = np.sqrt(g_background_std**2 + sum(all_deltas_std_sqrt) + sum([pi_std**2 for pi_std in all_pairwise_interactions_std]))
    
    #observed
    g_observed = get_fitness_in_environment(df, env_with_resources)
    g_observed_std = get_fitness_std_in_environment(df, env_with_resources)
    
    #interaction
    epsilon = g_observed - g_expected
    epsilon_std = np.sqrt(g_observed_std**2+g_expected_std**2)
    
    return(g_observed,g_observed_std,g_expected,g_expected_std,epsilon,epsilon_std)

def possible_environments(n_resources):
    all_envs = [comb for comb in product(['0','1'], repeat=n_resources)]
    return(all_envs)

def environment_to_string(env_tuple):
    env_string = ''.join(env_tuple)
    return(env_string)
    
def environment_in_data(df, env):
    missing_data = df.loc[df['Resources']==env].empty
    if missing_data == False:
        env_in_data = 'True'
    else:
        env_in_data = 'False'
        
    return(env_in_data)

def get_FEEs(f_background_all_mean, f_effect_ri_all_mean):
    #fit
    slope, intercept, r_value, p_value, std_err = stats.linregress(f_background_all_mean, f_effect_ri_all_mean)
    return(slope, intercept, r_value, p_value, std_err)

def linear_func(p, x):
    m, c = p
    return m*x + c

def fit_using_total_least_squares(x,y):
    linear_model = Model(linear_func)
    
    data = RealData(x, y)
    
    odr = ODR(data, linear_model, beta0=[1, 0])
    out = odr.run()
    
    #out.pprint()
    slope = out.beta[0]
    intercept = out.beta[1]
    
    beta_0 = 0  # test if slope is significantly different from zero
    t_stat = (slope - beta_0) / out.sd_beta[0]  # t statistic for the slope parameter
    df = out.iwork[10] # degrees of freedom (n_sample-n_parameters)
    p_value = stats.t.sf(np.abs(t_stat), df) * 2

    y_predicted = [x_val*slope+intercept for x_val in x]
    r2_score_value = r2_score(y, y_predicted)
    
    return(slope, intercept, p_value, r2_score_value)
    
# def backgrounds(n_species):
#     #get species combinations in binary numbers
#     combinations = species_combinations_in_binary(n_species)
    
#     d_backgrounds = {}
#     for s in range(n_species):
#         #get all combinations with the element s in
#         these_ones = [c for c in combinations if c[s]=='1']#128
#         d_backgrounds[s] = these_ones
        
#     return(d_backgrounds)
    
# def get_deltaF_per_replicate(species,d_data,d_plate):
#     #calculate deltaF for each species and replicate
    
#     n_species = len(species)
#     d_backgrounds = backgrounds(n_species)
    
#     #dictionary to store the results
#     d_deltaF = {}
#     for r in d_data.keys():
#         d_deltaF[r] = {}
    
#         for s in range(n_species):
            
#             d_deltaF[r][species[s]] = {}
            
#             #get all combinations with the element s in
#             these_ones = d_backgrounds[s]
            
#             F_background_all = []
#             F_background_with_ri = []
#             deltaF = []
#             deltaF_reverse = []
#             backgrounds_for_s = []

#             for t in these_ones:
#                 #
#                 background = list(t)
#                 background[s] = '0'
#                 background = ''.join(background)
#                 #
#                 well_background = d_plate[background]
#                 well_with_s = d_plate[t]
#                 #
#                 if (well_background in d_data[r]) and (well_with_s in d_data[r]):
#                     F_background = d_data[r][well_background]
#                     F_background_all.append(F_background)
#                     F_s = d_data[r][well_with_s]
#                     F_background_with_ri.append(F_s)
#                     deltaF.append(F_s-F_background)
#                     deltaF_reverse.append(F_background-F_s)
#                     backgrounds_for_s.append(background)
  
#             d_deltaF[r][species[s]]['function_background'] = F_background_all
#             d_deltaF[r][species[s]]['function_background_with_ri'] = F_background_with_ri
#             d_deltaF[r][species[s]]['delta_function'] = deltaF
#             d_deltaF[r][species[s]]['delta_function_reverse'] = deltaF_reverse
#             d_deltaF[r][species[s]]['backgrounds'] = backgrounds_for_s
         
#     return(d_deltaF)


# def average_deltaF(d_deltaF,species):
#     #calculate average among all replicates
#     d_deltaF_mean = {}
#     for s in d_deltaF['replicate1']:
#         d_deltaF_mean[s] = {}
#         d_deltaF_mean[s]['function_background'] = {}
#         d_deltaF_mean[s]['function_background']['mean'] = []
#         d_deltaF_mean[s]['function_background']['std'] = []
#         d_deltaF_mean[s]['function_background_with_ri'] = {}
#         d_deltaF_mean[s]['function_background_with_ri']['mean'] = []
#         d_deltaF_mean[s]['function_background_with_ri']['std'] = []
#         d_deltaF_mean[s]['delta_function'] = {}
#         d_deltaF_mean[s]['delta_function']['mean'] = []
#         d_deltaF_mean[s]['delta_function']['std'] = []
#         d_deltaF_mean[s]['delta_function_reverse'] = {}
#         d_deltaF_mean[s]['delta_function_reverse']['mean'] = []
#         d_deltaF_mean[s]['delta_function_reverse']['std'] = []
        
#         for i in range(len(d_deltaF['replicate1'][s]['delta_function'])):
#             fun_back_i = [d_deltaF[rep][s]['function_background'][i] for rep in d_deltaF]
#             d_deltaF_mean[s]['function_background']['mean'].append(np.mean(fun_back_i))
#             d_deltaF_mean[s]['function_background']['std'].append(np.std(fun_back_i))
            
#             fun_back_i_with_r = [d_deltaF[rep][s]['function_background_with_ri'][i] for rep in d_deltaF]
#             d_deltaF_mean[s]['function_background_with_ri']['mean'].append(np.mean(fun_back_i_with_r))
#             d_deltaF_mean[s]['function_background_with_ri']['std'].append(np.std(fun_back_i_with_r))

#             deltafun_i = [d_deltaF[rep][s]['delta_function'][i] for rep in d_deltaF]
#             d_deltaF_mean[s]['delta_function']['mean'].append(np.mean(deltafun_i))
#             d_deltaF_mean[s]['delta_function']['std'].append(np.std(deltafun_i))

#             deltafunrev_i = [d_deltaF[rep][s]['delta_function_reverse'][i] for rep in d_deltaF]
#             d_deltaF_mean[s]['delta_function_reverse']['mean'].append(np.mean(deltafunrev_i))
#             d_deltaF_mean[s]['delta_function_reverse']['std'].append(np.std(deltafunrev_i))

#     return(d_deltaF_mean)
        
# def fit_deltaF_backgroundF(d_deltaF,species_name):
#     F_background_all_replicates = []
#     deltaF_replicates = []
    
#     for r in d_deltaF:
#         F_background_all_replicates = F_background_all_replicates + d_deltaF[r][species_name]['function_background']
#         deltaF_replicates = deltaF_replicates + d_deltaF[r][species_name]['delta_function']
#     #fit
#     slope, intercept, r_value, p_value, std_err = stats.linregress(F_background_all_replicates, deltaF_replicates)
#     return(slope, intercept, r_value, p_value, std_err)      

# def linear_func(p, x):
#    m, c = p
#    return m*x + c

# def fit_using_total_least_squares(x,y):
#     linear_model = Model(linear_func)
    
#     data = RealData(x, y)
    
#     odr = ODR(data, linear_model, beta0=[1, 0])
#     out = odr.run()
    
#     #out.pprint()
#     slope = out.beta[0]
#     intercept = out.beta[1]
    
#     beta_0 = 0  # test if slope is significantly different from zero
#     t_stat = (slope - beta_0) / out.sd_beta[0]  # t statistic for the slope parameter
#     df = out.iwork[10] # degrees of freedom (n_sample-n_parameters)
#     p_value = stats.t.sf(np.abs(t_stat), df) * 2

#     return(slope, intercept, p_value)

# def fit_F_background_with_ri_vs_backgroundF(d_deltaF,species_name,method):
#     F_background_all_replicates = []
#     F_background_with_ri_all_replicates = []
    
#     for r in d_deltaF:
#         F_background_all_replicates = F_background_all_replicates + d_deltaF[r][species_name]['function_background']
#         F_background_with_ri_all_replicates = F_background_with_ri_all_replicates + d_deltaF[r][species_name]['function_background_with_ri']
    
#     if method=='ordinary': #ordinary least squares fit
#         slope, intercept, r_value, p_value, std_err = stats.linregress(F_background_all_replicates, F_background_with_ri_all_replicates)
#     elif method=='total': #total least squares fit
#         slope, intercept, p_value = fit_using_total_least_squares(F_background_all_replicates, F_background_with_ri_all_replicates)
    
#     return(slope, intercept, p_value)      

    
# def min_max_landscape_values(d_deltaF):
#     F_background_all = []
#     F_background_with_ri_all = []
#     deltaF_all = []
    
#     for r in d_deltaF:
#         for s in d_deltaF[r]:
#             F_background_all = F_background_all + d_deltaF[r][s]['function_background']
#             F_background_with_ri_all = F_background_with_ri_all + d_deltaF[r][s]['function_background_with_ri']
#             deltaF_all = deltaF_all + d_deltaF[r][s]['delta_function']

#     min_background = min(F_background_all)
#     max_background = max(F_background_all)
#     min_background_with_ri = min(F_background_with_ri_all)
#     max_background_with_ri = max(F_background_with_ri_all)
#     min_deltaF = min(deltaF_all)
#     max_deltaF = max(deltaF_all)
    
#     return(min_background,max_background,min_background_with_ri,max_background_with_ri,min_deltaF,max_deltaF)


    