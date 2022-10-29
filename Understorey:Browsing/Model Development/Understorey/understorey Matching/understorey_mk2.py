import pandas as pd
import numpy as np
import nltk


######################################################
# Join with flora info incorporating manual matching #
######################################################


flora = pd.read_excel('./Reference/Flora and Fauna species 2008.xls')
flora = flora[['Origin',
               'Scientific Name',
               'Common Name',
               'Family Name',
               'Life Form']]

observ = pd.read_csv('../dataset/understoreydata2020.csv')
observ = observ[['Year_monitoring',
                 'Plot_number',
                 'Plot_treatment',
                 'Quadrat_number',
                 'Quadrat_fenced',
                 'Quadrat_gap_or_no',
                 'Record_ID',
                 'T002_Flora_Species_name',
                 'Score']]
    
manual_matching = pd.read_csv('./Manual Matching/matching_round_2.csv')
manual_matching.index = manual_matching['0']
manual_matching = manual_matching['new_key']



observ['key'] = observ['T002_Flora_Species_name']
for i in observ.index:
    if observ.loc[i, 'key'] in manual_matching.index:
        observ.loc[i, 'key'] = manual_matching[observ.loc[i, 'key']]



observ = observ.merge(flora, left_on='key', right_on='Scientific Name')

observ.drop(['key'], axis=1, inplace=True)


#######################################
# Complete in gap / not in gap fields #
#######################################

observ.reset_index(drop=True, inplace=True)
observ.rename(columns={'Quadrat_gap_or_no': 'Quadrat_gap'}, inplace=True)

# CHECK AND DISCUSS THIS LOGIC
observ.at[observ['Quadrat_gap']=='control', 'Quadrat_gap'] = 'Not in Gap'
observ.at[observ['Quadrat_gap']=='natural gap', 'Quadrat_gap'] = 'Not in Gap'
observ.at[observ['Quadrat_gap']=='Edge of Gap', 'Quadrat_gap'] = 'In Gap'

# Plots 4 and 11 were never actually fenced due to budged (but records at year 0 reflect the intention)
observ.at[(observ['Plot_number']==4) | (observ['Plot_number']==11), 'Quadrat_fenced'] = False

for i in pd.unique(observ['Plot_number']):

    quadrats = pd.unique(observ.loc[observ['Plot_number']==i]['Quadrat_number'])
    
    for q in quadrats:
        
        a = observ.loc[(observ['Plot_number']==i) & (observ['Quadrat_number']==q)]
        
        g = a['Quadrat_gap']
        g = pd.unique(g[~g.isna()])
        #if len(g) != 1:
        #    print('GAP CONFLICT: Plot', i, 'Quadrat', q)
        #    print(a[['Year_monitoring', 'Quadrat_gap']].sort_values('Year_monitoring'))
        
        # FOR THE TIME BEING - WHERE THERE'S A CONFLICT, BREAK TIE ARBITRARILY
        g = g[0]
        
        f = a['Quadrat_fenced']
        f = pd.unique(f[~f.isna()])
        if len(f) != 1: print('FENCED CONFLiCT: Plot', i, 'Quadrat', q)
        f = f[0]
        
        observ.at[a.index, 'Quadrat_fenced'] = f
        observ.at[a.index, 'Quadrat_gap'] = g

observ.at[observ['Quadrat_gap']=='In Gap', 'Quadrat_gap'] = True
observ.at[observ['Quadrat_gap']=='Not in Gap', 'Quadrat_gap'] = False


observ['Life Form'].fillna('Unknown', inplace=True)

inter_tussock = ['CWD',
                 'Leaf Litter',
                 'Moss/Lichen',
                 'Bare Earth',
                 'Dead Standing Timber',
                 'MoSs/Lichen',
                 'Unknown']

observ = observ.loc[~observ['T002_Flora_Species_name'].isin(inter_tussock)]

observ.to_csv('./outputs/joined_understorey_observations.csv', index=False)
