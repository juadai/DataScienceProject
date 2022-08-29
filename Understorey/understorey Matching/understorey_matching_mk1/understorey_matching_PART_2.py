import pandas as pd
import numpy as np
import nltk

######################################
# Import data and discard extra rows #
######################################

flora1 = pd.read_excel('./Reference/Flora and Fauna species 2008.xls')
flora2 = pd.read_csv('./Reference/VBA Species-Checklist 100921.csv')

observ = pd.read_csv('../dataset/understoreydata2020.csv')


SCORES = {0.5: 'Trace <1%',
          1: '1-5%',
          2: '5-10%',
          3: '10-20%',
          4: '20-30%',
          5: '30-40%',
          6: '40-50%',
          7: '50-60%',
          8: '60-70%',
          9: '70-80%',
          10: '80-90%'}

flora1 = flora1[['Scientific Name',
                 'Common Name',
                 'Family Name',
                 'Division Name',
                 'Life Form']]

# Keep Flora only
flora2 = flora2.loc[flora2['PRIMARY_DISCIPLINE']=='Flora']
# Drop rows from flora2 where the scientific name is in flora1
flora2 = flora2.loc[~flora2['SCIENTIFIC_NAME'].isin(flora1['Scientific Name'])]


flora2 = flora2[['SCIENTIFIC_NAME',
                 'COMMON_NAME',
                 'NVIS_GROWTHFORM',
                 'ORIGIN']]

observ = observ[['Year_monitoring',
                 'Plot_number',
                 'Plot_treatment',
                 'Quadrat_number',
                 'Quadrat_fenced',
                 'Quadrat_gap_or_no',
                 'Record_ID',
                 'T002_Flora_Species_name',
                 'Score']]



#############################################################
# Join and consolidate flora data under a consistent schema #
#############################################################

flora2.columns = ['Scientific Name', 'Common Name', 'Life Form', 'Origin']
cols = ['Scientific Name', 'Common Name', 'Family Name',
        'Division Name', 'Life Form', 'Origin']
for c in cols:
    for f in [flora1, flora2]:
        if c not in f.columns:
            f[c] = np.NaN

flora = pd.concat([flora1, flora2])
flora = flora.loc[~flora['Life Form'].isna()]

#################################
# Load Manually Matched Species #
#################################

matched = pd.read_csv('./Manual Matching/to_match.csv')
matched.index = matched['Observation']
matched = matched['Candidate 1']
matched.dropna(inplace=True)

matched_extra = pd.read_csv('./Manual Matching/matching_round_2.csv')
matched_extra.index = matched_extra['0']
matched_extra = matched_extra['new_key']


######################################
# Join observations with flora table #
######################################

observ['key'] = None

for i in observ.index:
    # Check if observation is a perfect match
    if observ.loc[i, 'T002_Flora_Species_name'] in list(flora['Scientific Name']):
        observ.loc[i, 'key'] = observ.loc[i, 'T002_Flora_Species_name']
    # If not, check if we've manuallly matched it
    elif observ.loc[i, 'T002_Flora_Species_name'] in list(matched.index):
        observ.loc[i, 'key'] = matched[observ.loc[i, 'T002_Flora_Species_name']]

    if observ.loc[i, 'key'] in list(matched_extra.index):
        observ.at[i, 'key'] = matched_extra[observ.loc[i, 'key']]



observ = observ.merge(flora1, left_on='key', right_on='Scientific Name', how='left')

observ.drop(['key'], inplace=True, axis=1)
observ.reset_index(drop=True, inplace=True)

######################################
# Join observations with flora table #
######################################

#weeds = pd.read_excel('./RGW species list for weeds.xlsx')




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
                 'MoSs/Lichen']

observ = observ.loc[~observ['T002_Flora_Species_name'].isin(inter_tussock)]

# observ.to_csv('./joined_understorey_observations.csv', index=False)


###########################################
# Try to narrow # of life form categories #
###########################################

observ_to_fix = observ.loc[observ['Life Form'].isin(pd.unique(flora2['Life Form']))]

species_to_rematch = pd.unique(observ_to_fix['Scientific Name'])


