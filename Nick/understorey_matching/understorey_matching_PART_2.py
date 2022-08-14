import pandas as pd
import numpy as np
import nltk

######################################
# Import data and discard extra rows #
######################################

flora1 = pd.read_excel('./Flora and Fauna species 2007.xls')
flora2 = pd.read_csv('./VBA Species-Checklist 100921.csv')

observ = pd.read_csv('./understoreydata2020.csv')

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

matched = pd.read_csv('./to_match.csv', index_col=0)
matched = matched['Candidate 1']
matched.dropna(inplace=True)


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

observ = observ.merge(flora, left_on='key', right_on='Scientific Name')

observ.drop(['key'], inplace=True, axis=1)

######################################
# Join observations with flora table #
######################################

weeds = pd.read_excel('./RGW species list for weeds.xlsx')



#observ.to_csv('./joined_understorey_observations.csv')



