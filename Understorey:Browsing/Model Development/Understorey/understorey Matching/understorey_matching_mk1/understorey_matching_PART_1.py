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

############################################
# Define function for approximate matching #
############################################

def lev(word, candidates, n_options=15):
    scores = []
    
    for c in candidates:
        score = nltk.edit_distance(word.lower(), c)
        scores.append((score, c))
    return(sorted(scores)[:n_options])

#####################################
# Find items without an exact match #
#####################################

to_match = []

for s in observ['T002_Flora_Species_name'].unique(): 
    if s not in list(flora['Scientific Name']):
        to_match.append(s)

print(len(to_match), 'records to match\n')

match_df = pd.DataFrame(None, columns=['Top candidate', 'Other candidates', 'Options'])
j = 0
for i in to_match:
    match_df.loc[i, 'Options'] = str(lev(i, flora['Scientific Name']))
    j += 1
    if not j%5: print(j)


print("Lev done!")

match_df.to_csv('./to_match.csv')
