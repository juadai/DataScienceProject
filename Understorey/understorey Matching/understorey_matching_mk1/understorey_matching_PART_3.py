import pandas as pd
import numpy as np
import nltk

SCORES = {0.5: [0,1],
          1: [1,5],
          2: [5,10],
          3: [10,20],
          4: [20,30],
          5: [30,40],
          6: [40,50],
          7: [50,60],
          8: [60,70],
          9: [70,80],
          10: [80,90]}

######################################
# Import data and discard extra rows #
######################################

df = pd.read_csv('../../dataset/joined_understorey_observations.csv')

life_forms = pd.unique(df['Life Form'])


for i in df.index:
    df.loc[i, 'Abundance'] = sum(SCORES[df.loc[i,'Score']])/2

df.drop(['Record_ID', 'Score'], inplace=True, axis=1)


df = df.groupby(by=['Plot_treatment',
                    'Plot_number',
                    'Quadrat_fenced',
                    'Quadrat_gap',
                    'Quadrat_number',
                    'Year_monitoring',
                    'Life Form']).sum()

df.reset_index('Life Form', inplace=True)


matrix = pd.DataFrame(None, index=df.index, columns=pd.unique(df['Life Form']))

for i in df.index:
    items = df.loc[i,:]
    for j in range(0, len(items)):
        lf = items.iloc[j, 0]
        val = items.iloc[j, 1]
        matrix.loc[i, lf] = val

matrix.drop_duplicates(inplace=True)

matrix.fillna(0, inplace=True)

matrix.to_csv('../outputs/abundances_sum.csv')

# Relative life form abundances
# Relative abundance of a species =
# its abundance as a proportion of total abundance in the quadrat
# (Original score = abundance relative to the size of the quadrat)
totals = matrix.sum(axis=1)
matrix = matrix.divide(totals, axis='rows')


matrix.to_csv('../outputs/observations_rel_abundances.csv')
