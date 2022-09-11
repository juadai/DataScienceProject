import numpy as np

def Shannons (vector):

    p = [i/100 for i in vector]
    p = [i/sum(p) for i in p]
    p = [np.log(i)*i for i in p]
    
    return(-sum(p))

import pandas as pd
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


df = pd.read_csv('./dataset/joined_understorey_observations.csv')

life_forms = pd.unique(df['Life Form'])


for i in df.index:
    df.loc[i, 'Abundance'] = sum(SCORES[df.loc[i,'Score']])/2


#Keeps the treatment, gap, fencing, plot number quadrat number
'''
df.drop(['Record_ID', 'Score'], inplace=True, axis=1)

df = df.groupby(by=['Plot_treatment',
                    'Plot_number', 
                    'Quadrat_fenced',
                    'Quadrat_gap',
                    'Quadrat_number',
                    'Year_monitoring',
                    'Scientific Name']).sum()
'''


# Keeps the treatment, gap, fencing

df.drop(['Record_ID', 'Score', 'Plot_number', 'Quadrat_number'], inplace=True, axis=1)

df = df.groupby(by=['Plot_treatment',
                    'Quadrat_fenced',
                    'Quadrat_gap',
                    'Year_monitoring',
                    'Scientific Name']).sum()



# Just keeps the treatment and the year
'''
df.drop(['Record_ID', 'Score', 'Plot_number', 'Quadrat_number',
                    'Quadrat_fenced',
                    'Quadrat_gap',], inplace=True, axis=1)


df = df.groupby(by=['Plot_treatment',
                    'Year_monitoring',
                    'Scientific Name']).sum()
'''



new_index = [i[0:-1] for i in df.index]
new_index = list(set(new_index))

df['Shannon']=0

for i in new_index:
    df.loc[i,'Shannon'] = Shannons(df.loc[i]['Abundance'].values)

df.index = df.index.droplevel('Scientific Name')
df = df['Shannon']
df = df.drop_duplicates()

df = pd.DataFrame(df)
df.reset_index('Year_monitoring', inplace=True)


matrix = pd.DataFrame(None, index=df.index, columns=[0,3,6])

for i in df.index:
    items = df.loc[i,:]
    for j in range(0,len(items)):
        year = items.iloc[j,0]
        val = items.iloc[j,1]
        matrix.loc[i,year] = val

matrix = matrix.drop_duplicates()




