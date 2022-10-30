import pandas as pd
import numpy as np

df = pd.read_csv('../dataset/joined_understorey_observations.csv')

def Shannons (vector):

    p = [i/100 for i in vector]
    p = [i/sum(p) for i in p]
    p = [np.log(i)*i for i in p]
    
    return(-sum(p))

SCORES = {0.5: [0.5,0.5],
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

for i in df.index:
    df.loc[i, 'Abundance'] = sum(SCORES[df.loc[i, 'Score']])/2


df = df[['Plot_treatment',
         'Plot_number',
         'Quadrat_number',
         'Quadrat_fenced',
         'Quadrat_gap',
         'Year_monitoring',
         'T002_Flora_Species_name',
         'Abundance']]

df.rename(columns={'Plot_treatment': 'Treatment',
                   'Plot_number': 'Plot Number',
                   'Quadrat_number': 'Quadrat Number',
                   'Quadrat_fenced': 'Fenced',
                   'Quadrat_gap': 'Gap',
                   'Year_monitoring': 'Year',
                   'T002_Flora_Species_name':'Species Name'}, inplace=True)

species_divers = df.groupby(['Treatment',
                             'Plot Number',
                             'Quadrat Number',
                             'Fenced',
                             'Gap',
                             'Year',
                             'Species Name']).sum()


###########################
# To aggregate abundances #
###########################

# Aggregate life form abundance over quadrats?
quads = True

# Aggregate life form abundance over plots as well?
plots = True



if quads and not plots:
    species_divers = species_divers.groupby(['Treatment',
                                             'Plot Number',
                                             'Fenced',
                                             'Gap',
                                             'Year',
                                             'Species Name']).sum()


if quads and plots:
    species_divers = species_divers.groupby(['Treatment',
                                             'Fenced',
                                             'Gap',
                                             'Year',
                                             'Species Name']).sum()






rows = species_divers.index.droplevel(species_divers.index.names[-2:])
rows = rows.drop_duplicates()

cols = [0, 3, 6]

matrix_species_divers = pd.DataFrame(index=rows, columns=cols)

to_calculate = species_divers.index.droplevel(['Species Name']).drop_duplicates()

for i in to_calculate:
    
    vector = species_divers.loc[i, 'Abundance']
    SDI = Shannons(vector)
    
    matrix_species_divers.loc[i[:-1], i[-1:]] = SDI

    
#matrix_species_divers['Change 0-6'] = matrix_species_divers[6] - matrix_species_divers[0]

matrix_species_divers.to_csv('./output/midpoint/diversity/species_divers_treatment.csv')
