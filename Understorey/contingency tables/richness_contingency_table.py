import pandas as pd
import numpy as np

df = pd.read_csv('../dataset/joined_understorey_observations.csv')


SCORES = {0.5: [0.5],
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
    df.loc[i, 'Abundance'] = min(SCORES[df.loc[i, 'Score']])


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

species_rich = df.groupby(['Treatment',
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
quads = False

# Aggregate life form abundance over plots as well?
plots = False



if quads and not plots:
    species_rich = species_rich.groupby(['Treatment',
                                         'Plot Number',
                                         'Fenced',
                                         'Gap',
                                         'Year',
                                         'Species Name']).sum()


if quads and plots:
    species_rich = species_rich.groupby(['Treatment',
                                         'Fenced',
                                         'Gap',
                                         'Year',
                                         'Species Name']).sum()




rows = species_rich.index.droplevel(species_rich.index.names[-2:])
rows = rows.drop_duplicates()

cols = [0, 3, 6]

matrix_species_richness = pd.DataFrame(index=rows, columns=cols)

to_calculate = species_rich.index.droplevel(['Species Name']).drop_duplicates()

for i in to_calculate:
    
    vector = species_rich.loc[i, 'Abundance'].index
    
    matrix_species_richness.loc[i[:-1], i[-1:]] = len(pd.unique(vector))

    
matrix_species_richness['Change 0-6'] = matrix_species_richness[6] - matrix_species_richness[0]

matrix_species_richness.to_csv('./output/spcies_richness_quadrats.csv')
