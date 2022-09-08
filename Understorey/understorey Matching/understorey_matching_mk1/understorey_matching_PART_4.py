import pandas as pd
import numpy as np
import nltk


df = pd.read_csv('/Users/nickjolly/Desktop/tfn/DataScienceProject/Understorey/understorey Matching/outputs/abundances_sum.csv')


# Average over quadrats
df = df.groupby(['Plot_treatment',
                 'Plot_number',
                 'Quadrat_fenced',
                 'Quadrat_gap',
                 'Year_monitoring',]).mean()

df.drop(['Quadrat_number'], axis=1, inplace=True)


# And now average over plots
df = df.groupby(['Plot_treatment',
                 'Quadrat_fenced',
                 'Quadrat_gap',
                 'Year_monitoring']).mean()
