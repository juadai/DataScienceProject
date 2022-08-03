import pandas as pd
import nltk

# Load dataframes
observations = pd.read_csv('understoreydata2020.csv')
flora = pd.read_excel('./Flora and Fauna species 2007.xls')
matching = pd.read_csv('understorey_MATCHING.csv')

# Drop extra columns
to_drop = ['Observer1', 'Observer2', 'Plot_percent_burgan', 'Quadrat_Can_Cov', 'Percentage_range']
observations.drop(to_drop, axis=1, inplace=True)

to_drop = ['Origin', 'EPBC', 'VROTS', 'FFG', 'Species Number']
flora.drop(to_drop, axis=1, inplace=True)


# Remove the observations related to inter-tussock gap
inter_tussock = ['CWD', 'Leaf Litter', 'Moss/Lichen', 'Bare Earth']
observations = observations.loc[~observations['T002_Flora_Species_name'].isin(inter_tussock)]


# Set up key on observations for join
matching.index = matching['Entry']
matching = matching['Candidate']
matching.dropna(inplace=True)

observations['Key'] = None
for i in observations.index:
    if observations.loc[i, 'T002_Flora_Species_name'] in matching.index:
        observations.loc[i, 'Key'] = matching[observations.loc[i, 'T002_Flora_Species_name']]
    else:
        observations.loc[i, 'Key'] = observations.loc[i, 'T002_Flora_Species_name']
    

# Left join to preserve all observations
df = observations.merge(flora, how='left',
                        left_on='Key',
                        right_on='Scientific Name')



