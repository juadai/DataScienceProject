import pandas as pd

df = pd.read_csv('../dataset/joined_understorey_observations.csv')
df['Shade Tolerant'] = False

shade_df = pd.read_csv('../../Reference/shade_tolerant_species.csv')
shade_df = shade_df.loc[~shade_df['shade tolerant'].isna(), 'Observed Species']
shade_df = list(shade_df)


for i in df.index:
    if df.loc[i, 'T002_Flora_Species_name'] in shade_df:
        df.loc[i, 'Shade Tolerant'] = True

df.to_csv('../dataset/joined_understorey_observations.csv', index=False)
