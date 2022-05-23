import pandas as pd
import os

TREATMENTS = {1: 'T1: Gap',
              2: 'Control',
              3: 'T2: Radial',
              4: 'Control',
              5: 'T2: Radial',
              6: 'Control',
              7: 'T2: Radial',
              8: 'T1: Gap',
              9: 'T2: Radial',
              10: 'T1: Gap',
              11: 'Control',
              12: 'T1: Gap'}


def process(path):
    
    plot = int(path[-7:-5])
    print('Plot', plot, '...')

    df_input = pd.read_excel(path)

    df_input_0 = df_input.loc[df_input['Year']==0]
    df_input_6 = df_input.loc[df_input['Year']==6]



    # Build output dataframe and set tree numbers as index
    df_tree = pd.DataFrame(columns=['Tree Number',
                                    'Plot',
                                    'Treatment',
                                    'Large',
                                    'Alive at year 6',
                                    'DBH year 0',
                                    'DBH year 6',
                                    'Class year 0',
                                    'Class year 6',
                                    'Hollows Class 1 year 0',
                                    'Hollows Class 2 year 0',
                                    'Hollows Class 3 year 0',
                                    'Hollows Class 1 year 6',
                                    'Hollows Class 2 year 6',
                                    'Hollows Class 3 year 6'])

    df_tree['Tree Number'] = df_input['Tree Number'].unique()
    df_tree.set_index('Tree Number', inplace=True)

    df_tree['Plot'] = plot
    df_tree['Treatment'] = TREATMENTS[plot]

    # Iterate over year 0 and year 6 records, filling in the tree dataframe
    for i in df_input_0['Tree Number']:

        try:
            df_tree.at[i, 'DBH year 0'] = df_input_0.loc[df_input_0['Tree Number']==i, 'DBH'].item()
            df_tree.at[i, 'Large'] = df_input_0.loc[df_input_0['Tree Number']==i, 'Large?'].item()

            df_tree.at[i, 'Hollows Class 1 year 0'] = df_input_0.loc[df_input_0['Tree Number']==i, 'Hollows Class 1'].item()
            df_tree.at[i, 'Hollows Class 2 year 0'] = df_input_0.loc[df_input_0['Tree Number']==i, 'Hollows Class 2'].item()
            df_tree.at[i, 'Hollows Class 3 year 0'] = df_input_0.loc[df_input_0['Tree Number']==i, 'Hollows Class 3'].item()
        except:
            print(i)

    for i in df_input_6['Tree Number']:

        try:
            df_tree.at[i, 'DBH year 6'] = df_input_6.loc[df_input_6['Tree Number']==i, 'DBH'].item()
            df_tree.at[i, 'Alive at year 6'] = df_input_6.loc[df_input_6['Tree Number']==i, 'Alive'].item()
            df_tree.at[i, 'Large'] = df_input_6.loc[df_input_6['Tree Number']==i, 'Large?'].item()
            
            df_tree.at[i, 'Hollows Class 1 year 6'] = df_input_6.loc[df_input_6['Tree Number']==i, 'Hollows Class 1'].item()
            df_tree.at[i, 'Hollows Class 2 year 6'] = df_input_6.loc[df_input_6['Tree Number']==i, 'Hollows Class 2'].item()
            df_tree.at[i, 'Hollows Class 3 year 6'] = df_input_6.loc[df_input_6['Tree Number']==i, 'Hollows Class 3'].item()
        except:
            print(i)

    outpath = 'output/processed_'+str(plot)+'.csv'
    df_tree.to_csv(outpath)
    
    print('done\n')
    return df_tree

paths = os.listdir('./input')
paths = ['./input/' + i for i in paths]


df = pd.DataFrame(None)
for i in paths:
    df = pd.concat([df, process(i)])

df.to_csv('output/all_plots.csv')
