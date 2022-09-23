import pandas as pd

df = pd.read_csv('../dataset/joined_understorey_observations.csv')

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
         'Life Form',
         'T002_Flora_Species_name',
         'Abundance']]

df.rename(columns={'Plot_treatment': 'Treatment',
                   'Plot_number': 'Plot Number',
                   'Quadrat_number': 'Quadrat Number',
                   'Quadrat_fenced': 'Fenced',
                   'Quadrat_gap': 'Gap',
                   'Year_monitoring': 'Year',
                   'T002_Flora_Species_name':'Species Name'}, inplace=True)

lf_abund = df.drop('Species Name', axis=1).groupby(
                            ['Treatment',
                             'Plot Number',
                             'Quadrat Number',
                             'Fenced',
                             'Gap',
                             'Life Form',
                             'Year']).sum()


###########################
# To aggregate abundances #
###########################

# Aggregate life form abundance over quadrats?
quads = True

# Aggregate life form abundance over plots as well?
plots = True


if quads and not plots:
    lf_abund = lf_abund.groupby(['Treatment',
                                 'Plot Number',
                                 'Fenced',
                                 'Gap',
                                 'Life Form',
                                 'Year']).mean()


if quads and plots:
    lf_abund = lf_abund.groupby(['Treatment',
                                 'Fenced',
                                 'Gap',
                                 'Life Form',
                                 'Year']).mean()

####################
# Convert to table #
####################

lf_abund = lf_abund['Abundance']

rows = lf_abund.index.droplevel(lf_abund.index.names[-2:])
rows = rows.drop_duplicates()

cols = lf_abund.index.droplevel(lf_abund.index.names[:-2])
cols = cols.drop_duplicates()

matrix_lf_abund = pd.DataFrame(index=rows, columns=cols)

for i in lf_abund.index:
    matrix_lf_abund.at[i[:-2],i[-2:]] = lf_abund[i]

matrix_lf_abund.fillna(0, inplace=True)

# Add any missing columns and sort
for i in matrix_lf_abund.columns.levels[0]:
    for j in [0, 3, 6]:
        if (i, j) not in matrix_lf_abund.columns:
            matrix_lf_abund.loc[:,(i, j)] = 0
			
matrix_lf_abund.sort_index(axis=1, inplace=True)
matrix_lf_abund.sort_index(inplace=True)


# Caclulate change from year 0 to year 6
#for i in matrix_lf_abund.columns.levels[0]:
#    matrix_lf_abund.loc[:,(i, 'Change 0-6')] = matrix_lf_abund[(i, 6)] - matrix_lf_abund[(i, 0)]



# Calculate relative life form abundances
matrix_lf_abund_rel = matrix_lf_abund.divide(matrix_lf_abund.sum(axis=1), axis=0)


matrix_lf_abund.to_csv('./output/midpoint/abundance/lf_abund_treatment.csv')


# Format for analysis
analysis_lf_abund = matrix_lf_abund.stack(-2)
analysis_lf_abund_rel = matrix_lf_abund_rel.stack(-2)

analysis_lf_abund = analysis_lf_abund.loc[~((analysis_lf_abund[0] == 0) &
                                            (analysis_lf_abund[3] == 0) &
                                            (analysis_lf_abund[6] == 0)),:]

analysis_lf_abund_rel = analysis_lf_abund_rel.loc[~((analysis_lf_abund_rel[0] == 0) &
                                                    (analysis_lf_abund_rel[3] == 0) &
                                                    (analysis_lf_abund_rel[6] == 0)),:]

analysis_lf_abund = analysis_lf_abund.stack()
analysis_lf_abund_rel = analysis_lf_abund_rel.stack()

# analysis_lf_abund.to_csv('./output/midpoint/abundance/for analysis/lf_abund_analysis.csv')
# analysis_lf_abund_rel.to_csv('./output/midpoint/abundance/for analysis/lf_abund_rel_analysis.csv')
