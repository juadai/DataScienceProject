import pandas as pd

s1 = pd.read_csv('./species.csv')

df = pd.read_excel('./Flora and Fauna species 2007.xls')
s2 = list(df['Scientific Name'])

i = 0
for s in list(s1['x'].values):
  if s not in s2:
    print(s)
    i += 1

print(i, "out of", len(s1), "species to be matched")

