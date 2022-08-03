import pandas as pd
import nltk

s1 = pd.read_csv('./species.csv')

df = pd.read_excel('./Flora and Fauna species 2007.xls')
s2 = list(df['Scientific Name'])
s3 = [s.lower() for s in s2]

not_matched = []

i = 0
for s in list(s1['x'].values):
  if s.lower() not in s3:
    not_matched.append(s)
    i += 1

print(i, "out of", len(s1), "species still to be matched")


'''
j = 0
df_lev = pd.DataFrame(None, columns=['Entry', 'Options'])

for word in not_matched:
  df_lev.loc[j, 'Entry'] = word
  scores = []
  for s in s2:
    score = nltk.edit_distance(word.lower(), s)
    scores.append((score, s))
  df_lev.loc[j, 'Options'] = str(sorted(scores)[:15])
  j += 1
  print(j, 'entries done')
'''

