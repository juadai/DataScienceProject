import pandas as pd

def read(filename):
    df = pd.read_excel(filename, sheet_name=None, engine="openpyxl")
    for name, sheet in df.items(): # iterate over all sheets and contents
        print(name)
        print(sheet.head(3))

read('Gannawarra.xlsx')