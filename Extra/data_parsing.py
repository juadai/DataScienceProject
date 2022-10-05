import os
import pandas as pd

def read(filename):
    df = pd.read_excel(filename, sheet_name="Final Deliverable", engine="openpyxl")
    print(df.head(3))


def main():
    list_subfolders_paths = [f.path for f in os.scandir('.') if f.is_dir()]
    print(list_subfolders_paths)

    for path in list_subfolders_paths:
        for f in os.scandir(path):
            if f.is_file():
                read(f)

main()