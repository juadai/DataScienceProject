import os
import pandas as pd

def read(filename):
    df = pd.read_excel(filename, sheet_name="Final Deliverable", engine="openpyxl")
    df.rename(columns={ df.columns[0]: "Statistic" }, inplace = True)
    df.rename(columns = lambda x: x.strip(), inplace = True)

    idx_parcel = df.index[df['Statistic'].str.contains("parcel statistics", na=False, case=False)].tolist()
    idx_veg = df.index[df['Statistic'].str.contains("Native vegationtion per parcel", na=False, case=False)].tolist()
    idx_covenant = df.index[df['Statistic'].str.contains("covenant statistics", na=False, case=False)].tolist()
    idx_empty_rows = df[df.isnull().all(axis=1)].index.to_list()
    end_of_first_block = idx_empty_rows[0]

    return df.loc[:, 'Statistic':'total study area'], idx_parcel[0], idx_veg[0], idx_covenant[0], end_of_first_block

def concat_df(lga_name, df, new_data, idx_parcel, idx_veg, idx_covenant, end_of_first_block):
    temp = new_data.iloc[idx_parcel:end_of_first_block, :]

    col_names = new_data.columns.to_list()

    new_df = pd.DataFrame(columns=col_names, index=range(1))
    print(new_df)

    return df


def main():
    list_subfolders_paths = [f.path for f in os.scandir('.') if f.is_dir()]
    print(list_subfolders_paths)

    result_df = []
    for path in list_subfolders_paths:
        for f in os.scandir(path):
            if f.is_file():
                lga_name = os.path.basename(f).split('.')[0]
                new_data, idx_parcel, idx_veg, idx_covenant, end_of_first_block = read(f)
                result_df = concat_df(lga_name, result_df, new_data, idx_parcel, idx_veg, idx_covenant, end_of_first_block)

main()