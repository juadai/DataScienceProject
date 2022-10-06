import os
import pandas as pd

def read(filename):
    df = pd.read_excel(filename, sheet_name=None, engine="openpyxl")
    for name, sheet in df.items():
        if name.lower() == 'final deliverable':
            sheet.rename(columns={ sheet.columns[0]: "Statistic" }, inplace = True)
            sheet.rename(columns = lambda x: x.strip(), inplace = True)

            idx_parcel = sheet.index[sheet['Statistic'].str.contains("parcel statistics", na=False, case=False)].tolist()
            idx_veg = sheet.index[sheet['Statistic'].str.contains("Native vegationtion per parcel", na=False, case=False)].tolist()
            idx_covenant = sheet.index[sheet['Statistic'].str.contains("covenant statistics", na=False, case=False)].tolist()
            idx_empty_rows = sheet[sheet.isnull().all(axis=1)].index.to_list()
            end_of_first_block = idx_empty_rows[0]

            return sheet.loc[:, 'Statistic':'total study area'], idx_parcel[0], idx_veg[0], idx_covenant[0], end_of_first_block

def concat_df(lga_name, df, new_data, idx_parcel, idx_veg, idx_covenant, end_of_first_block):
    truncated_data = new_data.iloc[idx_parcel:end_of_first_block, :]

    col_names = ['LGA'] + new_data.columns.to_list()
    lga_row = pd.DataFrame(columns=col_names, index=range(1))
    lga_row.iloc[0,0] = lga_name

    new_df = pd.concat([lga_row, truncated_data], ignore_index=True)
    if df is not None:
        new_df = pd.concat([df, new_df], ignore_index=True)
    return new_df

def get_lga_name(filename):
    return os.path.basename(filename).split('.')[0]

def write_to_excel(df):
    df.rename(columns={'bioregion': 'Bioregion', 
                        'sample size': 'Sample Size',
                        'min': 'Min',
                        'max': 'Max',
                        'mean': 'Mean',
                        'median': 'Median',
                        'total study area': 'Total Study Area'}, inplace=True)
    
    df.to_excel("output.xlsx", sheet_name = "Final Deliverable") 

def main():
    list_subfolders_paths = [f.path for f in os.scandir('.') if f.is_dir()]
    # print(list_subfolders_paths)

    if len(list_subfolders_paths) > 0:
        first_path = list_subfolders_paths[0]
        for f in os.scandir(first_path):
            if f.is_file():
                lga_name = get_lga_name(f)
                new_data, idx_parcel, idx_veg, idx_covenant, end_of_first_block = read(f)
                final_df = concat_df(lga_name, None, new_data, idx_parcel, idx_veg, idx_covenant, end_of_first_block)

        if len(list_subfolders_paths) > 1:
            for path in list_subfolders_paths[1:]:
                for f in os.scandir(path):
                    if f.is_file():
                        lga_name = get_lga_name(f)
                        new_data, idx_parcel, idx_veg, idx_covenant, end_of_first_block = read(f)
                        final_df = concat_df(lga_name, final_df, new_data, idx_parcel, idx_veg, idx_covenant, end_of_first_block)
    
    write_to_excel(final_df)

if __name__ == "__main__":
    main()