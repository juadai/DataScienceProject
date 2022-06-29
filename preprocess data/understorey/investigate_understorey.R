library(readxl)
df <- read_excel('/Users/nickjolly/Documents/2021:22/Unimelb/MAST90106/TfN Dataset/Ecothining_understoreydata2020.xlsx')
species <- unique(df$T002_Flora_Species_name)

# names(df)
# keep = c("Year_monitoring", "Date_monitoring", "Plot_number", "Plot_treatment", "Quadrat_number", 
#          "Quardat_fenced", "Quadrat_gap_or_no", "Record_ID", "T002_Flora_Species_name", "Score")
# df <- subset(df, select=keep)


