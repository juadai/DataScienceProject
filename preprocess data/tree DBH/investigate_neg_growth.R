df = read.csv('/Users/nickjolly/Desktop/tfn/DataScienceProject/datasets/tree_DBH_hollows/output/all_plots.csv')
str(df)

df <- df[df$Tree.Number!='',]
df$Tree.Number <- as.numeric(df$Tree.Number)

df <- df[is.na(df$DBH.year.0)==F,]
df <- df[is.na(df$DBH.year.6)==F,]

df$growth <- df$DBH.year.6 - df$DBH.year.0


neg_trees <- df[df$growth<=0,]

keep = c("Tree.Number", "Plot", "Treatment", "DBH.year.0", "DBH.year.6", "growth")
neg_trees <- subset(neg_trees, select=keep)
neg_trees <- neg_trees[order(neg_trees$Plot),]
rownames(neg_trees) <- 1:length(rownames(neg_trees))

neg_trees
