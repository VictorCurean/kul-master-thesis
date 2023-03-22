setwd("C:\\Users\\curea\\Documents\\Thesis Code\\kul-master-thesis\\kul-master-thesis\\kul-thesis-ab")

aff_summary_1 = as.data.frame(read.csv("datasets_aff//20230321_0047670//20230321_0047670_summary.tsv", header = T, sep="\t"))
aff_summary_2 = as.data.frame(read.csv("datasets_aff//20230321_0975319//20230321_0975319_summary.tsv", header = T, sep="\t"))

df <- rbind(aff_summary_1, aff_summary_2)

df = df[!duplicated(df$pdb),]

plot(log(df$affinity), df$delta_g)

hist(df$delta_g)

write.csv(df, "datasets_aff\\Summary.tsv", row.names = F)