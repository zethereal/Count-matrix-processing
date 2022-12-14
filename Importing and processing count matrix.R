#importing and processing RNA seq count matrix

#############FILE importing, processing, and cleaning

setwd("C:/Users/zkart/Desktop/Working_folder")

#list files in working directory
list.files()

#import .txt dataset into new object 
dataset <- read.table("GSE153960.txt", header = T, sep = '')

#importing meta data labels
#from GEO -> SRA run selector -> metadata
meta <- read.delim("SraRunTable.txt",header=TRUE,sep=',')

library(tidyverse) 

#make new column that sums existing columns
ds2  <-  rowSums(dataset[ , 2:1660])  

#bind sum column to existing data frame
dataset2 <- cbind(dataset, ds2)

#remove genes that sum below 10,000 reads across 1660 samples
dataset3 <- subset(dataset3, ds2 > 10000) 

#now that genes have been filtered
#making tidy data with pivot_longer
tdata <- pivot_longer(dataset3, starts_with("CGND"), names_to = "sample", values_to = "counts")

#subsetting metadata table for only few columns
source_key <- meta %>% select(sample_id_alt,source_name,Group)
#'.' is separator in dataset, '-' in in metadata

#substituting out '.' character for '-'
source_key$sample_id_alt <- gsub('-', '.', source_key$sample_id_alt)

#adding new column that specifies class of disease
key <- source_key %>% 
  mutate(disease = case_when(str_detect(Group, "Other" ) ~ "other",
                             str_detect(Group, "ALS" ) ~ "ALS",
                             str_detect(Group, "Control" ) ~ "control"))

#renaming columns
colnames(key) <- c("sample", "source", "group", "status")

colnames(tdata) <- c("gene", "total", "sample", "counts")

#merging two data frames by 'sample' column
merge <- left_join(tdata, key) 

#removing NA value rows
merged <- merge %>% drop_na()

#removing rows that contain 'Other' disease in status column
merged <- merged[- grep("other", merged$status),]

#removing decimals following gene ID
merged$gene <- sub("*\\.[0-9]", "", merged$gene)

#subsetting tdataframe
merge_clean <- select(merged, gene, counts, source, status)

#pivot_wider
mc2 <- merge_clean %>% pivot_wider(names_from = source, values_from = counts)

#removing totaled column
ds <- dataset3[,-1661]

#removing decimals following gene ID
ds$EnsemblID <- sub("*\\.[0-9]", "", ds$EnsemblID)

#replacing '-' with '.' 
meta$sample_id_alt <- gsub('-', '.', meta$sample_id_alt)

#making columns of meta table into char
meta_samp_id <- meta_c$sample_id_alt

#colnames of dataset
data_samp_id <- colnames(ds)

#dropping replicates in metadata by specific columns
meta_uni <- meta2[!duplicated(meta2$sample_id_alt), ]

#finding unique labels
unique <- dplyr::setdiff(meta_samp_id, data_samp_id)
unique2 <- dplyr::setdiff(data_samp_id,meta_samp_id)

#repeating for ALS vs control groups
uni_als_meta_id <- coldata_comp$sample_id_alt

#select multiple columns by name
#small example
counts_sub <- counts[,c("CGND.HRA.01123", "CGND.HRA.01124")]

counts.als <- counts[,c("insert many colnames")]

#remove rows that contain unique strings
#many indicates ~200 unique samples that were not present in both
meta_c <- meta_uni %>% filter(!grepl("CGND.HRA.01132| many", sample_id_alt))

#remove columns that contain unique strings
ds_c <- subset(ds, select=-c(CGND.HRA.00440,CGND.HRA.00441,CGND.HRA.00442,CGND.HRA.00443,
                             CGND.HRA.00444,CGND.HRA.00445,CGND.HRA.00446,CGND.HRA.00447,
                             CGND.HRA.00448,CGND.HRA.00449,CGND.HRA.00450,CGND.HRA.00451,
                             CGND.HRA.00452,CGND.HRA.00454,CGND.HRA.00455,CGND.HRA.00456,
                             CGND.HRA.00457,CGND.HRA.00458,CGND.HRA.00459,CGND.HRA.01504))

#reordering dataset object using meta character object
dsc2 <- ds[, meta_samp_id]

#adding gene ID column back to new DF
dsc <- cbind(dsc2,ds_c$EnsemblID)

#reordering gene ID column to 1st column
dsc <- dsc[, c(1640,1:1639)]

#renaming gene ID column
colnames(dsc2)[1] <- "gene"

#reordering sample_alt_id column to 1st column
metadata <- meta_c[, c(24,1:23,25:30)]

#making col1 into rownames
rownames(metadata) <- metadata$sample_id_alt
rownames(dsc2) <- dsc2$gene

#removing first column of countsdata
colData <- metadata

counts <- dsc2[,-1]

#do row names in coldata match column names in counts data
all(colnames(counts) %in% rownames(coldata))
#are they in the same order?
all(colnames(counts) == rownames(coldata))

length(colnames(counts))
length(rownames(colData))

#remove first row of coldata to match missing id20
coldata <- metadata[-1,]


###############
##############testing combining and filtering through common names

#creating a new object that specifies control or disease states
key <- source_key %>%
  mutate(disease= if_else(str_detect(Group, "Control"), "control", "disease"))

#testing binary classification of disease state based on Group column
key <- key$treatment= if_else(str_detect(Group, "ALS"), "control", "disease")

#replacing description with 1 word class in Group column
key <- key %>% 
  mutate(disease = case_when(str_detect(Group, "Other" ) ~ "other",
                             str_detect(Group, "ALS" ) ~ "ALS",
                             str_detect(Group, "Control" ) ~ "control"))

#key2 removes duplicate rows
key <- unique(key)

#renaming columns
colnames(key) <- c("sample", "source", "group", "disease", "status")

#replacing library name
key$sample<-gsub("CGND.HRA.","",as.character(key$sample))

#remove rows containing string
key <- key[- grep("internal", key$sample),]

#convert sample column to numeric
key$sample <- as.numeric(key$sample)

#drop NAs
key_clean <- key %>% drop_na()
#so now
any(is.na(key)) #returns TRUE
any(is.na(key_clean)) # returns FALSE

#removing decimals following gene ID
dataset$EnsemblID <- sub("*\\.[0-9]", "", dataset$EnsemblID)

#removing prefix of sample name on counts matrix
tdata_clean <- tdata
tdata_clean$sample<-gsub("CGND.HRA.","",as.character(tdata_clean$sample))
tdc2 <- tdata_clean
tdc2$sample <- as.numeric(tdc2$sample)

#join two data frames
tdata_cdf <- as.data.frame(tdc2)

#combining the metadata with the counts matrix
merge <- left_join(tdata_cdf,key_clean)

#dropping na containing rows
merged <- merge %>% drop_na()

#making dataset column names into data frame, then removing the first entry
data_id <- colnames(dataset)
data_id_df <- as.data.frame(data_id)
data_id_df <- (data_id_df[-1,])
data_id_df <- as.data.frame(data_id_df)

samp_id_df <- as.data.frame(samp_id)
#tissue classes
source_name_list <- key$source_name

#finding unique values between two columns
id_unique <- unique(c(colnames(dataset,key2$sample_id_alt)))

samp_id <- key2$sample_id_alt

#unique
unique_ids <- dplyr::setdiff(samp_id_df, data_id_df)
unique_ids <- unique(samp_id_df, data_id_df)

#making tidy data with pivot_longer
tdata <- pivot_longer(dataset, starts_with("CGND"), names_to = "sample", values_to = "counts")

#extremely large dataset, did not complete
merge <- merge(tdata, key2, by="sample")




#can pass these two objects COUNTS and COLDATA onto DESeq2 for differential gene expression analysis

