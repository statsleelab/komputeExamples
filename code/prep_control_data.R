##################################
# Prep Behavioral phenotype data
##################################
behav.pheno.file = "~/Google Drive Miami/Miami_IMPC/data/v10.1/AllBehaviourControls.Rdata"
load(file=behav.pheno.file) # df
dim(df) # 686399 74
table(df$data_type)
df %>% filter(data_type=="FLOAT") %>% select(data_point) %>% summary()
df %>% filter(data_type=="IMAGE") %>% select(data_point) %>% summary()
df %>% filter(data_type=="INT") %>% select(data_point) %>% summary()
df %>% filter(data_type=="TEXT") %>% select(data_point) %>% summary()

## data filtering - filter out text, image and all rows with data_point == "NA"
dim(df)
df = df %>% filter(data_type!="TEXT" & data_type!="IMAGE") %>% drop_na(data_point)
dim(df)

behav <- df[, c("biological_sample_id","procedure_stable_id","procedure_name","parameter_stable_id","parameter_name",
                "data_point","sex","age_in_weeks","weight","phenotyping_center","date_of_experiment","strain_name",
                "developmental_stage_name","observation_type","data_type")]
behav$parameter_name <- trimws(behav$parameter_name)
save(behav, file="~/Google Drive Miami/Miami_IMPC/data/v10.1/AllBehavControls_small.Rdata")


######################################
# Prep Non-behavioral phenotype data
######################################
nonbehav.pheno.file = "~/Google Drive Miami/Miami_IMPC/data/v10.1/AllNonBehaviourControls.Rdata"
load(file=nonbehav.pheno.file) # df
dim(df) # 17400869 74
table(df$data_type)
df %>% filter(data_type=="FLOAT") %>% select(data_point) %>% summary()
df %>% filter(data_type=="IMAGE") %>% select(data_point) %>% summary()
df %>% filter(data_type=="INT") %>% select(data_point) %>% summary()
df %>% filter(data_type=="TEXT") %>% select(data_point) %>% summary()

## data filtering - filter out text, image and all rows with data_point == "NA"
dim(df)
df = df %>% filter(data_type!="TEXT" & data_type!="IMAGE") %>% drop_na(data_point)
dim(df)

## find out why there are multiple measurements for a biological_ID(mouse id) & phenotype pair.
#IPGTT_bgc <- subset(df, procedure_name=="Intraperitoneal glucose tolerance test (IPGTT)"&parameter_name=="Blood glucose concentration")
#IPGTT_bgc <- IPGTT_bgc[order(IPGTT_bgc$biological_sample_id, IPGTT_bgc$parameter_name),]
#dim(IPGTT_bgc)
#head(IPGTT_bgc, 20)


nonbehav <- df[, c("biological_sample_id","procedure_stable_id","procedure_name","parameter_stable_id","parameter_name",
                   "data_point","sex","age_in_weeks","weight","phenotyping_center","date_of_experiment","strain_name",
                   "developmental_stage_name","observation_type","data_type")]
nonbehav$parameter_name <- trimws(nonbehav$parameter_name)
save(nonbehav, file="~/Google Drive Miami/Miami_IMPC/data/v10.1/AllNonBehavControls_small.Rdata")


## check both datasets are different
#PPI.behav <- subset(behav, procedure_name=="Acoustic Startle and Pre-pulse Inhibition (PPI)")
#table(PPI.behav$parameter_name)
#PPI.nonbehav <- subset(nonbehav, procedure_name=="Acoustic Startle and Pre-pulse Inhibition (PPI)")
#table(PPI.nonbehav$parameter_name)

#########################
# Combine both data sets
#########################
allpheno <- rbind (behav, nonbehav)
save(allpheno, file="~/Google Drive Miami/Miami_IMPC/data/v10.1/AllControls_small.Rdata")
