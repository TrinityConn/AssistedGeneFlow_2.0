###AGF_Genotyping_Results##

#Required R-packages for multi-locus genotype calling
library(adegenet)
library(ape)
library(data.table)
library(dbplyr)
library(dplyr)
library(ggplot2)
library("knitr")
library("maps")
library("mapproj")
library(poppr)
library(RColorBrewer)
library("SNPRelate")
library("tidyr")
library(vcfR)
library(vegan)
library(yarrr) #color palettes
install.packages("BiocManager")
BiocManager::install("SNPRelate")
theme_set(theme_bw())


#import_vcf
e
vcf <- read.vcfR("darpa_final.vcf")

# Convert VCF file into a genind for the Poppr package.
genind_obj <- vcfR2genind(vcf)

# Add population information to the genind object.
population_info_data_table <- read.table("popinfo.txt",
                                         check.names=FALSE, header=T, 
                                         na.strings=c("", "NA"), stringsAsFactors=FALSE, sep="\t",quote="")
colnames(population_info_data_table) <- c("affy_id","user_specimen_id","tube_id", "cross")

colnames(AGF_Norepeatsamples)<-c("affy_id")

pop_info<-subset(population_info_data_table, affy_id %in% AGF_Norepeatsamples$affy_id)

genind_obj@pop <- as.factor(evenpop$cross)
strata(genind_obj) <- data.frame(pop(genind_obj))

# Convert genind object to a genclone object.
genind_clone <- as.genclone(genind_obj)

#################################
####Filter duplicate samples#####
#################################





######################
### Calculate MLGs ###
######################

# Calculate the bitwise distance between individuals.
bitwise_distance <- bitwise.dist(genind_clone)

genDist<-as.matrix(bitwise_distance)

xy<-t(combn(colnames(genDist), 2))
df<-data.frame(xy, dist=genDist[xy])

# Multilocus genotypes (threshold of 3.2%).
mlg.filter(genind_clone, distance=bitwise_distance) <- 0.032
m <- mlg.table(genind_clone, background=TRUE, color=TRUE)

# Create list of MLGs.
mlg_ids <- mlg.id(genind_clone)


mlg_ids_data_table <- data.table(mlg_ids, keep.rownames=TRUE)
# Rename the mlg_ids column.
setnames(mlg_ids_data_table, c("mlg_ids"), c("affy_id"))

sample_mlg_tibble <- mlg_ids_data_table %>%
  dplyr::group_by(row_number()) %>%
  dplyr::rename(group="row_number()") %>%
  unnest (affy_id) %>%
  # Join with mlg table.
  left_join(population_info_data_table %>%
              select("affy_id", "cross"),
            by="affy_id")


###############################
###Calculate Heterozygosity####
###############################

  # Heterozygous alleles of all SNPs.
  gt <- extract.gt(vcf, element="GT", as.numeric=FALSE)
  heterozygous_alleles_all <- apply(gt, MARGIN=2, function(x) {sum(lengths(regmatches(x, gregexpr("0/1", x))))})
  heterozygous_alleles_all <- (heterozygous_alleles_all / nrow(gt)) * 100
  heterozygous_alleles_all_data_frame <- data.frame(heterozygous_alleles_all)
  heterozygous_alleles_all_data_table <- setDT(heterozygous_alleles_all_data_frame, keep.rownames=TRUE)[]
View(heterozygous_alleles_all_data_table)  



het<-setnames(heterozygous_alleles_all_data_table, c("affy_id"), c("Heterozygosity"))


############################################
####Combine Heterozygosity and Clonal IDS###
############################################
colnames(heterozygous_alleles_all_data_table)<-c("affy_id", "heterozygosity")
colnames(sample_mlg_tibble)<-c("affy_id", "Clonal_Group", "Population")

agf_report<-merge(sample_mlg_tibble, heterozygous_alleles_all_data_table, by=("affy_id"))

View(agf_report)

sampleinfo<-population_info_data_table[c("affy_id", "tube_id")]

write.csv(agf_report, "AGF_MLG_Report.csv")
agf_report<-merge(agf_report, sampleinfo, by=c("affy_id"))


cur<-agf_report%>%
  filter(Population=="CURxCUR")

pr<-agf_report%>%
  filter(Population=="CURxPR")


write.csv(cur, file="CUR_individuals.csv")
write.csv(pr, file="PR_individuals.csv")


############################################
####### Calculate Relatedness ##############
############################################

vcf.fn <- "AGF_NoRepeats.vcf"
snpgdsVCF2GDS(vcf.fn, "AGF.gds", method="biallelic.only")
genofile <- snpgdsOpen(filename="AGF.gds", readonly=FALSE)


ibd<-snpgdsIBDKING(genofile, autosome.only = FALSE, type=c("KING-robust", "KING-homo"))


ibd<-ibd$kinship
names<-ibd.robust$sample.id
names<-youngcolonies_id
rownames(ibd)<-names
colnames(ibd)<-names
closefn.gds(genofile)

ibddist<-as.matrix(ibd)
xy<-t(combn(colnames(ibd), 2))
df<-data.frame(xy, dist=ibd[xy])


write.csv(kin,  file="AGF2_Kinship2.csv")


kin<-AGF2_Kinship2

##create new column 

kin$relate<-paste(kin$pop1, kin$pop2)

ggplot(kin, aes(x=pop1, y=kinship, fill=type))+geom_boxplot()+geom_hline(yintercept=0, linetype="dashed")+theme_bw()+xlab("Cross")+ggtitle("Relatedness within and between Admixed and Non-Admixed Populations")


## heat map of kinship 

heatmap(ibd)


##try it in ggplot 
#convert to tidy data 

ibd_2<-ibd%>%
  as.data.frame()%>%
  rownames_to_column("indiv1")%>%
  pivot_longer(-c(indiv1), names_to="indiv2", values_to="kinship")

ggplot(ibd_2, aes(x=indiv1, y=indiv2, fill=kinship))+geom_raster()+scale_fill_viridis_c()
###############################################
#########Principal Component Analysis##########
###############################################

###plotting pca of SNP data 



#plot pca
pca<-snpgdsPCA(genofile,autosome.only = FALSE, num.thread = 2)

#variance proportion
pc.percent<-pca$varprop*100
head(round(pc.percent,2))

#make a data frame

tab<-data.frame(affy_id=pca$sample.id,
                EV1=pca$eigenvect[,1],
                EV2=pca$eigenvect[,2],
                stringsAsFactors = FALSE)

plot(tab$EV2, tab$EV1, xlab="Eigenvector2", ylab="Eigenvector 1")
pca_pop<-merge(tab, pop_info, by=c("affy_id"))

View(pca_pop)
ggplot(pca_pop, aes(x=EV2, y=EV1, color=cross))+geom_point(size=2)+theme_bw()+xlab("PC2 (5.27%)")+
  ylab("PC1(12.12%)")+
  ggtitle("PCA AGF2.0")+
  labs(color="Population")+stat_ellipse(type='norm', linetype=2)+stat_ellipse(type="t")

  

############################################################
#####Within and Between Population Diversity Statistics#####
############################################################

ggplot(AGF2_Kinship2, aes(x=pop1, y=kinship, fill=type))+geom_violin()+ggtitle("kinship values within each population")+theme_bw()+geom_hline(yintercept = 0, linetype="dashed")
##FL x FL ##


#subset FLxFL 

FL<- 


  
  
#############################################################
####### Plotting Admixture Analysis ##########################
#############################################################


##load in required tools 
library(devtools)
#install_github('royfrancis/pophelper')
library(pophelper)

afiles<-list.files(pattern="*.Q", full.names=T)



alist<-readQ(files=afiles, indlabfromfile=T)

#load in pop file 
popfile<-Fam_CEL_pop
#change to character type 
popfile$V1<-as.character(popfile$V1)
popfile$Pop<-popfile$V1
popfile<-subset(popfile, select=-c(V1))

#repeat for just fam file
afilesfam<-list.files(path="~/Documents/PSU/AGF_Project/AGF_2.0/SNP_Plates/VCF_Files/", 
                      pattern="*Q", full.names=T)
alistfam<-readQ(files=afilesfam, indlabfromfile=T)

#loadinpopfile
popfile<-AGF2_pop
popfile<-popfile[popfile$V1=="FloridaParent"]<-"FL"
popfile<-popfile[popfile$V1=="CuracaoParent"]<-"C"

outputlabs<-c("FamPlot_Final1", "FamPlot_Final2", "FamPlot_Final3", "FamPlot_Final4", "FamPlot_final5", "FamPlot_final6")
spnames<-c("K=2", "K=3", "K=4", "K=5", "K=6", "K=7")
outputlabs2<-c("DPlot_1", "Dplot_2", "Dplot_3", "Dplot_4", "Dplot_5", "Dplot_6")
titlea<-c("Admixture Results of Assisted Gene Flow Recruits", "Admixture Results of Assisted Gene Flow Recruits", "Admixture Results of Assisted Gene Flow Recruits", "Admixture Results of Assisted Gene Flow Recruits", "Admixture Results of Assisted Gene Flow Recruits", "Admixture Results of Assisted Gene Flow Recruits")
#plot files of diagnostic loci
a<-plotQ(alist, imgoutput="sep", showindlab=T, grplab=popfile,
         showtitle=F, showsubtitle=F, 
         , height=1.6, indlabsize=1.4, indlabheight=.08, 
        ,grplabsize=1.4, subtitlesize=4, indlabspacer=-1, barbordercolour="white")

#plot files of just fam individuals with full genotyping set 
b<-plotQ(alistfam, imgoutput="sep", showindlab=T,
         showtitle=F, showsubtitle=F, titlelab=titlea, height=1.6, indlabsize=1.4, splab=spnames, 
         indlabheight=.08, grplabsize=1.4, indlabspacer=-1, barbordercolour="white", subtitlesize = 4, 
         outputfilename=outputlabs, imgtype="png")

library(cowplot)

plot_grid(a,b)


f<-as.dataframe('12'13'')


##load in png files 
library(png)
install.packages("EBImage")
library(EBImage)
png<-readPNG(system.file('Dplot_2.png'), package = png)



###representing families###


##load in kinship values ##

kin<-AGF2_Kinship2

##assign relationships 

sibs<-kin%>%
  filter(kinship>=.18)

sibs<-sibs%>%
  filter(kinship<=.4)

write.csv(sibs, file="AGF_Sib.csv")



AGF2_Kinship2

mean<-within%>%
  group_by(pop1)%>%
  summarize(mean=mean(kinship))


V##merge in gendist information 