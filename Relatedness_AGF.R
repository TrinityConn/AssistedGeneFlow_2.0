### Import VCF ###
##################

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

# Read in VCF input file.
#vcf <- read.vcfR("D:/PSU/NOAA/PRO100175_PSU175_SAX_b02/Plate9SR10074_allSamples.vcf")
vcf <- read.vcfR("AGF_NoClones.vcf")


######################################
### Convert to Genind and genclone ###
######################################
# Convert VCF file into a genind for the Poppr package.
genind_obj <- vcfR2genind(vcf)

# Add population information to the genind object.
population_info_data_table <- read.table("darpasamples.txt",
                                         check.names=FALSE, header=F, 
                                         na.strings=c("", "NA"), stringsAsFactors=FALSE, sep="\t",quote="")
colnames(population_info_data_table) <- c("sample.id", "Population")



genind_obj@pop <- as.factor(population_info_data_table$region)
strata(genind_obj) <- data.frame(pop(genind_obj))

# Convert genind object to a genclone object.
genind_clone <- as.genclone(genind_obj)


#####################
### Calculate MLGs ###
######################

# Calculate the bitwise distance between individuals.
bitwise_distance <- bitwise.dist(genind_clone)
genDist<-as.matrix(bitwise_distance)

xy<-t(combn(colnames(genDist), 2))
df<-data.frame(xy, dist=genDist[xy])
genDist2<-rownames_to_column(as.data.frame(genDist),"affy_id")
# Multilocus genotypes (threshold of 3.2%).
mlg.filter(genind_clone, distance=bitwise_distance) <- 0.032
m <- mlg.table(genind_clone, background=TRUE, color=TRUE)

# Create list of MLGs.
mlg_ids <- mlg.id(genind_clone)


#PCA
vcf2<-read.vcfR('darpa_final.vcf.gz')
vcf.fn <- "darpa_final.vcf.gz"
snpgdsVCF2GDS(vcf.fn, "pca3.gds", method="biallelic.only")
snpgdsSummary("pca3.gds")
genofile<-snpgdsOpen(filename="pca3.gds", readonly = FALSE)


pca<-snpgdsPCA(genofile, autosome.only = FALSE)


#plot with pop info 
popinfo<-population_info_data_table
colnames(popinfo)<-c('sample.id', 'pop')

sample.id<-read.gdsn(index.gdsn(genofile, "sample.id"))

head(cbind(sample.id, popinfo))

tab<-data.frame(sample.id=pca$sample.id,
                EV1=pca$eigenvect[,1], 
                EV2=pca$eigenvect[,2],
                stringsAsFactors = FALSE)

plot(tab$EV2, tab$EV1, xlab="Eigenvector 2", ylab="Eigenvector 1")

View(tab)

pca2<-merge(tab, popinfo, by=c("sample.id"))


ggplot(pca2, aes(x=EV2, y=EV1, color=pop))+geom_point(size=2)+theme_bw()+xlab("PC2 (6.73 %)")+ ylab("PC1 (12.41 %)")+
  ggtitle('Principal Component Analysis of individuals in AGF 2.0')+labs(color="Population")+stat_ellipse(type='norm', linetype=2)+stat_ellipse(type="t")


#relatedness 

ibd.robust<-snpgdsIBDKING(genofile, autosome.only = FALSE, type=c("KING-robust", "KING-homo"))
class(genofile)

KinshipCoeff<-ibd.robust$kinship
sample.id<-ibd.robust$sample.id
Ibs<-ibd.robust$IBS0

Kinship<-rbind(KinshipCoeff, sample.id, Ibs)

ibd<-snpgdsIBDKING(genofile, sample.id=sample.id, maf=0.05, missing.rate=0.05, autosome.only=FALSE, 
                   type=c("KING-robust", "KING-homo"))
heatmap(KinshipCoeff)
ibd.coeff<-snpgdsIBDSelection(ibd)

head(ibd.coeff)

plot(ibd.coeff$IBS0, ibd.coeff$kinship)

write.csv(ibd.coeff, file="ibd_coeff.csv")
#heterozygosity 




#Neighbor joining tree 


ibs <- snpgdsIBS(genofile, num.thread=2, autosome.only=FALSE)
ibs$sample.id <-gds_data_table_join$user_specimen_id

# Cluster analysis on the genome-wide IBS pairwise distance matrix.
set.seed(100)
ibs.hc <- snpgdsHCluster(ibs)

cols <- piratepal("basel")
set.seed(999)
par(cex=0.6, cex.lab=1, cex.axis=1.8,cex.main=2)
rv <- snpgdsCutTree(ibs.hc, col.list=cols, pch.list=15)
snpgdsDrawTree(rv, main="Color by Cluster", leaflab="perpendicular", yaxis.kinship=FALSE)

abline(h = 0.032, lty = 2)
legend("topleft", legend=levels(rv$samp.group), xpd=T, col=cols[1:nlevels(rv$samp.group)], pch=15, cex=1.2)

# Color cluster by region.
race <- as.factor(popinfo)


rv2 <- snpgdsCutTree(ibs.hc, samp.group=race,col.list=cols, pch.list=15, label.H=TRUE, label.Z=TRUE)
snpgdsDrawTree(rv2, main="Color by Region", leaflab="perpendicular", yaxis.kinship=FALSE)
abline(h = 0.032, lty = 2)
legend("topleft", legend=levels(race), xpd=T, col=cols[1:nlevels(race)], pch=15, ncol=4, cex=1.2)

snpgdsClose(genofile)




###density plot of kinship values 
within<-AGF2_Kinship2%>%
  filter(type=="within")
aov1<-aov(kinship~pop1, within)

ggplot(data=AGF2_Kinship2, aes(x=kinship, group=type, fill=type))+geom_density(adjust=1.5, alpha=.4)+labs(fill="Population")+
  theme_bw()+xlab("Kinship Coefficient")+geom_vline(xintercept = 0, linetype="dashed")+facet_wrap(~pop1)

##make heatmap of kinship values 

##upload kinship data 

kin<-AGF2_kinship_matrix

#make kinship data a heatmap 

data2<-data.matrix(kin)

kin2<-as.matrix(kin)

library(tidyr)
kin3<-pivot_wider(kin, names_from=samp2, values_from=kinship)
View(kin3)

ggplot(CURFLProportions, aes(x=Ancestry, y=Proportion, fill=Ancestry))+geom_boxplot()+theme_bw()
a

##calculating heterozygosities 
tuk<-TukeyHSD(anova1)
het<-AGF_MLG_Report%>%
  group_by(Population)%>%
  summarise(heto=mean(heterozygosity), sd=sd(heterozygosity))

anova1<-aov(heterozygosity~Population, data=AGF_MLG_Report)
vcf2
gen2<-vcfR2genlight(vcf)
library(dartR)
install.packages("dartR")
library(dartR)
pop<-pop_final
colnames(pop)[1]<-"ind"
colnames(pop)[2]<-"popu"

pop_FL<-pop%>%
  filter(popu=="FLxFL")_
pop_fl_names<-pop_FL$ind
pop_C<-pop%>%
  filter(popu=="CURxCUR")
pop(gen2)<-pop$pop
gen2<-gl.compliance.check(gen2)

gl.define.pop(gen1, pop_fl_names, "FLxFL")
gen<-gl.filter.monomorphs(gen1)
gl.report.heterozygosity(gen2, method="ind")
x