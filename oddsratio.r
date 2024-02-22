#Odds ratio and correlation analysis using the "GeneOverlap" package 

library('GeneOverlap')
# BiocManager::install("GeneOverlap")


#Comparing differentially expressed genes that have a p <0.05  
#opening the table and extracting column subset through the $
AMPMDE = read.csv('./bmal1/New_RNAseq_020623/oddsratio/files/fdr_files/AMWTvPMWT_new.csv', header=T)
AMPMUP = AMPMDE[AMPMDE$log2FoldChange>0,] 
AMPMDOWN = AMPMDE[AMPMDE$log2FoldChange<0,] 


SNADE = read.csv('./bmal1/New_RNAseq_020623/oddsratio/files/fdr_files/SNAvsSham.csv', header=T) 
SNAUP = SNADE[SNADE$log2FoldChange>0,] 
SNADown = SNADE[SNADE$log2FoldChange<0,]


SNAEEDE = read.csv('./bmal1/New_RNAseq_020623/oddsratio/files/fdr_files/SNAEEvsSH.csv', header=T) 
SNAEEUP = SNAEEDE[SNAEEDE$log2FoldChange>0,] 
SNAEEDown = SNAEEDE[SNAEEDE$log2FoldChange<0,]


EEDE = read.csv('./bmal1/New_RNAseq_020623/oddsratio/files/fdr_files/EEvsSH.csv', header=T) 
EEUP = EEDE[EEDE$log2FoldChange>0,] 
EEDown = EEDE[EEDE$log2FoldChange<0,]


DCADE=read.csv('./bmal1/New_RNAseq_020623/oddsratio/files/fdr_files/DCAvLAM.csv', header=T)
DCAUP=DCADE[DCADE$log2FoldChange>0,]
DCADown=DCADE[DCADE$log2FoldChange<0,]


IFDE = read.csv('./bmal1/New_RNAseq_020623/oddsratio/files/fdr_files/IFvAL.csv', header=T) 
IFUP = IFDE[IFDE$log2FoldChange>0,] 
IFDown = IFDE[IFDE$log2FoldChange<0,]

LTDE = read.csv('./bmal1/New_RNAseq_020623/oddsratio/files/fdr_files/lithiumvscontrol.csv', header=T) 
LTUP = LTDE[LTDE$log2FoldChange>0,] 
LTDown = LTDE[LTDE$log2FoldChange<0,]

#Creating the lists

#My list of interest 
CircadianList = list(AMPMDE=AMPMDE$gene_name,
                    AMPMUP=AMPMUP$gene_name,
                    AMPMDOWN=AMPMDOWN$gene_name)
                     
CircadianList = list(AMPMUP=AMPMUP$gene_name,
                     AMPMDOWN=AMPMDOWN$gene_name)
#Second list
DatasetList = list(LTUP=LTUP$gene_name,
                   LTDown=LTDown$gene_name,
                   SNAUP=SNADown$gene_name,
                   SNADown=SNADE$gene_name,
                   SNAEEUP=SNAEEUP$gene_name,
                   SNAEEDown=SNAEEDown$gene_name,
                   EEUP=EEUP$gene_name,
                   EEDown=EEDown$gene_name,
                   IFUP=IFUP$gene_name,
                   IFDown=IFDown$gene_name,
                   DCAUP=DCAUP$gene_name,
                   DCADown=DCADown$gene_name)

#DE only
DatasetList = list(LTDE=LTDE$gene_name,
                   SNADE=SNAUP$gene_name,
                   SNAEEDE=SNAEEDE$gene_name,
                   EEDE=EEDE$gene_name,
                   IFDE=IFDE$gene_name,
                   DCADE=DCADE$gene_name)

#total
DatasetList = list(LTDE=LTDE$gene_name,
                   LTUP=LTUP$gene_name,
                   LTDown=LTDown$gene_name,
                   SNADE=SNAUP$gene_name,
                   SNAUP=SNADown$gene_name,
                   SNADown=SNADE$gene_name,
                   SNAEEDE=SNAEEDE$gene_name,
                   SNAEEUP=SNAEEUP$gene_name,
                   SNAEEDown=SNAEEDown$gene_name,
                   EEDE=EEDE$gene_name,
                   EEUP=EEUP$gene_name,
                   EEDown=EEDown$gene_name,
                   IFDE=IFDE$gene_name,
                   IFUP=IFUP$gene_name,
                   IFDown=IFDown$gene_name,
                   DCADE=DCADE$gene_name,
                   DCAUP=DCAUP$gene_name,
                   DCADown=DCADown$gene_name)


#Now we visualise the lists and check them. 
sapply(CircadianList, length)
sapply(DatasetList, length)


#fixing the background size
#Insert background of your largest list 
# in this case it was AMPM
numGenes=23525


#Now create the GOM object 
gom.obj=newGOM(CircadianList,DatasetList,numGenes)
pdf('./bmal1/New_RNAseq_020623/oddsratio/AMPM_split_FDR.pdf',height=40,width=40)
drawHeatmap(gom.obj)
dev.off()
gom.obj
#This will now create a png file that shows your comparisons



