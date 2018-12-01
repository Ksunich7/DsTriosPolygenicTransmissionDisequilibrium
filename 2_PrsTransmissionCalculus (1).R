
rm(list=ls(all=TRUE))
getwd()

setwd("C:/Users/Ksuha/Documents/set")

ListOfFiles = list.files();

for (CounterOfFiles in 1:length(ListOfFiles))
{  
  GWAS = read.table(ListOfFiles[CounterOfFiles], sep = '\t', header = TRUE, quote = '', comment.char = "")
  FamilyId = gsub(".genotypes",'',ListOfFiles[CounterOfFiles])
  
  ### change names of columns
  VecNames = names(GWAS)
  for (i in 1:length(VecNames))
    { # i = 2
    if (length(grep("FatherGenotype",VecNames[i]))) {FatherId = gsub("FatherGenotype.",'',VecNames[i]); FatherId = gsub("\\.",'',FatherId); VecNames[i] = 'FatherGenotype';}
    if (length(grep("MotherGenotype",VecNames[i]))) {MotherId = gsub("MotherGenotype.",'',VecNames[i]); MotherId = gsub("\\.",'',MotherId); VecNames[i] = 'MotherGenotype';}
    if (length(grep("OffsprinGenotype",VecNames[i]))) {OffsprinId = gsub("OffsprinGenotype.",'',VecNames[i]); OffsprinId = gsub("\\.",'',OffsprinId); VecNames[i] = 'OffsprinGenotype';}
    }
  
  GWAS$Summary = paste(GWAS$FatherGenotype,GWAS$MotherGenotype,GWAS$OffsprinGenotype,GWAS$A1.major.allele.,GWAS$A2.minor.allele., sep = '_')
  GWAS$Summary = gsub('/','',GWAS$Summary)

  Prs1 = function(x)
  {
  temp = unlist(strsplit(as.character(x),'_'))
  AssesedAllele = as.character(temp[4])
  Father = unlist(strsplit(as.character(temp[1]),''))
  Mother = unlist(strsplit(as.character(temp[2]),''))
  Kid = unlist(strsplit(as.character(temp[3]),''))
  NumberOfAccessedAllelesInFather = length(Father[Father %in% AssesedAllele])  # 0, 1 or 2
  NumberOfAccessedAllelesInMother = length(Mother[Mother %in% AssesedAllele])  # 0, 1 or 2
  NumberOfAccessedAllelesInKid = length(Kid[Kid %in% AssesedAllele])  # 0, 1 or 2
  res = paste(NumberOfAccessedAllelesInFather,NumberOfAccessedAllelesInMother,NumberOfAccessedAllelesInKid,sep = '_')
  return(res)
  }
  GWAS$res = apply(as.matrix(GWAS$Summary), 1, FUN = Prs1)
  Extract = function(x) {unlist(strsplit(as.character(x),'_'))[1]}; GWAS$NumberOfAccessedAllelesInFather = apply(as.matrix(GWAS$res), 1, FUN = Extract)
  Extract = function(x) {unlist(strsplit(as.character(x),'_'))[2]}; GWAS$NumberOfAccessedAllelesInMother = apply(as.matrix(GWAS$res), 1, FUN = Extract)
  Extract = function(x) {unlist(strsplit(as.character(x),'_'))[3]}; GWAS$NumberOfAccessedAllelesInKid = apply(as.matrix(GWAS$res), 1, FUN = Extract)

  GWAS$NumberOfAccessedAllelesInFather = as.numeric(GWAS$NumberOfAccessedAllelesInFather)*GWAS$Beta
  GWAS$NumberOfAccessedAllelesInMother = as.numeric(GWAS$NumberOfAccessedAllelesInMother)*GWAS$Beta
  GWAS$NumberOfAccessedAllelesInKid = as.numeric(GWAS$NumberOfAccessedAllelesInKid)*GWAS$Beta

  PrsFather = sum(GWAS$NumberOfAccessedAllelesInFather)
  PrsMother = sum(GWAS$NumberOfAccessedAllelesInMother)
  PrsKid = sum(GWAS$NumberOfAccessedAllelesInKid)

  OneLine = data.frame(FamilyId,FatherId,MotherId,OffsprinId,PrsFather,PrsMother,PrsKid)  # add - family ID and GWAS name
  if (CounterOfFiles == 1) {Final = OneLine}
  if (CounterOfFiles >  1) {Final = rbind(Final,OneLine)}
}

Final$MidParent = (Final$PrsFather+Final$PrsMother)/2
Final$MidParentMinusKid = Final$MidParent - Final$PrsKid
summary(Final$MidParentMinusKid)
t.test(Final$MidParent,Final$PrsKid, paired = TRUE)
wilcox.test(Final$MidParent,Final$PrsKid, paired = TRUE)
wilcox.test(Final$MidParentMinusKid, mu = 0)
hist(Final$MidParentMinusKid, breaks = 100)



