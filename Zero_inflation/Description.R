library(JM);
library(numDeriv)
library(survminer)
library(numDeriv)
library(survival);library(splines)
library(cplm)
rm(list=ls())
Long=read.table("/Users/taban/Desktop/ZI_DP_Oct/Application/Long.txt",header=TRUE)
Surv=read.table("/Users/taban/Desktop/ZI_DP_Oct/Application/surv.txt",header=TRUE)
attach(Long)
attach(Surv)
names(Surv)
names(Long)
SubjectID
table(SubjectID)
table(SubjectIDs)


names(Surv)
ID=unique(SubjectID)
m=table(SubjectID)
n=dim(Surv)[1]
y=obstime=LibrarySize1=HH=matrix(NA,n,max(m)) 
for(i in 1:n){
  y[i,1:length(Prevotella[SubjectID==ID[i]])]=Prevotella[SubjectID==ID[i]]
  obstime[i,1:length(GWColl[SubjectID==ID[i]])]=GWColl[SubjectID==ID[i]]
  LibrarySize1[i,1:length(LibrarySize[SubjectID==ID[i]])]=LibrarySize[SubjectID==ID[i]]
  HH[i,1:length(History_of_preterm_delivery[SubjectID==ID[i]])]=History_of_preterm_delivery[SubjectID==ID[i]]
  
  
}

History=HH[,1]

#############################################################
library(vcdExtra) ## P vs ZIP
zero.test(Prevotella)
#############################################################
n2=length(SubjectIDs)
n1=length(SubjectID)
z=rep(0,n1)
z[Prevotella==0]=1
Z=1-z
Z[Z==0]="zero"
Z[Z==1]="non-zero"

df <- data.frame(table(Z))
rownames(df)<-c("0","n0")
head(df)

library(ggplot2)
# Basic barplot
p<-ggplot(data=df, aes(x=Z , y=Freq)) +
  geom_bar(stat="identity", color="salmon", fill="khaki1")
p

# Horizontal bar plot
#p + coord_flip()
##########################################
Prevotellabar=Prevotella
Prevotellabar[Prevotella==0]="0"
Prevotellabar[Prevotella>0 & Prevotella <=10]="1-10"
Prevotellabar[Prevotella>10 & Prevotella <21]="10-20"
Prevotellabar[Prevotella>20 & Prevotella <101]="21-100"
Prevotellabar[Prevotella>100 & Prevotella <500]="101-500"
#Prevotellabar[Prevotella>501 & Prevotella <2000]="501-20"
Prevotellabar[Prevotella>500]=">500"



df <- data.frame(table(Prevotellabar))
head(df)

df$Prevotellabar <- factor(df$Prevotellabar,levels = c("0", "1-10","10-20",  "21-100", "101-500",">500"))


# Basic barplot
p<-ggplot(data=df, aes(x=Prevotellabar, y=Freq)) +
  geom_bar(stat="identity", color="darkgoldenrod3", fill="darkgoldenrod2")+
  xlab("Prevotella Abundances")+ylab("Frequency")
p


############################################################

library("survminer")


require("survival")



AA=coxph(Surv(GWDels, Deliverys) ~ Preeclampsias)
summary(AA)
fit <- survfit(Surv(GWDels,event=Deliverys)~1, data=Surv)
ggsurvplot(fit,data=Surv,
           color="darkgoldenrod3",
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme#palette = c("#37004D"),show.legend=FALSE)+
           palette = c("Turquoise"),show.legend=FALSE)+
  xlab("Time to delivery (week)")

