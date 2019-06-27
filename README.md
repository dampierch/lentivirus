# Lentivirus RNA-seq project

## Codes

### Samples
320542    218L Par
320543    218L CC1 LV4
320544    218L CC2 LV9
320545    218L SCR LV8
320546    210R Par
320547    210R CC1 LV1
320548    210R CC2 LV1
320549    210R SCR LV1
320550    129L Par
320551    129L CC1 LV1
320552    129L CC2 LV1
320553    129L SCR LV1

```
# found 144 files in master vs 38 in Lenti
cd /nv/vol326/cphg_caseylab/FirstExposures/casey_grc_rnaseq_2/releaseFiles
for i in {320542..320553}; do
  idx=`find . -type f -name "${i}*"`
  for j in ${idx}; do
    cp $j ../../../lentivirus/raw_data
  done
done
```

### Glossary
```
# set vars
nwgc_ep=178d2980-769b-11e9-8e59-029d279f7e24
uva_main=uva#main-DTN

# login and get uva#main-DTN endpoint id
globus login
globus endpoint search ${uva_main}
uva_ep=c4d80096-7612-11e7-8b5e-22000b9923ef
globus ls ${uva_ep}:~/nv/vol326/cphg_caseylab/nwgc/

# confirm nwgc endpoint id
globus endpoint search ${nwgc_ep}

# exlore nwgc endpoint and transfer glossary
globus ls ${nwgc_ep}:~
globus ls ${nwgc_ep}:~/casey_grc_rnaseq_2/
globus transfer ${nwgc_ep}:~/casey_grc_rnaseq_2/RNASeq_Files_Glossary.pdf ${uva_ep}:~/nv/vol326/cphg_caseylab/nwgc/RNASeq_Files_Glossary.pdf --label "CLI_nwgc_glossary"
```

### Genes

__COLCA 1__
* ENSG00000196167
* ENST00000620864.1
* ENST00000540738.3
* ENST00000355430.4
* ENST00000532918.4
* ENST00000526150.1

__COLCA 2__
* ENSG00000214290
* ENST00000398035.6
* ENST00000526216.1
* ENST00000614153.4
* ENST00000610738.5
* ENST00000638573.1
* ENST00000528846.5
* ENST00000639470.1

## Counts extraction
id    218L Par    218L CC1 LV4    218L CC2 LV9    218L SCR LV8    210R Par    210R CC1 LV1    210R CC2 LV1    210R SCR LV1    129L Par    129L CC1 LV1    129L CC2 LV1    129L SCR LV1
ENSG00000196167    204    97    312    172    91    54    102    100    162    151    344    109
ENSG00000214290    447    196    186    112    71    115    107    158    177    1130    1139    365

## Annotations
module load anaconda
python
import os
import pandas as pd
os.chdir("/nv/vol326/cphg_caseylab/FirstExposures/casey_grc_rnaseq_2/")
ann = pd.read_excel("lookup_casey_grc_rnaseq_2.xlsx")
os.chdir("/home/chd5n/projects/lentivirus/")
ann_dict = pd.read_csv("dict.tsv")

## Code
```
# bash
cd /nv/vol326/cphg_caseylab/lentivirus
awk -v FS="," -v OFS="," 'NR==1 {$1="gene_id"; print $0} $1=="ENSG00000196167" {print $0} $1=="ENSG00000214290" {print $0}' Lenti_counts.csv > lenti_colca_counts.csv

# R

## prep data
x <- read.table("/nv/vol326/cphg_caseylab/lentivirus/lenti_colca_counts.csv", sep=",", header=TRUE)
y <- data.frame(t(x))
colnames(y) <- y$gene_id
colnames(y) <- c("COLCA1","COLCA2")
z <- y[2:nrow(y),]
z$COLCA1 <- as.numeric(z$COLCA1)
z$COLCA2 <- as.numeric(z$COLCA2)
z$line <- as.factor(substr(rownames(z),2,4))
z$condition <- as.factor(substr(rownames(z),7,9))

## base plots
plot(COLCA1~line, data=z)
plot(COLCA1~condition:line, data=z)

x <- factor(z[z[,"condition"]!="CC2",c("line")])
y <- z[z[,"condition"]!="CC2",c("COLCA1")]
boxplot(y~x, outline=FALSE, main=paste0("Expression of COLCA1 across cell lines"), ylab="Raw counts")
stripchart(y~x, col=c("blue","red","purple"), vertical=TRUE, pch=21, method="jitter", add=TRUE)

x <- factor(z[z[,"condition"]!="CC2",c("condition")])
y <- z[z[,"condition"]!="CC2",c("COLCA1")]
boxplot(y~x, outline=FALSE, main=paste0("Expression of COLCA1 across conditions"), ylab="Raw counts")
stripchart(y~x, col=c("blue","red","purple"), vertical=TRUE, pch=21, method="jitter", add=TRUE)

## ggplot2
library(ggplot2)

setwd("~/projects/lentivirus/")
file <- "target-expr-plots.pdf"
pdf(file=file)

a <- z[z[,"condition"]!="CC2",c("COLCA1","line","condition")]
ggplot(a, aes(x=condition, y=COLCA1)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), cex=5, aes(col=line)) +
  labs(title="Expression of COLCA1 across conditions", x="condition", y="raw counts") +
  theme_classic() +
  theme(axis.text=element_text(size=12), axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), plot.title=element_text(hjust=0.5, size=16, face="bold"))

a <- z[z[,"condition"]!="CC1",c("COLCA2","line","condition")]
ggplot(a, aes(x=condition, y=COLCA2)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), cex=5, aes(col=line)) +
  labs(title="Expression of COLCA2 across conditions", x="condition", y="raw counts") +
  theme_classic() +
  theme(axis.text=element_text(size=12), axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), plot.title=element_text(hjust=0.5, size=16, face="bold"))

dev.off()

# python
import pandas as pd
import os
x <- pd.read_csv
```
