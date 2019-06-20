# Results for organoid lentivirus pilot

## Preliminary gene expression analysis

### Parental line (Par) vs scramble (SCR)

__PCA:__
![parental vs scramble](images/parScrPca.png)

__Results summary from DESeq2:__
```
> summary(res)

out of 18236 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 7, 0.038%
LFC < 0 (down)     : 3, 0.016%
outliers [1]       : 0, 0%
low counts [2]     : 3536, 19%
```

__LFC vs mean norm counts (MA plot):__
![parental vs scramble](images/parScrMa.png)

### CC1 treatment vs SCR

__PCA:__
![cc1 vs scramble](images/cc1ScrPca.png)

__Results summary from DESeq2:__
```
> summary(res)

out of 17791 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 23, 0.13%
LFC < 0 (down)     : 88, 0.49%
outliers [1]       : 0, 0%
low counts [2]     : 3794, 21%
```

__LFC vs mean norm counts (MA plot):__
![cc1 vs scramble](images/cc1ScrMa.png)

### CC2 treatment vs SCR

__PCA:__
![cc2 vs scramble](images/cc2ScrPca.png)

__Results summary from DESeq2:__
```
out of 18064 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 1150, 6.4%
LFC < 0 (down)     : 1811, 10%
outliers [1]       : 0, 0%
low counts [2]     : 1401, 7.8%
```

__LFC vs mean norm counts (MA plot):__
![cc2 vs scramble](images/cc2ScrMa.png)
