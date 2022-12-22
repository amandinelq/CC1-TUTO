CC1-TUTO
================
2022-12-22

\#charger toutes les librairies

``` r
library(dada2)
```

    ## Loading required package: Rcpp

``` r
library(phyloseq)
library(DECIPHER)
```

    ## Loading required package: Biostrings

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
    ##     table, tapply, union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:phyloseq':
    ## 
    ##     distance

    ## Loading required package: XVector

    ## Loading required package: GenomeInfoDb

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

    ## Loading required package: RSQLite

    ## Loading required package: parallel

``` r
library(phangorn)
```

    ## Loading required package: ape

    ## 
    ## Attaching package: 'ape'

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     complement

``` r
library(ggplot2)
library(gridExtra)
```

    ## 
    ## Attaching package: 'gridExtra'

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     combine

\#créer variable miseq_path = variable pour aller dans le chemin de
miseq, et lister les ficher qu’il y a dans miseq_path

``` r
miseq_path<-"/home/rstudio/MiSeq_SOP"
list.files(miseq_path)
```

    ##  [1] "F3D0_S188_L001_R1_001.fastq"   "F3D0_S188_L001_R2_001.fastq"  
    ##  [3] "F3D1_S189_L001_R1_001.fastq"   "F3D1_S189_L001_R2_001.fastq"  
    ##  [5] "F3D141_S207_L001_R1_001.fastq" "F3D141_S207_L001_R2_001.fastq"
    ##  [7] "F3D142_S208_L001_R1_001.fastq" "F3D142_S208_L001_R2_001.fastq"
    ##  [9] "F3D143_S209_L001_R1_001.fastq" "F3D143_S209_L001_R2_001.fastq"
    ## [11] "F3D144_S210_L001_R1_001.fastq" "F3D144_S210_L001_R2_001.fastq"
    ## [13] "F3D145_S211_L001_R1_001.fastq" "F3D145_S211_L001_R2_001.fastq"
    ## [15] "F3D146_S212_L001_R1_001.fastq" "F3D146_S212_L001_R2_001.fastq"
    ## [17] "F3D147_S213_L001_R1_001.fastq" "F3D147_S213_L001_R2_001.fastq"
    ## [19] "F3D148_S214_L001_R1_001.fastq" "F3D148_S214_L001_R2_001.fastq"
    ## [21] "F3D149_S215_L001_R1_001.fastq" "F3D149_S215_L001_R2_001.fastq"
    ## [23] "F3D150_S216_L001_R1_001.fastq" "F3D150_S216_L001_R2_001.fastq"
    ## [25] "F3D2_S190_L001_R1_001.fastq"   "F3D2_S190_L001_R2_001.fastq"  
    ## [27] "F3D3_S191_L001_R1_001.fastq"   "F3D3_S191_L001_R2_001.fastq"  
    ## [29] "F3D5_S193_L001_R1_001.fastq"   "F3D5_S193_L001_R2_001.fastq"  
    ## [31] "F3D6_S194_L001_R1_001.fastq"   "F3D6_S194_L001_R2_001.fastq"  
    ## [33] "F3D7_S195_L001_R1_001.fastq"   "F3D7_S195_L001_R2_001.fastq"  
    ## [35] "F3D8_S196_L001_R1_001.fastq"   "F3D8_S196_L001_R2_001.fastq"  
    ## [37] "F3D9_S197_L001_R1_001.fastq"   "F3D9_S197_L001_R2_001.fastq"  
    ## [39] "filtered"                      "HMP_MOCK.v35.fasta"           
    ## [41] "Mock_S280_L001_R1_001.fastq"   "Mock_S280_L001_R2_001.fastq"  
    ## [43] "mouse.dpw.metadata"            "mouse.time.design"            
    ## [45] "stability.batch"               "stability.files"

\#créer nouv variables qui reçoivent tous les noms de fichiers qui se
terminent par -R1 ou -R2 et les tries par ordre alphabetique

``` r
fnFs <- sort(list.files(miseq_path, pattern="_R1_001.fastq"))
fnRs <- sort(list.files(miseq_path, pattern="_R2_001.fastq"))
#extraire un echantillon
sampleNames <- sapply(strsplit(fnFs, "_"), `[`, 1)
#specifier la voie complete vers fnFs et FnRs
fnFs <- file.path(miseq_path, fnFs)
fnRs <- file.path(miseq_path, fnRs)
#afficher les 3 premieres sequences de fnFs
fnFs[1:3]
```

    ## [1] "/home/rstudio/MiSeq_SOP/F3D0_S188_L001_R1_001.fastq"  
    ## [2] "/home/rstudio/MiSeq_SOP/F3D1_S189_L001_R1_001.fastq"  
    ## [3] "/home/rstudio/MiSeq_SOP/F3D141_S207_L001_R1_001.fastq"

``` r
#afficher les 3 premieres sequences de fnRs
fnRs[1:3]
```

    ## [1] "/home/rstudio/MiSeq_SOP/F3D0_S188_L001_R2_001.fastq"  
    ## [2] "/home/rstudio/MiSeq_SOP/F3D1_S189_L001_R2_001.fastq"  
    ## [3] "/home/rstudio/MiSeq_SOP/F3D141_S207_L001_R2_001.fastq"

\#profils qualité des lectures : obtient graph

``` r
plotQualityProfile(fnFs[1:2])
```

    ## Warning: The `<scale>` argument of `guides()` cannot be `FALSE`. Use "none" instead as
    ## of ggplot2 3.3.4.
    ## ℹ The deprecated feature was likely used in the dada2 package.
    ##   Please report the issue at <]8;;https://github.com/benjjneb/dada2/issueshttps://github.com/benjjneb/dada2/issues]8;;>.

![](CC1-TUTO_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
filt_path <- file.path(miseq_path, "filtered")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))
```

\#etape de filtration de qualité

``` r
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE)
head(out)
```

    ##                               reads.in reads.out
    ## F3D0_S188_L001_R1_001.fastq       7793      7113
    ## F3D1_S189_L001_R1_001.fastq       5869      5299
    ## F3D141_S207_L001_R1_001.fastq     5958      5463
    ## F3D142_S208_L001_R1_001.fastq     3183      2914
    ## F3D143_S209_L001_R1_001.fastq     3178      2941
    ## F3D144_S210_L001_R1_001.fastq     4827      4312

``` r
#truncLen veut dire qu'on coupe à 160 (déterminé jusqu'à où le score de qualité est acceptable )sur le R1=forward et 240 sur le R2=reverse, faut bien regarder la longueur des fragments pour garder une superposition des deux lors de l'alignement (overlap) si on coupe trop court on en aura pas 
#maxN=0 quand séquenceur sait pas quelle pb c'est il met un N, donc on dit que si il y a au moins 1 N dans la seq on l'enlève car sera de mauvaise qualité 
#truncQ : a chaque fois que le long d'une sequence on voit apparaitre un score de qualié qui est inférieur à Q20 il coupe la séquence à ce niveau
#Trimleft : enlever les amorces à gauches (18 premiers nucléotides)
#filter and trim : fonction qui permet de faire la filtration quelité des séquences 
#obtient read.in : nbr de séquences qu'il avait avant et read.out : nbr de séquences qu'il obtient après les avoir filtré 
```

\#dada2: calcul des erreurs de séquençage \#première étape :
déréplication –\> pour qu’il reste que les séquences uniques avec au
bout le nombre de fois où elles apparaissaient

``` r
derepFs <- derepFastq(filtFs, verbose=TRUE)
```

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D0_F_filt.fastq.gz

    ## Encountered 1979 unique sequences from 7113 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D1_F_filt.fastq.gz

    ## Encountered 1639 unique sequences from 5299 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D141_F_filt.fastq.gz

    ## Encountered 1477 unique sequences from 5463 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D142_F_filt.fastq.gz

    ## Encountered 904 unique sequences from 2914 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D143_F_filt.fastq.gz

    ## Encountered 939 unique sequences from 2941 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D144_F_filt.fastq.gz

    ## Encountered 1267 unique sequences from 4312 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D145_F_filt.fastq.gz

    ## Encountered 1756 unique sequences from 6741 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D146_F_filt.fastq.gz

    ## Encountered 1438 unique sequences from 4560 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D147_F_filt.fastq.gz

    ## Encountered 3590 unique sequences from 15637 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D148_F_filt.fastq.gz

    ## Encountered 2762 unique sequences from 11413 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D149_F_filt.fastq.gz

    ## Encountered 3021 unique sequences from 12017 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D150_F_filt.fastq.gz

    ## Encountered 1566 unique sequences from 5032 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D2_F_filt.fastq.gz

    ## Encountered 3707 unique sequences from 18075 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D3_F_filt.fastq.gz

    ## Encountered 1479 unique sequences from 6250 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D5_F_filt.fastq.gz

    ## Encountered 1195 unique sequences from 4052 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D6_F_filt.fastq.gz

    ## Encountered 1832 unique sequences from 7369 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D7_F_filt.fastq.gz

    ## Encountered 1183 unique sequences from 4765 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D8_F_filt.fastq.gz

    ## Encountered 1382 unique sequences from 4871 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D9_F_filt.fastq.gz

    ## Encountered 1709 unique sequences from 6504 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/Mock_F_filt.fastq.gz

    ## Encountered 897 unique sequences from 4314 total sequences read.

``` r
derepRs <- derepFastq(filtRs, verbose=TRUE)
```

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D0_R_filt.fastq.gz

    ## Encountered 1660 unique sequences from 7113 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D1_R_filt.fastq.gz

    ## Encountered 1349 unique sequences from 5299 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D141_R_filt.fastq.gz

    ## Encountered 1335 unique sequences from 5463 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D142_R_filt.fastq.gz

    ## Encountered 853 unique sequences from 2914 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D143_R_filt.fastq.gz

    ## Encountered 880 unique sequences from 2941 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D144_R_filt.fastq.gz

    ## Encountered 1286 unique sequences from 4312 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D145_R_filt.fastq.gz

    ## Encountered 1803 unique sequences from 6741 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D146_R_filt.fastq.gz

    ## Encountered 1265 unique sequences from 4560 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D147_R_filt.fastq.gz

    ## Encountered 3414 unique sequences from 15637 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D148_R_filt.fastq.gz

    ## Encountered 2522 unique sequences from 11413 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D149_R_filt.fastq.gz

    ## Encountered 2771 unique sequences from 12017 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D150_R_filt.fastq.gz

    ## Encountered 1415 unique sequences from 5032 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D2_R_filt.fastq.gz

    ## Encountered 3290 unique sequences from 18075 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D3_R_filt.fastq.gz

    ## Encountered 1390 unique sequences from 6250 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D5_R_filt.fastq.gz

    ## Encountered 1134 unique sequences from 4052 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D6_R_filt.fastq.gz

    ## Encountered 1635 unique sequences from 7369 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D7_R_filt.fastq.gz

    ## Encountered 1084 unique sequences from 4765 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D8_R_filt.fastq.gz

    ## Encountered 1161 unique sequences from 4871 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D9_R_filt.fastq.gz

    ## Encountered 1502 unique sequences from 6504 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/Mock_R_filt.fastq.gz

    ## Encountered 732 unique sequences from 4314 total sequences read.

``` r
# Name the derep-class objects by the sample names
names(derepFs) <- sampleNames
names(derepRs) <- sampleNames
```

\#model d’erreur

``` r
errF <- learnErrors(filtFs, multithread=TRUE)
```

    ## 33514080 total bases in 139642 reads from 20 samples will be used for learning the error rates.

``` r
errR <- learnErrors(filtRs, multithread=TRUE)
```

    ## 22342720 total bases in 139642 reads from 20 samples will be used for learning the error rates.

``` r
plotErrors(errF)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](CC1-TUTO_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
plotErrors(errR)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](CC1-TUTO_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->

\#corriger les erreurs

``` r
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
```

    ## Sample 1 - 7113 reads in 1979 unique sequences.
    ## Sample 2 - 5299 reads in 1639 unique sequences.
    ## Sample 3 - 5463 reads in 1477 unique sequences.
    ## Sample 4 - 2914 reads in 904 unique sequences.
    ## Sample 5 - 2941 reads in 939 unique sequences.
    ## Sample 6 - 4312 reads in 1267 unique sequences.
    ## Sample 7 - 6741 reads in 1756 unique sequences.
    ## Sample 8 - 4560 reads in 1438 unique sequences.
    ## Sample 9 - 15637 reads in 3590 unique sequences.
    ## Sample 10 - 11413 reads in 2762 unique sequences.
    ## Sample 11 - 12017 reads in 3021 unique sequences.
    ## Sample 12 - 5032 reads in 1566 unique sequences.
    ## Sample 13 - 18075 reads in 3707 unique sequences.
    ## Sample 14 - 6250 reads in 1479 unique sequences.
    ## Sample 15 - 4052 reads in 1195 unique sequences.
    ## Sample 16 - 7369 reads in 1832 unique sequences.
    ## Sample 17 - 4765 reads in 1183 unique sequences.
    ## Sample 18 - 4871 reads in 1382 unique sequences.
    ## Sample 19 - 6504 reads in 1709 unique sequences.
    ## Sample 20 - 4314 reads in 897 unique sequences.

``` r
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
```

    ## Sample 1 - 7113 reads in 1660 unique sequences.
    ## Sample 2 - 5299 reads in 1349 unique sequences.
    ## Sample 3 - 5463 reads in 1335 unique sequences.
    ## Sample 4 - 2914 reads in 853 unique sequences.
    ## Sample 5 - 2941 reads in 880 unique sequences.
    ## Sample 6 - 4312 reads in 1286 unique sequences.
    ## Sample 7 - 6741 reads in 1803 unique sequences.
    ## Sample 8 - 4560 reads in 1265 unique sequences.
    ## Sample 9 - 15637 reads in 3414 unique sequences.
    ## Sample 10 - 11413 reads in 2522 unique sequences.
    ## Sample 11 - 12017 reads in 2771 unique sequences.
    ## Sample 12 - 5032 reads in 1415 unique sequences.
    ## Sample 13 - 18075 reads in 3290 unique sequences.
    ## Sample 14 - 6250 reads in 1390 unique sequences.
    ## Sample 15 - 4052 reads in 1134 unique sequences.
    ## Sample 16 - 7369 reads in 1635 unique sequences.
    ## Sample 17 - 4765 reads in 1084 unique sequences.
    ## Sample 18 - 4871 reads in 1161 unique sequences.
    ## Sample 19 - 6504 reads in 1502 unique sequences.
    ## Sample 20 - 4314 reads in 732 unique sequences.

``` r
#pour print les data
dadaFs[[1]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 128 sequence variants were inferred from 1979 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

\#faire l’allignement des R1 et R2

``` r
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)
```

\#créer table d’observation des séquences :

``` r
seqtabAll <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])
table(nchar(getSequences(seqtabAll)))
```

    ## 
    ## 251 252 253 254 255 
    ##   1  85 186   5   2

``` r
#on importe toutes les séquences de la table sauf celle Mock (car est une séquence artificielle introduite pour vérifier que ça marche)
#deuxième ligne = nombre de caractères 
#troisième ligne = nombre de séquences qui ont ce nombre de caractères 
#permet de vérifier que l'allignement est bien fait 
```

\#enlever les chimères = séquences avec un bout de séquence d’une
bactérie et un bout d’une autre bactérie , se produit pendant la PCR
lorsque l’ARNpol se décroche avant la fin

``` r
seqtabNoC <- removeBimeraDenovo(seqtabAll)
```

\#annotation taxonomique

\#1-télécharger fichier

``` bash
cd ~
wget https://zenodo.org/record/4310151/files/rdp_train_set_18.fa.gz
```

    ## --2022-12-22 13:23:12--  https://zenodo.org/record/4310151/files/rdp_train_set_18.fa.gz
    ## Resolving zenodo.org (zenodo.org)... 188.185.124.72
    ## Connecting to zenodo.org (zenodo.org)|188.185.124.72|:443... connected.
    ## HTTP request sent, awaiting response... 200 OK
    ## Length: 5760850 (5.5M) [application/octet-stream]
    ## Saving to: ‘rdp_train_set_18.fa.gz.9’
    ## 
    ##      0K .......... .......... .......... .......... ..........  0% 13.9M 0s
    ##     50K .......... .......... .......... .......... ..........  1% 15.3M 0s
    ##    100K .......... .......... .......... .......... ..........  2% 8.23M 0s
    ##    150K .......... .......... .......... .......... ..........  3% 85.3M 0s
    ##    200K .......... .......... .......... .......... ..........  4% 16.6M 0s
    ##    250K .......... .......... .......... .......... ..........  5% 13.5M 0s
    ##    300K .......... .......... .......... .......... ..........  6%  106M 0s
    ##    350K .......... .......... .......... .......... ..........  7% 13.9M 0s
    ##    400K .......... .......... .......... .......... ..........  7% 15.4M 0s
    ##    450K .......... .......... .......... .......... ..........  8% 86.2M 0s
    ##    500K .......... .......... .......... .......... ..........  9% 16.0M 0s
    ##    550K .......... .......... .......... .......... .......... 10% 61.4M 0s
    ##    600K .......... .......... .......... .......... .......... 11% 16.5M 0s
    ##    650K .......... .......... .......... .......... .......... 12% 77.1M 0s
    ##    700K .......... .......... .......... .......... .......... 13% 16.3M 0s
    ##    750K .......... .......... .......... .......... .......... 14% 15.4M 0s
    ##    800K .......... .......... .......... .......... .......... 15% 70.7M 0s
    ##    850K .......... .......... .......... .......... .......... 15% 15.3M 0s
    ##    900K .......... .......... .......... .......... .......... 16% 91.1M 0s
    ##    950K .......... .......... .......... .......... .......... 17% 15.9M 0s
    ##   1000K .......... .......... .......... .......... .......... 18% 80.5M 0s
    ##   1050K .......... .......... .......... .......... .......... 19% 15.4M 0s
    ##   1100K .......... .......... .......... .......... .......... 20%  104M 0s
    ##   1150K .......... .......... .......... .......... .......... 21% 15.4M 0s
    ##   1200K .......... .......... .......... .......... .......... 22% 73.4M 0s
    ##   1250K .......... .......... .......... .......... .......... 23% 16.8M 0s
    ##   1300K .......... .......... .......... .......... .......... 23% 59.5M 0s
    ##   1350K .......... .......... .......... .......... .......... 24% 15.2M 0s
    ##   1400K .......... .......... .......... .......... .......... 25% 16.5M 0s
    ##   1450K .......... .......... .......... .......... .......... 26% 77.8M 0s
    ##   1500K .......... .......... .......... .......... .......... 27% 15.5M 0s
    ##   1550K .......... .......... .......... .......... .......... 28% 15.2M 0s
    ##   1600K .......... .......... .......... .......... .......... 29% 68.8M 0s
    ##   1650K .......... .......... .......... .......... .......... 30% 16.7M 0s
    ##   1700K .......... .......... .......... .......... .......... 31% 49.2M 0s
    ##   1750K .......... .......... .......... .......... .......... 31% 17.1M 0s
    ##   1800K .......... .......... .......... .......... .......... 32% 71.9M 0s
    ##   1850K .......... .......... .......... .......... .......... 33% 16.4M 0s
    ##   1900K .......... .......... .......... .......... .......... 34% 15.1M 0s
    ##   1950K .......... .......... .......... .......... .......... 35% 55.9M 0s
    ##   2000K .......... .......... .......... .......... .......... 36% 16.2M 0s
    ##   2050K .......... .......... .......... .......... .......... 37% 61.9M 0s
    ##   2100K .......... .......... .......... .......... .......... 38% 16.0M 0s
    ##   2150K .......... .......... .......... .......... .......... 39% 16.9M 0s
    ##   2200K .......... .......... .......... .......... .......... 39% 46.2M 0s
    ##   2250K .......... .......... .......... .......... .......... 40% 16.5M 0s
    ##   2300K .......... .......... .......... .......... .......... 41% 69.8M 0s
    ##   2350K .......... .......... .......... .......... .......... 42% 16.9M 0s
    ##   2400K .......... .......... .......... .......... .......... 43% 49.2M 0s
    ##   2450K .......... .......... .......... .......... .......... 44% 17.0M 0s
    ##   2500K .......... .......... .......... .......... .......... 45% 14.3M 0s
    ##   2550K .......... .......... .......... .......... .......... 46% 74.2M 0s
    ##   2600K .......... .......... .......... .......... .......... 47% 16.9M 0s
    ##   2650K .......... .......... .......... .......... .......... 47% 85.3M 0s
    ##   2700K .......... .......... .......... .......... .......... 48% 15.2M 0s
    ##   2750K .......... .......... .......... .......... .......... 49% 69.3M 0s
    ##   2800K .......... .......... .......... .......... .......... 50% 16.2M 0s
    ##   2850K .......... .......... .......... .......... .......... 51% 78.9M 0s
    ##   2900K .......... .......... .......... .......... .......... 52% 15.9M 0s
    ##   2950K .......... .......... .......... .......... .......... 53% 57.8M 0s
    ##   3000K .......... .......... .......... .......... .......... 54% 16.6M 0s
    ##   3050K .......... .......... .......... .......... .......... 55% 76.2M 0s
    ##   3100K .......... .......... .......... .......... .......... 55% 15.8M 0s
    ##   3150K .......... .......... .......... .......... .......... 56% 15.2M 0s
    ##   3200K .......... .......... .......... .......... .......... 57% 83.3M 0s
    ##   3250K .......... .......... .......... .......... .......... 58% 16.1M 0s
    ##   3300K .......... .......... .......... .......... .......... 59% 41.5M 0s
    ##   3350K .......... .......... .......... .......... .......... 60% 18.1M 0s
    ##   3400K .......... .......... .......... .......... .......... 61% 81.6M 0s
    ##   3450K .......... .......... .......... .......... .......... 62% 15.4M 0s
    ##   3500K .......... .......... .......... .......... .......... 63% 49.7M 0s
    ##   3550K .......... .......... .......... .......... .......... 63% 13.5M 0s
    ##   3600K .......... .......... .......... .......... .......... 64% 86.6M 0s
    ##   3650K .......... .......... .......... .......... .......... 65% 12.9M 0s
    ##   3700K .......... .......... .......... .......... .......... 66% 96.2M 0s
    ##   3750K .......... .......... .......... .......... .......... 67% 13.2M 0s
    ##   3800K .......... .......... .......... .......... .......... 68%  102M 0s
    ##   3850K .......... .......... .......... .......... .......... 69% 12.9M 0s
    ##   3900K .......... .......... .......... .......... .......... 70% 90.3M 0s
    ##   3950K .......... .......... .......... .......... .......... 71% 13.7M 0s
    ##   4000K .......... .......... .......... .......... .......... 71% 15.2M 0s
    ##   4050K .......... .......... .......... .......... .......... 72% 16.6M 0s
    ##   4100K .......... .......... .......... .......... .......... 73% 52.3M 0s
    ##   4150K .......... .......... .......... .......... .......... 74% 15.3M 0s
    ##   4200K .......... .......... .......... .......... .......... 75% 54.5M 0s
    ##   4250K .......... .......... .......... .......... .......... 76% 16.5M 0s
    ##   4300K .......... .......... .......... .......... .......... 77% 94.7M 0s
    ##   4350K .......... .......... .......... .......... .......... 78% 14.4M 0s
    ##   4400K .......... .......... .......... .......... .......... 79% 17.6M 0s
    ##   4450K .......... .......... .......... .......... .......... 79% 72.1M 0s
    ##   4500K .......... .......... .......... .......... .......... 80% 16.6M 0s
    ##   4550K .......... .......... .......... .......... .......... 81% 44.2M 0s
    ##   4600K .......... .......... .......... .......... .......... 82% 17.2M 0s
    ##   4650K .......... .......... .......... .......... .......... 83% 82.1M 0s
    ##   4700K .......... .......... .......... .......... .......... 84% 16.0M 0s
    ##   4750K .......... .......... .......... .......... .......... 85% 13.9M 0s
    ##   4800K .......... .......... .......... .......... .......... 86% 66.9M 0s
    ##   4850K .......... .......... .......... .......... .......... 87% 16.9M 0s
    ##   4900K .......... .......... .......... .......... .......... 87% 49.9M 0s
    ##   4950K .......... .......... .......... .......... .......... 88% 16.3M 0s
    ##   5000K .......... .......... .......... .......... .......... 89% 75.9M 0s
    ##   5050K .......... .......... .......... .......... .......... 90% 16.5M 0s
    ##   5100K .......... .......... .......... .......... .......... 91% 70.4M 0s
    ##   5150K .......... .......... .......... .......... .......... 92% 16.0M 0s
    ##   5200K .......... .......... .......... .......... .......... 93% 57.6M 0s
    ##   5250K .......... .......... .......... .......... .......... 94% 17.2M 0s
    ##   5300K .......... .......... .......... .......... .......... 95% 63.9M 0s
    ##   5350K .......... .......... .......... .......... .......... 95% 14.7M 0s
    ##   5400K .......... .......... .......... .......... .......... 96% 79.0M 0s
    ##   5450K .......... .......... .......... .......... .......... 97% 16.7M 0s
    ##   5500K .......... .......... .......... .......... .......... 98% 63.8M 0s
    ##   5550K .......... .......... .......... .......... .......... 99% 14.9M 0s
    ##   5600K .......... .......... .....                           100%  105M=0.2s
    ## 
    ## 2022-12-22 13:23:13 (23.1 MB/s) - ‘rdp_train_set_18.fa.gz.9’ saved [5760850/5760850]

\#assigner taxo

``` r
fastaRef <- "/home/rstudio/rdp_train_set_18.fa.gz"
taxTab <- assignTaxonomy(seqtabNoC, refFasta = fastaRef, multithread=TRUE)
unname(head(taxTab))
```

    ##      [,1]       [,2]            [,3]          [,4]            [,5]            
    ## [1,] "Bacteria" "Bacteroidetes" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [2,] "Bacteria" "Bacteroidetes" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [3,] "Bacteria" "Bacteroidetes" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [4,] "Bacteria" "Bacteroidetes" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [5,] "Bacteria" "Bacteroidetes" "Bacteroidia" "Bacteroidales" "Bacteroidaceae"
    ## [6,] "Bacteria" "Bacteroidetes" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ##      [,6]         
    ## [1,] "Duncaniella"
    ## [2,] "Duncaniella"
    ## [3,] "Muribaculum"
    ## [4,] "Muribaculum"
    ## [5,] "Bacteroides"
    ## [6,] NA

\#construction arbre phylogenetique

``` r
seqs <- getSequences(seqtabNoC)
names(seqs) <- seqs # nomme les séquences de l'arbre 
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)#aligne les sequences 
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
        rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)
```

``` r
samdf <- read.csv("https://raw.githubusercontent.com/spholmes/F1000_workflow/master/data/MIMARKS_Data_combined.csv",header=TRUE)
#table à faire pour le projet, on la crée sur excel et on l'importe dans R avec cette fonction 
```

``` r
#gsub : première variable à rentrer est l'objet qu'il doit touver, la deuxième est ce par quoi il doit le remplacer --> ici trouve le pattern 00 et le remplace par rien 
samdf$SampleID <- paste0(gsub("00", "", samdf$host_subject_id), "D", samdf$age-21)
samdf <- samdf[!duplicated(samdf$SampleID),] #enlever les duplicats
rownames(seqtabAll) <- gsub("124", "125", rownames(seqtabAll)) #corriger l'écart
all(rownames(seqtabAll) %in% samdf$SampleID)
```

    ## [1] TRUE

``` r
rownames(samdf) <- samdf$SampleID
keep.cols <- c("collection_date", "biome", "target_gene", "target_subfragment",
"host_common_name", "host_subject_id", "age", "sex", "body_product", "tot_mass",
"diet", "family_relationship", "genotype", "SampleID") 
samdf <- samdf[rownames(seqtabAll), keep.cols]
```

``` r
ps <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxTab),phy_tree(fitGTR$tree))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # enlever les échantillons Mock
ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 218 taxa and 19 samples ]
    ## sample_data() Sample Data:       [ 19 samples by 14 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 218 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 218 tips and 216 internal nodes ]

``` r
#OTU: en ligne l'échantillon et en colonne sa séquence et on note le nombre de fois où elle apparait = sa fréuence 
#Sample data : échantillons en igne et en colonne les caractéristiques (mâle, femelle, age...)
#Phylogenique tree : arbe phylogenet
#taxTable : séquence en ligne et en colonne le phylium
```

``` r
#télécharge un ps.rds : rds= un objet R et dedans on à toutes les données OTU,Sample data...
#si on veut regarder la table OTU on fait : ps@OTU_table
ps_connect <-url("https://raw.githubusercontent.com/spholmes/F1000_workflow/master/data/ps.rds")
ps = readRDS(ps_connect)
ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 389 taxa and 360 samples ]
    ## sample_data() Sample Data:       [ 360 samples by 14 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 389 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 389 tips and 387 internal nodes ]

\#filtration taxonomique

``` r
#afficher les classements disponibles dans les données (kingdom, phylium, classe...)
rank_names(ps)
```

    ## [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"

``` r
#créer un tableau pour chaque phylium
table(tax_table(ps)[, "Phylum"], exclude = NULL)
```

    ## 
    ##              Actinobacteria               Bacteroidetes 
    ##                          13                          23 
    ## Candidatus_Saccharibacteria   Cyanobacteria/Chloroplast 
    ##                           1                           4 
    ##         Deinococcus-Thermus                  Firmicutes 
    ##                           1                         327 
    ##                Fusobacteria              Proteobacteria 
    ##                           1                          11 
    ##                 Tenericutes             Verrucomicrobia 
    ##                           1                           1 
    ##                        <NA> 
    ##                           6

``` r
#on enlève tt les séquences pour lequel pylum c'est NA ou rien ou uncharacterized
#les entités avec des annotations de phylum ambiguës sont supprimées
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
```

\#explorer la prévalence =nbr d’ech dans lequel apparaissent ces taxa

``` r
# Compile toute les données qui sont prévalentes = définit ici comme le nombre d'échantillons dans lesquels un taxon apparaît au moins une fois
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Ajouter une taxonomie et le nombre totale de read
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))
```

``` r
#Calculez les prévalences totales et moyennes des caractéristiques dans chaque phylum 
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
```

    ##                         Phylum         1     2
    ## 1               Actinobacteria 120.15385  1562
    ## 2                Bacteroidetes 265.52174  6107
    ## 3  Candidatus_Saccharibacteria 280.00000   280
    ## 4    Cyanobacteria/Chloroplast  64.25000   257
    ## 5          Deinococcus-Thermus  52.00000    52
    ## 6                   Firmicutes 179.24771 58614
    ## 7                 Fusobacteria   2.00000     2
    ## 8               Proteobacteria  59.09091   650
    ## 9                  Tenericutes 234.00000   234
    ## 10             Verrucomicrobia 104.00000   104

``` r
# Définir les phyliums à filtrer 
filterPhyla = c("Fusobacteria", "Deinococcus-Thermus")
# Filtrer les entrées avec Phylum non identifié
ps1 = subset_taxa(ps, !Phylum %in% filterPhyla)
ps1
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 381 taxa and 360 samples ]
    ## sample_data() Sample Data:       [ 360 samples by 14 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 381 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 381 tips and 379 internal nodes ]

\#Filtration par prévalence

``` r
# Sous-ensemble des embranchements restants
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +
  # Inclure une supposition pour le paramètre
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
```

![](CC1-TUTO_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

``` r
# Définir le seuil de prévalence à 5 % du total des échantillons
prevalenceThreshold = 0.05 * nsamples(ps)
prevalenceThreshold
```

    ## [1] 18

``` r
# Exécuter le filtre de prévalence à l'aide de la fonction `prune_taxa()`
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps)
```

\#taxons agglomérés : regroupe tt les séquences du même taxon \#utile
quand on considère que tt les bactéries d’un même genre ou espèce ont
une même fonction

``` r
# Combien de genres seraient présents après filtrage ?
length(get_taxa_unique(ps2, taxonomic.rank = "Genus"))
```

    ## [1] 49

``` r
ps3 = tax_glom(ps2, "Genus", NArm = TRUE)
```

``` r
h1 = 0.4
ps4 = tip_glom(ps2, h = h1)
```

``` r
#Ici, la fonction plot_tree() de phyloseq compare les données originales non filtrées, l'arbre après agglomération taxonomique et l'arbre après agglomération phylogénétique. 
multiPlotTitleTextSize = 15
p2tree = plot_tree(ps2, method = "treeonly",
                   ladderize = "left",
                   title = "Before Agglomeration") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p3tree = plot_tree(ps3, method = "treeonly",
                   ladderize = "left", title = "By Genus") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p4tree = plot_tree(ps4, method = "treeonly",
                   ladderize = "left", title = "By Height") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
#sont ensuite mis ensemble dans un graphique combiné à l'aide de gridExtra::grid.arrange.
grid.arrange(nrow = 1, p2tree, p3tree, p4tree)
```

![](CC1-TUTO_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

``` r
#premier arrbre = données brutes 
#deuxième arbre= on a regroupé tt les séquences qui ont le même genre donc on a plus qu'on arbre de genre 
#troisième arbre= fait regroupement par distance dans l'arbre : on reprends le premier arbre et tt les séquences qui sont à une distance évolutive plus petite que 0,4 entre elles on les agglomères ensemble (en une seule branche)
```

\#Transformation de la valeur d’abondance : normalisation

``` r
#transformation des valeurs d'abondnace via la fonction trasform_sample_count(). Le premeier argument de cette fonction est l'objet phyloseq que l'on veut transformer, le segond est la fonction R qui définit la transformation
#c'est l'étape de normalisation
plot_abundance = function(physeq,title = "",
                          Facet = "Order", Color = "Phylum"){
  # Arbitrary subset, based on Phylum, for plotting
  p1f = subset_taxa(physeq, Phylum %in% c("Firmicutes"))
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "sex",y = "Abundance",
                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}
```

``` r
#convertit les comptages de chaque échantillon en leurs fréquences=abondances relative
ps3ra = transform_sample_counts(ps3, function(x){x / sum(x)})
```

``` r
#Maintenant, nous traçons les valeurs d'abondance avant et après la transformation
plotBefore = plot_abundance(ps3,"")
```

    ## Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
    ## ℹ Please use tidy evaluation ideoms with `aes()`

``` r
plotAfter = plot_abundance(ps3ra,"")
# on les mets dans un graphique 
grid.arrange(nrow = 2,  plotBefore, plotAfter)
```

![](CC1-TUTO_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

``` r
#c'est un graphiuque en violon: donne ici une abondnace autour de 30
#une fois que la matrice à été normalisé on est en % (deuxième graph)
```

\#Sous-ensemble par taxonomie

``` r
#trace uniquement ce sous-ensemble taxonomique des données, Pour cela, nous sous-ensembleons avec la fonction, puis spécifions un rang taxonomique plus précis à l'argument de la fonction que nous avons défini ci-dessus
psOrd = subset_taxa(ps3ra, Order == "Lactobacillales")
plot_abundance(psOrd, Facet = "Genus", Color = NULL)
```

![](CC1-TUTO_files/figure-gfm/unnamed-chunk-35-1.png)<!-- -->

``` r
library(shiny)
library(miniUI)
#library(caret)
library(pls)
```

    ## 
    ## Attaching package: 'pls'

    ## The following object is masked from 'package:ape':
    ## 
    ##     mvr

    ## The following object is masked from 'package:stats':
    ## 
    ##     loadings

``` r
library(e1071)
library(ggplot2)
library(randomForest)
```

    ## randomForest 4.7-1.1

    ## Type rfNews() to see new features/changes/bug fixes.

    ## 
    ## Attaching package: 'randomForest'

    ## The following object is masked from 'package:gridExtra':
    ## 
    ##     combine

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     margin

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     combine

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:randomForest':
    ## 
    ##     combine

    ## The following object is masked from 'package:gridExtra':
    ## 
    ##     combine

    ## The following objects are masked from 'package:Biostrings':
    ## 
    ##     collapse, intersect, setdiff, setequal, union

    ## The following object is masked from 'package:GenomeInfoDb':
    ## 
    ##     intersect

    ## The following object is masked from 'package:XVector':
    ## 
    ##     slice

    ## The following objects are masked from 'package:IRanges':
    ## 
    ##     collapse, desc, intersect, setdiff, slice, union

    ## The following objects are masked from 'package:S4Vectors':
    ## 
    ##     first, intersect, rename, setdiff, setequal, union

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(ggrepel)
library(nlme)
```

    ## 
    ## Attaching package: 'nlme'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     collapse

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     collapse

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     collapse

``` r
library(devtools)
```

    ## Loading required package: usethis

``` r
library(reshape2)
library(PMA)
#library(structSSI)
library(ade4)
```

    ## 
    ## Attaching package: 'ade4'

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     score

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     score

``` r
library(ggnetwork)
library(intergraph)
library(scales)
#library(jfukuyama/phyloseqGraphTest)
library(genefilter)
library(impute)
```

``` r
#on insert des colones en plus à notre tableau 
qplot(sample_data(ps)$age, geom = "histogram",binwidth=20) + xlab("age")
```

    ## Warning: `qplot()` was deprecated in ggplot2 3.4.0.

![](CC1-TUTO_files/figure-gfm/unnamed-chunk-37-1.png)<!-- -->

``` r
qplot(log10(rowSums(otu_table(ps))),binwidth=0.2) +
  xlab("Logged counts-per-sample")
```

![](CC1-TUTO_files/figure-gfm/unnamed-chunk-37-2.png)<!-- -->

\#PCoA

``` r
sample_data(ps)$age_binned <- cut(sample_data(ps)$age,
                          breaks = c(0, 100, 200, 400))
levels(sample_data(ps)$age_binned) <- list(Young100="(0,100]", Mid100to200="(100,200]", Old200="(200,400]")
sample_data(ps)$family_relationship=gsub(" ","",sample_data(ps)$family_relationship)
pslog <- transform_sample_counts(ps, function(x) log(1 + x))
out.wuf.log <- ordinate(pslog, method = "MDS", distance = "wunifrac")
```

    ## Warning in UniFrac(physeq, weighted = TRUE, ...): Randomly assigning root as --
    ## GCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGTGCGTAGGCGGCTTTGCAAGTCAGAAGTGAAATCCATGGGCTTAACCCATGAACTGCTTTTGAAACTGCAGAGCTTGAGTGGAGTAGAGGTAGGCGGAATTCCCGGTGTAGCGGTGAAATGCGTAGAGATCGGGAGGAACACCAGTGGCGAAGGCGGCCTGCTGGGCTCTAACTGACGCTGAGGCACGAAAGCGTGGGTAG
    ## -- in the phylogenetic tree in the data you provided.

``` r
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "age_binned") +
  labs(col = "Binned Age") +
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](CC1-TUTO_files/figure-gfm/unnamed-chunk-38-1.png)<!-- -->

``` r
rel_abund <- t(apply(otu_table(ps), 1, function(x) x / sum(x)))
qplot(rel_abund[, 12], geom = "histogram",binwidth=0.05) +
  xlab("Relative abundance")
```

![](CC1-TUTO_files/figure-gfm/unnamed-chunk-39-1.png)<!-- -->

\#Différentes ordinations

``` r
outliers <- c("F5D165", "F6D165", "M3D175", "M4D175", "M5D175", "M6D175")
ps <- prune_samples(!(sample_names(ps) %in% outliers), ps)
#on enlève les échantillons avec moins de 1000 reads
which(!rowSums(otu_table(ps)) > 1000)
```

    ## F5D145 M1D149   M1D9 M2D125  M2D19 M3D148 M3D149   M3D3   M3D5   M3D8 
    ##     69    185    200    204    218    243    244    252    256    260

``` r
ps <- prune_samples(rowSums(otu_table(ps)) > 1000, ps)
pslog <- transform_sample_counts(ps, function(x) log(1 + x))
```

``` r
#faire la PCoA
out.pcoa.log <- ordinate(pslog,  method = "MDS", distance = "bray")
evals <- out.pcoa.log$values[,1]
plot_ordination(pslog, out.pcoa.log, color = "age_binned",
                  shape = "family_relationship") +
  labs(col = "Binned Age", shape = "Litter")+
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](CC1-TUTO_files/figure-gfm/unnamed-chunk-41-1.png)<!-- -->

``` r
#DPCoA : double principal coordinates analysis
out.dpcoa.log <- ordinate(pslog, method = "DPCoA")
evals <- out.dpcoa.log$eig
plot_ordination(pslog, out.dpcoa.log, color = "age_binned", label= "SampleID",
                  shape = "family_relationship") +
  labs(col = "Binned Age", shape = "Litter")+
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](CC1-TUTO_files/figure-gfm/unnamed-chunk-42-1.png)<!-- -->

``` r
plot_ordination(pslog, out.dpcoa.log, type = "species", color = "Phylum") +
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](CC1-TUTO_files/figure-gfm/unnamed-chunk-43-1.png)<!-- -->

``` r
out.wuf.log <- ordinate(pslog, method = "PCoA", distance ="wunifrac")
```

    ## Warning in UniFrac(physeq, weighted = TRUE, ...): Randomly assigning root as --
    ## GCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTCAGCAAGTTGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCCAAAACTACTGAGCTAGAGTACGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGATACTGACACTGAGGTGCGAAAGCGTGGGGAG
    ## -- in the phylogenetic tree in the data you provided.

``` r
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "age_binned",
                  shape = "family_relationship") +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  labs(col = "Binned Age", shape = "Litter")
```

![](CC1-TUTO_files/figure-gfm/unnamed-chunk-44-1.png)<!-- -->

``` r
#Why are the ordination plots so far from square?
abund <- otu_table(pslog)
abund_ranks <- t(apply(abund, 1, rank))
abund_ranks <- abund_ranks - 329
abund_ranks[abund_ranks < 1] <- 1
```

``` r
library(dplyr)
library(reshape2)
abund_df <- melt(abund, value.name = "abund") %>%
  left_join(melt(abund_ranks, value.name = "rank"))
```

    ## Joining, by = c("Var1", "Var2")

``` r
colnames(abund_df) <- c("sample", "seq", "abund", "rank")
abund_df <- melt(abund, value.name = "abund") %>%
  left_join(melt(abund_ranks, value.name = "rank"))
```

    ## Joining, by = c("Var1", "Var2")

``` r
colnames(abund_df) <- c("sample", "seq", "abund", "rank")
sample_ix <- sample(1:nrow(abund_df), 8)
ggplot(abund_df %>%
         filter(sample %in% abund_df$sample[sample_ix])) +
  geom_point(aes(x = abund, y = rank, col = sample),
             position = position_jitter(width = 0.2), size = 1.5) +
  labs(x = "Abundance", y = "Thresholded rank") +
  scale_color_brewer(palette = "Set2")
```

![](CC1-TUTO_files/figure-gfm/unnamed-chunk-46-1.png)<!-- -->

``` r
library(ade4)
ranks_pca <- dudi.pca(abund_ranks, scannf = F, nf = 3)
row_scores <- data.frame(li = ranks_pca$li,
                         SampleID = rownames(abund_ranks))
col_scores <- data.frame(co = ranks_pca$co,
                         seq = colnames(abund_ranks))
tax <- tax_table(ps) %>%
  data.frame(stringsAsFactors = FALSE)
tax$seq <- rownames(tax)
main_orders <- c("Clostridiales", "Bacteroidales", "Lactobacillales",
                 "Coriobacteriales")
tax$Order[!(tax$Order %in% main_orders)] <- "Other"
tax$Order <- factor(tax$Order, levels = c(main_orders, "Other"))
tax$otu_id <- seq_len(ncol(otu_table(ps)))
row_scores <- row_scores %>%
  left_join(sample_data(pslog))
```

    ## Joining, by = "SampleID"

``` r
col_scores <- col_scores %>%
  left_join(tax)
```

    ## Joining, by = "seq"

``` r
evals_prop <- 100 * (ranks_pca$eig / sum(ranks_pca$eig))
ggplot() +
  geom_point(data = row_scores, aes(x = li.Axis1, y = li.Axis2), shape = 2) +
  geom_point(data = col_scores, aes(x = 25 * co.Comp1, y = 25 * co.Comp2, col = Order),
             size = .3, alpha = 0.6) +
  scale_color_brewer(palette = "Set2") +
  facet_grid(~ age_binned) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
       y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
  coord_fixed(sqrt(ranks_pca$eig[2] / ranks_pca$eig[1])) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))
```

![](CC1-TUTO_files/figure-gfm/unnamed-chunk-48-1.png)<!-- -->

``` r
#Canonical correspondence
ps_ccpna <- ordinate(pslog, "CCA", formula = pslog ~ age_binned + family_relationship)
```

``` r
library(ggrepel)
ps_scores <- vegan::scores(ps_ccpna)
sites <- data.frame(ps_scores$sites)
sites$SampleID <- rownames(sites)
sites <- sites %>%
  left_join(sample_data(ps))
```

    ## Joining, by = "SampleID"

``` r
species <- data.frame(ps_scores$species)
species$otu_id <- seq_along(colnames(otu_table(ps)))
species <- species %>%
  left_join(tax)
```

    ## Joining, by = "otu_id"

``` r
evals_prop <- 100 * ps_ccpna$CCA$eig[1:2] / sum(ps_ccpna$CA$eig)
ggplot() +
  geom_point(data = sites, aes(x = CCA1, y = CCA2), shape = 2, alpha = 0.5) +
  geom_point(data = species, aes(x = CCA1, y = CCA2, col = Order), size = 0.5) +
  geom_text_repel(data = species %>% filter(CCA2 < -2),
                    aes(x = CCA1, y = CCA2, label = otu_id),
            size = 1.5, segment.size = 0.1) +
  facet_grid(. ~ family_relationship) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
        y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
  scale_color_brewer(palette = "Set2") +
  coord_fixed(sqrt(ps_ccpna$CCA$eig[2] / ps_ccpna$CCA$eig[1])*0.45   ) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))
```

![](CC1-TUTO_files/figure-gfm/unnamed-chunk-50-1.png)<!-- -->
