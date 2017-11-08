Querying transcript length using Bioconductor

Overview
========

The goal of this analysis is to explore how Bioconductor handles transcript length. This is motivated by the realization that the functions which generate transcript database packages (e.g. `makeTxDbFromGFF`) result in different numbers for "transcript length" for intron-containing genes, compared with what is listed in the annotation TriTrypDB text files for the same gene.

For example, currently, the transcript length that gets computed uses the gene/transcript boundaries for [LmjF.02.0100](http://www.genedb.org/gene/LmjF.02.0100?actionName=%2FQuery%2FquickSearch&resultsSize=1&taxonNodeName=Root) which has the gene model:

    TABLE: Gene Model
    [Name]          [Type]  [Start] [End]   [Is Reversed]
    LmjF.02.0100    exon    35740   36189   1
    LmjF.02.0100    intron  35739   35739   1
    LmjF.02.0100    exon    34455   35738   1
    LmjF.02.0100    intron  34447   34454   1
    LmjF.02.0100    exon    33424   34446   1

The correct transcript length should be:

    > (36189-35740+1) + (35738-34455+1) + (34446-33424+1)
    [1] 2757

Which is what is listed in the annotation txt file (TriTrypDB-9.0\_LmajorFriedlinGene.txt). However, the TxDb generation logic does not account for introns and instead arrives at:

    > 36189 - 33424 + 1
    [1] 2766

This is a result of using the `useGenesAsTranscripts=TRUE` parameter in `makeTranscriptDbFromGFF()` (Bioconductor 3.0) or `gffTxName="gene"` parameter in `makeTxDbFromGFF()` (Bioconductor 3.1) to avoid excluding ncRNAs which do not have an `mRNA` row in the source GFF files.

**TriTrypDB-9.0\_LmajorFriedlin\_genes.gff**

    LmjF.02 TriTrypDB   gene    33424   36189   .   -   .   ID=LmjF.02.0100;Name=LmjF.02.0100;description=hypothetical+protein%2C+conserved+%28pseudogene%29;size=2766;web_id=LmjF.02.0100;locus_tag=LmjF.02.0100;size=2766;Alias=LmjF2.0100,LmjF02.0100,LmjF.02.0100,LmjF02.0100:pseudogenic_transcript,LmjF.02.0100:pseudogenic.transcript,LmjF02.0100:pseudogenic_transcript:pep,LmjF.02.0100:pseudogenic.transcript:pep
    LmjF.02 TriTrypDB   mRNA    33424   36189   .   -   .   ID=rna_LmjF.02.0100-1;Name=LmjF.02.0100-1;description=LmjF.02.0100-1;size=2766;Parent=LmjF.02.0100;Ontology_term=GO:0003676,GO:0008270;Dbxref=ApiDB:LmjF.02.0100,taxon:347515
    LmjF.02 TriTrypDB   CDS 33424   34446   .   -   0   ID=cds_LmjF.02.0100-3;Name=cds;description=.;size=1023;Parent=rna_LmjF.02.0100-1
    LmjF.02 TriTrypDB   CDS 34455   35738   .   -   0   ID=cds_LmjF.02.0100-2;Name=cds;description=.;size=1284;Parent=rna_LmjF.02.0100-1
    LmjF.02 TriTrypDB   CDS 35740   36189   .   -   0   ID=cds_LmjF.02.0100-1;Name=cds;description=.;size=450;Parent=rna_LmjF.02.0100-1
    LmjF.02 TriTrypDB   exon    35740   36189   .   -   .   ID=exon_LmjF.02.0100-1;Name=exon;description=exon;size=450;Parent=rna_LmjF.02.0100-1
    LmjF.02 TriTrypDB   exon    34455   35738   .   -   .   ID=exon_LmjF.02.0100-2;Name=exon;description=exon;size=1284;Parent=rna_LmjF.02.0100-1
    LmjF.02 TriTrypDB   exon    33424   34446   .   -   .   ID=exon_LmjF.02.0100-3;Name=exon;description=exon;size=1023;Parent=rna_LmjF.02.0100-1

*L. major* genes that are likely affected (i.e. contain introns):

``` r
library("rtracklayer")
```

    ## Loading required package: GenomicRanges

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, cbind, colMeans,
    ##     colnames, colSums, do.call, duplicated, eval, evalq, Filter,
    ##     Find, get, grep, grepl, intersect, is.unsorted, lapply,
    ##     lengths, Map, mapply, match, mget, order, paste, pmax,
    ##     pmax.int, pmin, pmin.int, Position, rank, rbind, Reduce,
    ##     rowMeans, rownames, rowSums, sapply, setdiff, sort, table,
    ##     tapply, union, unique, unsplit, which, which.max, which.min

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## Loading required package: GenomeInfoDb

``` r
# load gff
lmajor_gff_filepath = file.path(Sys.getenv('REF'), 'lmajor_friedlin', 
                                'annotation', 'TriTrypDB-32_LmajorFriedlin.gff')
lmajor_gff = import.gff3(lmajor_gff_filepath)

# get exons
lmajor_exons = lmajor_gff[lmajor_gff$type == 'exon']

# find genes with more than one exon
multiexons = substring(lmajor_exons[grepl('-2$', lmajor_exons$ID)]$ID, 6, 17)
print(multiexons)
```

    ## character(0)

*T. cruzi* genes that are likely affected (i.e. contain introns):

``` r
# load gff
tcruzi_gff_filepath = file.path(Sys.getenv('REF'), 'tcruzi_clbrener_esmeraldo-like/annotation', 
                                'TriTrypDB-32_TcruziCLBrenerEsmeraldo-like.gff')
tcruzi_gff = import.gff3(tcruzi_gff_filepath)

# get exons
tcruzi_exons = tcruzi_gff[tcruzi_gff$type == 'exon']

# find genes with more than one exon
multiexons = substring(tcruzi_exons[grepl('-2$', tcruzi_exons$ID)]$ID, 6)
print(multiexons)
```

    ## character(0)

Finally, it appears that what goseq calls "gene length" is the median of all transcripts for a gene. From the goseq vignette:

> Once you have a transcriptDb object, you can get a vector named by gene ID containing the median transcript length of each gene simply by using the command.
>
> > txsByGene=transcriptsBy(txdb,"gene")
>
> > lengthData=median(width(txsByGene))

Methods
=======

L. major
--------

First, let's look at what information is available for the generated *L. major* transcript databases.

**Note 2015/02/21**:

Tested both Bioconductor 3.0 and 3.1 (devel), using either the default arguments for `makeTranscriptDbFromGFF` / `makeTxDbFromGFF`, or by specifying either 'gffTxName="gene"' or 'useGenesAsTranscripts=TRUE'.

-   The ideal combination is to use Bioconductor 3.1 with the gffTxName=gene specified.
-   This results both in the transcripts being properly populated from multi-exon genes, and noncoding RNAs being parsed. The downside is that R-devel must be used and many basic packages such as rmarkdown are not yet available. (See note below though...)
-   When using Bioconductor 3.0, the `useGenesAsTranscripts` switch should be enabled to include ncRNAs, however, no settings tested will result in multi-exon genes being properly handled.

NOTE 2015/02/21 -- Previously it was possible to use bioc-devel to produce transcript databases with the CDSSTART and CDSEND fields properly populated, which could be used to determine the processed mRNA length. Currently, however, I am unable to reproduce this and can only get the TX values.

### Coding RNA

``` r
library(Leishmania.major.Friedlin)
```

    ## Loading required package: AnnotationDbi

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: OrganismDbi

    ## Loading required package: GenomicFeatures

    ## Loading required package: GO.db

    ## 

    ## Loading required package: org.LmjF.tritryp.db

    ## 

    ## Loading required package: TxDb.LmajorFriedlin.tritryp32.genes

``` r
orgdb = Leishmania.major.Friedlin

# total number of genes in database
length(keys(TxDb.LmajorFriedlin.tritryp32.genes))
```

    ## [1] 9378

``` r
# L. major gene with introns
gene_id = 'LmjF.02.0100'

# transcript boundaries
select(orgdb, keys=c(gene_id), keytype='GID', columns=c('TXSTART', 'TXEND'))
```

    ## 'select()' returned 1:1 mapping between keys and columns

    ##            GID TXSTART TXEND
    ## 1 LmjF.02.0100   33424 36189

``` r
# CDS boundaries
select(orgdb, keys=c(gene_id), keytype='GID', columns=c('CDSSTART', 'CDSEND'))
```

    ## 'select()' returned 1:many mapping between keys and columns

    ##            GID CDSSTART CDSEND
    ## 1 LmjF.02.0100    35740  36189
    ## 2 LmjF.02.0100    34455  35738
    ## 3 LmjF.02.0100    33424  34446

``` r
# Exon boundaries
select(orgdb, keys=c(gene_id), keytype='GID', columns=c('EXONSTART', 'EXONEND'))
```

    ## 'select()' returned 1:many mapping between keys and columns

    ##            GID EXONSTART EXONEND
    ## 1 LmjF.02.0100     35740   36189
    ## 2 LmjF.02.0100     34455   35738
    ## 3 LmjF.02.0100     33424   34446

``` r
# Transcript length excluding introns
# @TODO: Note that only in Bioconductor devel does the generated transcript
# database include proper TXSTART and TXEND entries; for the current stable
# (bioc 3.0) these fields appear as NAs.
```

### Non-coding RNA

Only TXSTART and TXEND are defined (not CDSSTART/CDSEND).

``` r
# transcript boundaries
gene_id = 'LmjF.02.ncRNA1'

select(orgdb, keys=c(gene_id), keytype='GID', columns=c('TXSTART', 'TXEND'))
```

    ## 'select()' returned 1:1 mapping between keys and columns

    ##              GID TXSTART TXEND
    ## 1 LmjF.02.ncRNA1      NA    NA

``` r
# CDS boundaries
select(orgdb, keys=c(gene_id), keytype='GID', columns=c('CDSSTART', 'CDSEND'))
```

    ## 'select()' returned 1:1 mapping between keys and columns

    ##              GID CDSSTART CDSEND
    ## 1 LmjF.02.ncRNA1       NA     NA

Human
-----

To ensure that the transcript lengths computed match what is expected by `goseq`, let's also compare the approach to the numbers in the gene length database used by goseq for a human gene.

As an example, we will look at the [HOXA10 gene](http://uswest.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000253293;r=7:27170591-27180261) which six known transcripts including [HOXA10-001](http://uswest.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=ENSG00000253293;r=7:27170591-27180261;t=ENST00000283921), which has two exons and is 2541bp long.

``` r
library("Homo.sapiens")
```

    ## Loading required package: org.Hs.eg.db

    ## 

    ## Loading required package: TxDb.Hsapiens.UCSC.hg19.knownGene

``` r
library("geneLenDataBase")
library("GenomeGraphs")
```

    ## Loading required package: biomaRt

    ## Loading required package: grid

``` r
# Load BioMart (used to draw transcript isoforms)
mart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

# target gene (HOXA10) GRCh38 identifiers
gene_id = 'ENSG00000253293'
tx_id   = 'ENST00000283921'

# GRCh37 (hg19) identifiers
gene_id_old = 'ENSG00000153807'

# Plot HOXA10 transcripts
gene = makeGene(id=gene_id, type="ensembl_gene_id", biomart=mart)
transcript = makeTranscript(id=gene_id, type="ensembl_gene_id", biomart=mart,
                            dp=DisplayPars(plotId=TRUE, cex=0.5))
gdPlot(list('Gene'=gene, 'Transcripts'=transcript))
```

![](README_files/figure-markdown_github/human_genelendb-1.png)

``` r
# Query all transcripts from gene length database
data(hg19.ensGene.LENGTH)
hg19.ensGene.LENGTH[hg19.ensGene.LENGTH$Gene == gene_id_old,]
```

    ##                  Gene      Transcript Length
    ## 40709 ENSG00000153807 ENST00000381834   2178
    ## 40710 ENSG00000153807 ENST00000421352   2491
    ## 40711 ENSG00000153807 ENST00000283921   2572
    ## 40712 ENSG00000153807 ENST00000396344   2196

``` r
# HXA10-001
# The length listed here is 2572, which differs from what is listed on the 
# Ensembl website...
hg19.ensGene.LENGTH[hg19.ensGene.LENGTH$Transcript == tx_id,]
```

    ##                  Gene      Transcript Length
    ## 40711 ENSG00000153807 ENST00000283921   2572

#### Gene

``` r
# Coordinates from Homo.sapiens database
orgdb = Homo.sapiens

tx = select(orgdb, keys=c(gene_id), keytype='ENSEMBL',
       columns=c('TXSTART', 'TXEND'))
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
print(tx)
```

    ##           ENSEMBL  TXSTART    TXEND
    ## 1 ENSG00000253293 27210210 27213955
    ## 2 ENSG00000253293 27210210 27219880

``` r
# CDS boundaries
select(orgdb, keys=c(gene_id), keytype='ENSEMBL',
       columns=c('CDSSTART', 'CDSEND'))
```

    ## 'select()' returned 1:many mapping between keys and columns

    ##           ENSEMBL CDSSTART   CDSEND
    ## 1 ENSG00000253293 27212968 27213925
    ## 2 ENSG00000253293 27211518 27211792

``` r
# Exon boundaries
select(orgdb, keys=c(gene_id), keytype='ENSEMBL',
       columns=c('EXONSTART', 'EXONEND'))
```

    ## 'select()' returned 1:many mapping between keys and columns

    ##           ENSEMBL EXONSTART  EXONEND
    ## 1 ENSG00000253293  27212968 27213955
    ## 2 ENSG00000253293  27210210 27211792
    ## 3 ENSG00000253293  27219265 27219880

``` r
# Median transcript length
lengths = c()
for (i in 1:nrow(tx)) {
    lengths = append(lengths, abs(tx[i,]$TXEND - tx[i,]$TXSTART) + 1)
}
print(sprintf("Median transcript length: %01f", median(lengths)))
```

    ## [1] "Median transcript length: 6708.500000"

#### Transcript

``` r
# Using transcript ID to query
select(orgdb, keys=c(tx_id), keytype='ENSEMBLTRANS',
       columns=c('TXSTART', 'TXEND'))
```

    ## 'select()' returned 1:many mapping between keys and columns

    ##      ENSEMBLTRANS  TXSTART    TXEND
    ## 1 ENST00000283921 27210210 27213955
    ## 2 ENST00000283921 27210210 27219880

``` r
# CDS boundaries
select(orgdb, keys=c(tx_id), keytype='ENSEMBLTRANS',
       columns=c('CDSSTART', 'CDSEND'))
```

    ## 'select()' returned 1:many mapping between keys and columns

    ##      ENSEMBLTRANS CDSSTART   CDSEND
    ## 1 ENST00000283921 27212968 27213925
    ## 2 ENST00000283921 27211518 27211792

``` r
# Exon boundaries
select(orgdb, keys=c(tx_id), keytype='ENSEMBLTRANS',
       columns=c('EXONSTART', 'EXONEND'))
```

    ## 'select()' returned 1:many mapping between keys and columns

    ##      ENSEMBLTRANS EXONSTART  EXONEND
    ## 1 ENST00000283921  27212968 27213955
    ## 2 ENST00000283921  27210210 27211792
    ## 3 ENST00000283921  27219265 27219880

System Info
-----------

``` r
sessionInfo()
```

    ## R version 3.4.2 (2017-09-28)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Arch Linux
    ## 
    ## Matrix products: default
    ## BLAS: /usr/lib/libblas.so.3.7.1
    ## LAPACK: /usr/lib/liblapack.so.3.7.1
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ##  [1] grid      parallel  stats4    stats     graphics  grDevices utils    
    ##  [8] datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] GenomeGraphs_1.38.0                     
    ##  [2] biomaRt_2.34.0                          
    ##  [3] geneLenDataBase_1.14.0                  
    ##  [4] Homo.sapiens_1.3.1                      
    ##  [5] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 
    ##  [6] org.Hs.eg.db_3.4.2                      
    ##  [7] Leishmania.major.Friedlin_32.0          
    ##  [8] TxDb.LmajorFriedlin.tritryp32.genes_32.0
    ##  [9] org.LmjF.tritryp.db_32.0                
    ## [10] GO.db_3.4.2                             
    ## [11] OrganismDbi_1.20.0                      
    ## [12] GenomicFeatures_1.30.0                  
    ## [13] AnnotationDbi_1.40.0                    
    ## [14] Biobase_2.38.0                          
    ## [15] rtracklayer_1.38.0                      
    ## [16] GenomicRanges_1.30.0                    
    ## [17] GenomeInfoDb_1.14.0                     
    ## [18] IRanges_2.12.0                          
    ## [19] S4Vectors_0.16.0                        
    ## [20] BiocGenerics_0.24.0                     
    ## [21] knitr_1.17                              
    ## [22] rmarkdown_1.6                           
    ## [23] nvimcom_0.9-40                          
    ## [24] colorout_1.1-3                          
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] SummarizedExperiment_1.8.0 progress_1.1.2            
    ##  [3] lattice_0.20-35            htmltools_0.3.6           
    ##  [5] yaml_2.1.14                blob_1.1.0                
    ##  [7] XML_3.98-1.9               RBGL_1.54.0               
    ##  [9] rlang_0.1.2                DBI_0.7                   
    ## [11] BiocParallel_1.12.0        bit64_0.9-7               
    ## [13] matrixStats_0.52.2         GenomeInfoDbData_0.99.1   
    ## [15] stringr_1.2.0              zlibbioc_1.24.0           
    ## [17] Biostrings_2.46.0          evaluate_0.10.1           
    ## [19] memoise_1.1.0              BiocInstaller_1.28.0      
    ## [21] Rcpp_0.12.13               backports_1.1.1           
    ## [23] DelayedArray_0.4.0         graph_1.56.0              
    ## [25] XVector_0.18.0             bit_1.1-12                
    ## [27] Rsamtools_1.30.0           RMySQL_0.10.13            
    ## [29] digest_0.6.12              stringi_1.1.5             
    ## [31] rprojroot_1.2              tools_3.4.2               
    ## [33] bitops_1.0-6               magrittr_1.5              
    ## [35] RCurl_1.95-4.8             RSQLite_2.0               
    ## [37] tibble_1.3.4               pkgconfig_2.0.1           
    ## [39] Matrix_1.2-11              prettyunits_1.0.2         
    ## [41] assertthat_0.2.0           R6_2.2.2                  
    ## [43] GenomicAlignments_1.14.0   compiler_3.4.2
