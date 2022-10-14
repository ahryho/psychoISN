# __Individual-Specific Network for Stress-Related Psychiatric Disorders__

_by Anastasiia Hryhorzhevska_


[Issues](#issues) 

[Brief introduction to ISNs](#brief-introduction-to-isns)

[1. Input data](#1-input-data)

## **1. Input data:**

Data are stored on the MPIP computational cluster:

- DNAm: 

The data are stored at `/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/input/dnam`.

- Genotype: 

  + __196 samples__
  + __3'957'337 SNPs__ after QC
  + __3'908'485 SNPs__ after QC and MAF filtering

The data are stored at `/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/input/snps`.

### QC Roadmap

__The roadmap for imputed genotype QC analysis__

1. Subset data:
    + MAF >= 5%
    + filter only SNPs
    + exclude MHC region