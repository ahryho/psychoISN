# __Uncovering relationships between stress and psychiatric traits via multi-omics individual-specific networks (ISNs)__

_by Anastasiia Hryhorzhevska_

[1. Brief introduction to ISNs](#brief-introduction-to-isns)

[2. Input data](#input-data)

[3. Research questions](#research-questions)

[4. Results and aditional information](#results-and-additional-information)

## **Brief introduction to ISNs**

**What?**

Individual-Specific Network (ISN) is a network which is based on multiple measurements for the same individual. In the present project, we concentrate in the individual-specific edges, i.e. the node information is available or not.  

**Why?**

ISNs are useful because 

1. A network derived from a collection of individuals can be used as a model for an “average” individual. Meaning that by doing ISN, one can translate the network interpretation strategies from the population to individual, i.e. perform extrapolation to the level of the individual. 

2. ISNs allow focusing on each individual and their specific associations and dynamics over time.

**How?**

The core of the approach is following: 

![isns_construction](https://github.com/ahryho/dex-stim-human-array-isns/blob/main/materials/figures/isns_construction.jpg)

For additional information on methodology please refer to the [slides](https://github.com/ahryho/dex-stim-human-array-isns/blob/main/materials/final_presentation.pdf).

## **Input data**

Participants comprised **196** Caucasian individuals of the Max Planck Institute of Psychiatry (MPIP) in Munich. Of the participants, 

+ 131 men and 65 women
+ 84 (50 men, 34 women) were treated for major dipressive disorders (MDD)
+ 112 (81 men, 31 women) were healthy controls with no history of a depressive disorder. 

Baseline whole blood samples were obtained at 6 pm after two hours of fasting and abstention from coffee and physical activity. Subjects then received 1.5 mg oral dexamethasone, and a second blood draw was performed at 9 pm, three hours after dexamethasone ingestion.

The available multiomics data for the currect project:

- Methylation
- Gene-expression (for further analysis)
- Genotype
- Phenotype 

For detailed data overview, please refer to the [data overview](https://github.com/ahryho/psychoGE/blob/master/code/integrative/data_overview/01_data_overview.html).

Data are stored on the MPIP computational cluster:

- DNAm: `/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/input/dnam`

- Genotype: `/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/input/snps`

- Phenotype: `/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/input/pheno`.

## Research questions

1. How are the individuals grouped according to the changes before and after stimulation with dexamethasone?

2. What are the main drivers of these changes?

## Results and additional information

Results are stored on the MPIP computational cluster: `/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/results`.

For plots which describe results, please refer to the [folder](https://github.com/ahryho/psychoISNs/tree/main/results/plots).

For additional information on methodology, results and limitations, please refer to the [slides](https://github.com/ahryho/dex-stim-human-array-isns/blob/main/materials/final_presentation.pdf).