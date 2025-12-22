<div align="center">
  <h1>Updated UDACHA-workflow</h1> 
</div>

Continuation of the [updated ADASTRA workflow](https://github.com/subpolare/ADASTRA) based on the [article at the link](https://www.nature.com/articles/s41467-024-55513-2). 

The code takes as input GATK/bcftools-derived per-sample bgzipped, tabix-indexed VCFs of heterozygous SNVs (with GT/AD and coverage filters, aligned to GRCh38/dbSNP151) and outputs an [updated ADASTRA database](https://adastra.autosome.org/) with recomputed genotype clusters, BAD maps, and allele-specific variant calls.

![python](https://img.shields.io/badge/python%20-%234584B6.svg?&style=for-the-badge&logo=python&logoColor=white) ![pandas](https://img.shields.io/badge/pandas%20-%23150458.svg?&style=for-the-badge&logo=pandas&logoColor=white) ![numpy](https://img.shields.io/badge/numpy-%23013243.svg?&style=for-the-badge&logo=numpy&logoColor=white) ![shell](https://img.shields.io/badge/shell-%234EAA25.svg?&style=for-the-badge&logo=gnu-bash&logoColor=white)

## 🌊 How to run 

1. Install the necessary packages: 

```
pip install babachi 
pip install mixalime  
conda install bioconda::plink2
pip install tqdm
```

2. Install [bcftools](https://www.htslib.org/download/) according to the instructions from the developers.

3. Copy resulting `.vcf.gz` files from the [updated ADASTRA workflow](https://github.com/subpolare/ADASTRA) to the `VCFs/`.

4. Update following variables in the beginning of the `run.sh`:

  - home: full path to the working directory
  - scripts: full path to the [scripts directory](https://github.com/subpolare/UDACHA/tree/main/scripts)
  - threads: number of worker threads 

6. Run `run.sh`
