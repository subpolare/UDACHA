<div align="center">
  <h1>Updated UDACHA-workflow</h1> 
</div>

Continuation of the [updated ADASTRA workflow](https://github.com/subpolare/ADASTRA) based on the [article at the link](https://www.nature.com/articles/s41467-024-55513-2). The script takes as input GATK/bcftools-derived per-sample bgzipped, tabix-indexed VCFs of heterozygous SNVs (with GT/AD and coverage filters, aligned to GRCh38/dbSNP151) and outputs an [updated ADASTRA database](https://adastra.autosome.org/mabel) with recomputed genotype clusters, BAD maps, and allele-specific variant calls.

## 🌊 How to run 

1. Install the necessary packages: 

```
pip install babachi 
pip install mixalime  
conda install bioconda::plink2
```

2. Install [bcftools](https://www.htslib.org/download/) according to the instructions from the developers.

3. Copy resulting `.vcf.gz` files from the [updated ADASTRA workflow](https://github.com/subpolare/ADASTRA) to the `VCFs/`.

4. Run `run.sh` with your `$home` (workdir) and `$threads` in the beginning of the file. 
