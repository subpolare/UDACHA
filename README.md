<div align="center">
  <h1>Updated UDACHA-workflow</h1> 
</div>

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
