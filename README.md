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

3. Create following directories: 

```
tmp/
VCFs/
BEDs/
BADs/
SNPs/
SNPScan/
mixalime/
hocomoco/v13/
```

4. Download [hocomococo v13](https://hocomoco13.autosome.org/downloads_v13) (PWMs from the H13RSNP subcollection) to the `hocomoco/v13/`.

5. Copy resulting `.vcf.gz` files from the [updated ADASTRA workflow](https://github.com/subpolare/ADASTRA) to the `VCFs/`.

6. Run `run.sh` with your `$home` (workdir) and `$threads` in the beginning of the file. 
