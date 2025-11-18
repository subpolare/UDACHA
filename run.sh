home='/sandbox/subpolare/adastra'
scripts='/home/subpolare/adastra-v7/scripts'
threads=50

# DO NOT EDIT BELOW

export home 
export scripts 
export threads

mkdir -p ${home}/VCFs/ ${home}/BEDs/ ${home}/clustering ${home}/BADs/ ${home}/SNPs/ ${home}/SNPScan/ ${home}/mixalime/ ${home}/mixalime/groups/
if [ ! -d ${home}/hocomoco/v13/pwm ]; then
    set -euo pipefail
    mkdir -p ${home}/hocomoco/v13/pwm
    curl -fsSL "https://hocomoco13.autosome.org/final_bundle/hocomoco13/H13RSNP/H13RSNP_pwm.tar.gz" \
        | tar -xz -C ${home}/hocomoco/v13/pwm --strip-components=2
fi

# 1. Merging files & Clustering

python3 ${scripts}/clustering/renamer.py 
for file in ${home}/VCFs/*.vcf.gz; do
    bcftools index -f $file
done

find ${home}/VCFs -maxdepth 1 -type f -name "*.vcf.gz" ! -name "*.without_MAF.vcf.gz" -print0 |
    while IFS= read -r -d '' f; do
        out="${f%.vcf.gz}.without_MAF.vcf.gz"
        bcftools annotate -x INFO/MAF -Oz -o "$out" "$f"
        bcftools index -f "$out"
    done

find ${home}/VCFs -type f -name "*.without_MAF.vcf.gz" | sort > ${home}/clustering/vcfs.list
echo -e "indiv_id\ttf\tcell\talgn_id\tgse\tpath" > ${home}/clustering/samples.meta.tsv
while IFS= read -r f; do
    b=$(basename "$f")
    base=${b%.vcf.gz}
    IFS=_ read -r tf cell fileid gse <<< "$base"
    sid=$(bcftools query -l "$f")
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$sid" "$tf" "$cell" "$fileid" "$gse" "$f" >> ${home}/clustering/samples.meta.tsv
done < ${home}/clustering/vcfs.list

bcftools merge -m none --threads $threads -Oz -o ${home}/VCFs/merged.without_MAF.vcf.gz -l ${home}/clustering/vcfs.list 
bcftools index -f ${home}/VCFs/merged.without_MAF.vcf.gz --threads $threads

# If there is an duplicate error in bcftools merge, use command below to find duplicates: 
# for f in VCFs/*.vcf.gz; do bcftools query -l "$f"; done | sort | uniq -d
# of all the duplicates, only one should be left 

plink2 --vcf ${home}/VCFs/merged.without_MAF.vcf.gz \
  --allow-extra-chr \
  --threads $threads \
  --make-king square \
  --out ${home}/clustering/king_all

python3 ${scripts}/clustering/clustering.py \
  --matrix     ${home}/clustering/king_all.king \
  --matrix-ids ${home}/clustering/king_all.king.id \
  --meta-file  ${home}/clustering/samples.meta.tsv \
  --outpath    ${home}/clustering

python3 ${scripts}/clustering/create_bed_clusters.py \
  --metadata ${home}/clustering/metadata.clustered.tsv \
  --work     ${home}

find ${home}/BEDs -type f -name '*.bed' -exec sh -c '
  for f do
    if [ "$(wc -l < "$f")" -eq 1 ]; then
      rm "$f"
    fi
  done
' sh {} +

# 2. BABACHI, https://github.com/autosome-ru/BABACHI 

run_babachi() { 
    name=$(basename $1 .bed)
    babachi ${home}/BEDs/${name}.bed -O ${home}/BADs/ >/dev/null 2>&1

    if [ "$(wc -l < "${home}/BADs/${name}.badmap.bed")" -gt 1 ]; then
        echo ${home}/BADs/${name}.badmap.bed
        babachi visualize ${home}/BEDs/${name}.bed -O ${home}/BADs/ -b ${home}/BADs/${name}.badmap.bed
        python3 ${scripts}/babachi/svg2png.py -d ${home}/BADs/${name}.badmap.visualization
        rm ${home}/BADs/${name}.badmap.visualization/*.svg
    fi    
    python3 ${scripts}/babachi/add_bad_to_bed.py \
        --bed    ${home}/BEDs/${name}.bed \
        --bad    ${home}/BADs/${name}.badmap.bed \
        --output ${home}/BEDs/${name}.with_bad.bed
    rm ${home}/BEDs/${name}.bed
}

export -f run_babachi
find ${home}/BEDs -name *.bed | parallel -j 10 run_babachi {} 

# 3. MixALiMe for TFs, http://mixalime.georgy.top/tutorial/quickstart.html

cut -f2 ${home}/clustering/metadata.clustered.tsv | tail -n +2 | sort -u > ${home}/mixalime/groups/factors.list
while read tf; do
    awk -F'\t' -v tf="$tf" 'NR > 1 && $2 == tf { print $7 }' ${home}/clustering/metadata.clustered.tsv \
        | sort -u \
        | awk -v home="$home" '{ printf "%s/BEDs/%s.with_bad.bed\n", home, $1 }' \
        > ${home}/mixalime/groups/factors_"$tf".list
done < ${home}/mixalime/groups/factors.list

cut -f3 ${home}/clustering/metadata.clustered.tsv | tail -n +2 | sort -u > ${home}/mixalime/groups/cell.list
while read cell; do
    awk -F'\t' -v cell="$cell" 'NR > 1 && $3 == cell { print $7 }' ${home}/clustering/metadata.clustered.tsv \
        | sort -u \
        | awk -v home="$home" '{ printf "%s/BEDs/%s.with_bad.bed\n", home, $1 }' \
        > ${home}/mixalime/groups/cell_${cell}.list
done < ${home}/mixalime/groups/cell.list

for model in MCNB NB BetaNB; do
    project=${home}/mixalime/${model}

    mixalime create $project ${home}/BEDs/*.with_bad.bed --no-snp-bad-check
    mixalime fit $project $model
    mixalime test $project

    while read tf; do
        mixalime combine \
            --subname "TF_${tf}" \
            --group ${home}/mixalime/groups/factors_${tf}.list \
            $project
    done < ${home}/mixalime/groups/factors.list

    while read cell; do
        mixalime combine \
            --subname "CELL_${cell}" \
            --group ${home}/mixalime/groups/cell_${cell}.list \
            $project
    done < ${home}/mixalime/groups/cell.list

    mixalime export all $project ${home}/mixalime/results_${model}
    mixalime plot all $project ${home}/mixalime/results_${model}
done










# NOT UPDATED || OLD VERSION

# 4. Creating tables for TFs

for model in MCNB NB BetaNB; do
    for TF in $(echo "$TFs" | sort -u); do 
        python3 ${scripts}/create_tf_tables.py \
            --mixalime ${home}/mixalime/results_${model}/pvalues/${TF}.tsv \
            --bed "${home}/BEDs/${TF}*.with_bad.bed" \
            --output ${home}/new-version/TF/${TF}_HUMAN_${model}.tsv
    done
done

# 5. Motif annotation of TF tables

for file in $(ls -1 ${home}/new-version/TF/*); do
    TF=$(basename $file | cut -f1 -d '.' | cut -f1 -d '_')
    python3 ${scripts}/make_snps_list.py \
        --genome '/home/subpolare/genome/GRCh38.primary_assembly.genome.fa' \
        --threads 20 \
        --input ${home}/new-version/TF/${TF}_HUMAN.tsv \
        > ${home}/SNPs/${TF}_HUMAN.snps
done

run_SNPScan() {
    home='/home/subpolare/adastra-v7'
    factor=$(basename $2 | cut -f1 -d '.')
    java -cp ${home}/scripts/ape-3.0.6.jar ru.autosome.perfectosape.SNPScan \
        ${home}/hocomoco/v12/pwm/${factor}.H12RSNP.${1}.*.pwm $2 \
        --single-motif -F 1 -P 1 -d 1000 \
        > ${home}/SNPScan/pwm_results_${1}/${factor}.perfectos
}
export -f run_SNPScan 
parallel -j $threads run_SNPScan ::: 0 1 2 3 ::: $(ls -1 ${home}/SNPs/*)
for file in /home/subpolare/adastra-v7/SNPScan/pwm_results_?/*; do sed -i '1s/^# //' "$file"; done
find ${home}/SNPScan/ -size 0 -delete

python3 ${home}/scripts/merge_snpscan_results.py
python3 ${home}/scripts/update_tf_tables.py

for TF in $(echo "$TFs" | sort -u); do
    python3 ${home}/scripts/add_raw_pvalue.py \
        --mixalime ${home}/mixalime/results_${model}/pvalues/${TF}.tsv \
        --adastra ${home}/new-version/TF/${TF}_HUMAN.tsv
done
