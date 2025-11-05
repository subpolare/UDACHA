home='/home/subpolare/adastra-v7'
threads=10

# DO NOT EDIT BELOW

mkdir -p tmp/ VCFs/ BEDs/ BADs/ SNPs/ SNPScan/ mixalime/
if [ ! -d hocomoco/v13/pwm ]; then
    set -euo pipefail
    mkdir -p hocomoco/v13/pwm
    curl -fsSL "https://hocomoco13.autosome.org/final_bundle/hocomoco13/H13RSNP/H13RSNP_pwm.tar.gz" \
        | tar -xz -C hocomoco/v13/pwm --strip-components=2
fi

# 1. Merging files

python3 ${home}/scripts/renamer.py 
for file in ${home}/VCFs/*.vcf.gz; do
    bcftools index -f $file
done

find ${home}/VCFs -maxdepth 1 -type f -name "*.vcf.gz" ! -name "*.without_MAF.vcf.gz" -print0 |
    while IFS= read -r -d '' f; do
        out="${f%.vcf.gz}.without_MAF.vcf.gz"
        bcftools annotate -x INFO/MAF -Oz -o "$out" "$f"
        bcftools index -f "$out"
    done

find ${home}/VCFs -type f -name "*.without_MAF.vcf.gz" | sort > ${home}/tmp/vcfs.list

echo -e "sample_id\ttf\tcell\talgn_id\tgse\tpath" > ${home}/tmp/samples.meta.tsv
while IFS= read -r f; do
    b=$(basename "$f")
    base=${b%.vcf.gz}
    IFS=_ read -r tf cell fileid gse <<< "$base"
    sid=$(bcftools query -l "$f")
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$sid" "$tf" "$cell" "$fileid" "$gse" "$f" >> ${home}/tmp/samples.meta.tsv
done < ${home}/tmp/vcfs.list

bcftools merge -m none --threads $threads -Oz -o ${home}/VCFs/merged.without_MAF.vcf.gz -l ${home}/tmp/vcfs.list

# If there is an duplicate error, use command below to find it: 
# for f in VCFs/*.vcf.gz; do bcftools query -l "$f"; done | sort | uniq -d
# than you need to remove duplicates 

bcftools index -f ${home}/VCFs/merged.without_MAF.vcf.gz











# NOT UPDATED || OLD VERSION

# 2. BABACHI, https://github.com/autosome-ru/BABACHI 

run_babachi() { 
    file=$1
    home='/home/subpolare/adastra-v7'
    name=$(basename $file .vcf.gz)

    if [ $(zcat $file | grep '^chr' | wc -l) -ne 0 ]; then
        python3 ${home}/scripts/create_bed.py -i $file -o ${home}/BEDs/${name}.bed
        sed -i '/\./d' ${home}/BEDs/${name}.bed
        babachi ${home}/BEDs/${name}.bed -O ${home}/BADs/
        find ${home}/BADs -type f -exec sh -c 'for f in "$@"; do if [ "$(wc -l < "$f")" -eq 1 ]; then rm "$f"; fi; done' sh {} +
        # babachi visualize ${home}/BEDs/${name}.bed -O ${home}/BADs/ -b ${home}/BADs/${name}.badmap.bed
        # python3 ${home}/scripts/svg2png.py -d ${home}/BADs/${name}.badmap.visualization
        if [ -e "${home}/BADs/${name}.badmap.bed" ]; then
            python3 ${home}/scripts/add_bad_to_bed.py \
                --bed ${home}/BEDs/${name}.bed \
                --bad ${home}/BADs/${name}.badmap.bed \
                --output ${home}/BEDs/${name}.with_bad.bed
        fi
        rm ${home}/BEDs/${name}.bed # ${home}/BADs/${name}.badmap.visualization/*.svg
    fi
}
export -f run_babachi
find ${home}/VCFs -name *.vcf.gz | parallel -j $threads run_babachi

# 3. MixALiMe for TFs, http://mixalime.georgy.top/tutorial/quickstart.html

for model in MCNB NB BetaNB; do
    mixalime create ${home}/mixalime/${model} ${home}/BEDs/*.with_bad.bed --no-snp-bad-check
    mixalime fit ${home}/mixalime/${model} $model
    mixalime test ${home}/mixalime/${model}

    TFs=$(for file in ${home}/BEDs/*.with_bad.bed; do
        basename "$file" | cut -d'_' -f1
    done)

    for TF in $(echo "$TFs" | sort -u); do
        mixalime combine ${home}/mixalime/${model} --group "m:${home}/BEDs/${TF}_*.with_bad.bed" --subname $TF
    done

    mixalime export all ${home}/mixalime/${model} ${home}/mixalime/results_${model}
    mixalime plot all ${home}/mixalime/${model} ${home}/mixalime/results_${model}
done

# 4. Creating tables for TFs

for model in MCNB NB BetaNB; do
    for TF in $(echo "$TFs" | sort -u); do 
        python3 ${home}/scripts/create_tf_tables.py \
            --mixalime ${home}/mixalime/results_${model}/pvalues/${TF}.tsv \
            --bed "${home}/BEDs/${TF}*.with_bad.bed" \
            --output ${home}/new-version/TF/${TF}_HUMAN_${model}.tsv
    done
done

# 5. Motif annotation of TF tables

for file in $(ls -1 ${home}/new-version/TF/*); do
    TF=$(basename $file | cut -f1 -d '.' | cut -f1 -d '_')
    python3 scripts/make_snps_list.py \
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
