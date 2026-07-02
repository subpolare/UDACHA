scripts='/home/subpolare/adastra-v7/scripts'
home='/sandbox/subpolare/adastra'
threads=50

mkdir -p ${home}/VCFs/ ${home}/BEDs/ ${home}/clustering ${home}/BADs/ ${home}/SNPs/ ${home}/SNPScan/ ${home}/mixalime/ ${home}/mixalime/groups/ ${home}/logs ${home}/mixalime/file_lists/cells_500K/
if [ ! -d ${home}/hocomoco/v13/pwm ]; then
    set -euo pipefail
    mkdir -p ${home}/hocomoco/v13/pwm
    curl -fsSL "https://hocomoco13.autosome.org/final_bundle/hocomoco13/H13RSNP/H13RSNP_pwm.tar.gz" \
        | tar -xz -C ${home}/hocomoco/v13/pwm --strip-components=2
fi

# Merging files 

python3 ${scripts}/clustering/renamer.py 
for file in ${home}/VCFs/*.vcf.gz; do
    bcftools index --threads $threads -f $file
done

find ${home}/VCFs -maxdepth 1 -type f -name "*.vcf.gz" ! -name "*.without_MAF.vcf.gz" -print0 |
    while IFS= read -r -d '' f; do
        out="${f%.vcf.gz}.without_MAF.vcf.gz"
        bcftools annotate --threads $threads -x INFO/MAF -Oz -o "$out" "$f"
        bcftools index --threads $threads -f "$out"
    done

find ${home}/VCFs -type f -name "*.without_MAF.vcf.gz" | sort > ${home}/clustering/merged.list
echo -e "indiv_id\ttf\tcell\talgn_id\tgse\tpath" > ${home}/clustering/samples.meta.tsv
while IFS= read -r f; do
    b=$(basename "$f")
    base=${b%.without_MAF.vcf.gz}
    IFS=_ read -r tf cell fileid gse <<< "$base"
    sid=$(bcftools query -l "$f")
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$sid" "$tf" "$cell" "$fileid" "$gse" "$f" >> ${home}/clustering/samples.meta.tsv
done < ${home}/clustering/merged.list

find ${home}/VCFs -maxdepth 1 -type f -name "*.without_MAF.vcf.gz" ! -name "merged*" -print0 \
    | xargs -0 -n 1 -P "$threads" bash -c '
        f="$1"
        n="$(bcftools view --no-version -v snps -H "$f" 2>/dev/null | head -n 100 | wc -l || true)"
        if [ "$n" -ge "$cutoff" ]; then
            printf "%s\n" "$f"
        fi
      ' _ \
    | LC_ALL=C sort -u > ${home}/clustering/merged.min100.list

bcftools merge -l ${home}/clustering/merged.min100.list \
    --missing-to-ref \
    -Oz -o ${home}/VCFs/merged.min100.vcf.gz \ 
    --threads $threads \
    -m none 

bcftools index --threads $threads -f ${home}/VCFs/merged.min100.vcf.gz

# 2. Hierarchical clustering using PLINK2 data 

plink2 --vcf ${home}/VCFs/merged.min100.vcf.gz \
  --allow-extra-chr \
  --threads $threads \
  --make-king square \
  --out ${home}/clustering/king_min100

python3 ${scripts}/clustering/clustering.py \
  --king ${home}/clustering/king_min100.king \
  --king-id ${home}/clustering/king_min100.king.id \
  --meta ${home}/clustering/samples.meta.tsv \
  --out ${home}/clustering/metadata.clustered.tsv \
  --floor 0.0 \
  --thr 0.8877 \
  --method complete

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

# BABACHI, https://github.com/autosome-ru/BABACHI 

find ${home}/BEDs -maxdepth 1 -name 'INDIV_*.bed' -print0 | while IFS= read -r -d '' file; do
    name=$(basename $file .bed)
    babachi ${home}/BEDs/${name}.bed -j 25 -p geometric -g 0.99 -s "1,4/3,3/2,2,5/2,3,4,5,6" -O ${home}/BADs/
    if [ "$(wc -l < "${home}/BADs/${name}.badmap.bed")" -gt 1 ]; then
        babachi visualize ${home}/BEDs/${name}.bed -O ${home}/BADs/ -b ${home}/BADs/${name}.badmap.bed
        # python3 ${scripts}/babachi/svg2png.py -d ${home}/BADs/${name}.badmap.visualization
        # rm ${home}/BADs/${name}.badmap.visualization/*.svg
    fi    
    python3 ${scripts}/babachi/add_bad_to_bed.py \
        --bed    ${home}/BEDs/${name}.bed \
        --bad    ${home}/BADs/${name}.badmap.bed \
        --output ${home}/BEDs/${name}.with_bad.bed
done 

# Filtration based on pooled samples, GSE and reads number

mkdir -p ${home}/mixalime/file_lists/
python3 ${scripts}/clustering/get_pooled_from_geo.py 
python3 ${scripts}/meta/join_meta.py \
    --meta6 ~/adastra-v7/meta/meta_6_may.tsv \
    --qc ~/adastra-v7/meta/full_qc_table.tsv \
    --clustered ${home}/clustering/GEO/metadata.clustered.pooled.tsv \
    --out ${home}/meta.tsv

awk -F'\t' '
NR==1 {for(i=1;i<=NF;i++) c[$i]=i; next}
function num(x){return x~/^[0-9.]+([eE][-+]?[0-9]+)?$/}
{
    id=$c["indiv_id"]; if(id=="" || id=="NA") next
    base=id; sub(/__CELL.*/,"",base)
    p=tolower($c["pooled"]); n=$c["NRF"]; r=$c["reads_num"]

    if (p=="true" || (num(n) && n<0.05) || (num(n) && num(r) && n<0.4 && r>0.75e9)) print base
}
' ${home}/meta.tsv | sort -u > ${home}/mixalime/file_lists/filtered_list.txt

# Find BADs with 500 000 SNPs at least and others

find ${home}/BEDs -maxdepth 1 -name 'INDIV_*.with_bad.bed' -print0 \
    | xargs -0 -I{} bash -c 'n=$(($(wc -l < "$1") - 1)); printf "%s\t%s\n" "$(basename "$1" .with_bad.bed)" "$n"' _ {} \
    | sort > ${home}/mixalime/file_lists/indiv_snps.tsv

find ${home}/BEDs -maxdepth 1 -name 'INDIV_*.with_bad.bed' -print0 \
    | xargs -0 -I{} bash -c 'n=$(($(wc -l < "$1") - 1)); [ "$n" -gt 500000 ] && basename "$1" .with_bad.bed' _ {} \
    | sort > ${home}/mixalime/file_lists/halfmillions.txt

find ${home}/BEDs -maxdepth 1 -name 'INDIV_*.with_bad.bed' -print0 \
    | xargs -0 -I{} basename {} .with_bad.bed \
    | sort | grep -Fvx -f ${home}/mixalime/file_lists/halfmillions.txt \
    > ${home}/mixalime/file_lists/not_halfmillions.txt

awk 'NR==FNR{bad[$1];next}{x=$1;sub(/__CELL.*/,"",x)}!(x in bad)&&!seen[$1]++' \
  ${home}/mixalime/file_lists/filtered_list.txt \
  ${home}/mixalime/file_lists/halfmillions.txt \
  > ${home}/mixalime/file_lists/halfmillions.filtered.txt

awk 'NR==FNR{bad[$1];next}{x=$1;sub(/__CELL.*/,"",x)}!(x in bad)&&!seen[$1]++' \
  ${home}/mixalime/file_lists/filtered_list.txt \
  ${home}/mixalime/file_lists/not_halfmillions.txt \
  > ${home}/mixalime/file_lists/not_halfmillions.filtered.txt

awk -F'\t' '
NR==FNR{keep[$1];next}
FNR==1{for(i=1;i<=NF;i++) c[$i]=i; next}
($c["indiv_id"] in keep){
    id=$c["indiv_id"]
    cell=$c["cell_id"]; if(cell==""||cell=="NA") cell=$c["cell"]
    if(cell!=""&&cell!="NA"&&!seen[id]++) print id,cell
}' OFS='\t' \
    ${home}/mixalime/file_lists/not_halfmillions.filtered.txt \
    ${home}/meta.tsv \
    | sort \
    > ${home}/mixalime/file_lists/not_halfmillions.indiv_cell.tsv

awk 'NR==FNR{s[$1]=$2;next}{print $1,$2,s[$1]}' OFS='\t' \
    ${home}/mixalime/file_lists/indiv_snps.tsv \
    ${home}/mixalime/file_lists/not_halfmillions.indiv_cell.tsv \
    > ${home}/mixalime/file_lists/not_halfmillions.indiv_cell_snps.tsv

awk -F'\t' '{s[$2]+=$3} END{for(c in s) if(s[c]>500000) print c}' \
    ${home}/mixalime/file_lists/not_halfmillions.indiv_cell_snps.tsv \
    | sort > ${home}/mixalime/file_lists/cells_500K.list

rm -f ${home}/mixalime/file_lists/cells_500K/CELL_*.txt
rm -f ${home}/mixalime/file_lists/less_than_500K.after_cells.txt

awk -F'\t' -v home="${home}" '
NR==FNR{big[$1];next}
{
    out = ($2 in big) ? home"/mixalime/file_lists/cells_500K/CELL_"$2".txt" : home"/mixalime/file_lists/less_than_500K.after_cells.txt"
    print $1 >> out
}' \
    ${home}/mixalime/file_lists/cells_500K.list \
    ${home}/mixalime/file_lists/not_halfmillions.indiv_cell_snps.tsv

# MixALiMe without final combine, for all clusters with at least 500 000 SNPs after filtration

rm -rf ${home}/mixalime/INDIV_???? ${home}/mixalime/less_than_500K_SNPs

while IFS= read -r indiv_id <&3; do
    mkdir -p ${home}/mixalime/${indiv_id}
    project=${home}/mixalime/${indiv_id}/${indiv_id}
    python3 ${scripts}/mixalime/limiter.py --threads $threads create $project ${home}/BEDs/${indiv_id}.with_bad.bed --no-snp-bad-check --max-cover 10000 
    python3 ${scripts}/mixalime/limiter.py --threads $threads fit $project NB
    python3 ${scripts}/mixalime/limiter.py --threads $threads test $project
    python3 ${scripts}/mixalime/limiter.py --threads $threads combine $project
    python3 ${scripts}/mixalime/limiter.py --threads $threads export all $project $project
    python3 ${scripts}/mixalime/limiter.py --threads $threads plot all $project $project
done 3< ${home}/mixalime/file_lists/halfmillions.filtered.txt

python3 ${scripts}/mixalime/indivs_sclices.py \
    --home ${home} --meta ${home}/meta.tsv \
    --cells-meta /home/subpolare/adastra-v7/meta/meta_cells_and_tissues.tsv

# MixALiMe without final combine, for all cell types with at least 500 000 SNPs for all clusters in this cell type

while IFS= read -r cell_id <&3; do
    list=${home}/mixalime/file_lists/cells_500K/CELL_${cell_id}.txt
    if [ ! -s "$list" ]; then
        continue
    fi

    mkdir -p ${home}/mixalime/CELL_${cell_id}
    project=${home}/mixalime/CELL_${cell_id}/CELL_${cell_id}
    mapfile -t files < <(
        awk -v home="${home}" '{print home "/BEDs/" $0 ".with_bad.bed"}' "$list"
    )

    if [ "${#files[@]}" -gt 0 ]; then
        python3 ${scripts}/mixalime/limiter.py --threads $threads create $project "${files[@]}" --no-snp-bad-check --max-cover 10000 
        python3 ${scripts}/mixalime/limiter.py --threads $threads fit $project NB
        python3 ${scripts}/mixalime/limiter.py --threads $threads test $project
        python3 ${scripts}/mixalime/limiter.py --threads $threads combine $project
        python3 ${scripts}/mixalime/limiter.py --threads $threads export all $project $project
        python3 ${scripts}/mixalime/limiter.py --threads $threads plot all $project $project
    fi
done 3< ${home}/mixalime/file_lists/cells_500K.list
python3 ${scripts}/mixalime/cells_sclices.py --home ${home} \
    --cells-meta /home/subpolare/adastra-v7/meta/meta_cells_and_tissues.tsv

# MixALiMe without final combine, for other clusters

mkdir -p ${home}/mixalime/less_than_500K_SNPs
project=${home}/mixalime/less_than_500K_SNPs/less_than_500K_SNPs
mapfile -t files < <(
    awk -v home="${home}" '{print home "/BEDs/" $0 ".with_bad.bed"}' \
        ${home}/mixalime/file_lists/less_than_500K.after_cells.txt
)
if [ "${#files[@]}" -gt 0 ]; then
    python3 ${scripts}/mixalime/limiter.py --threads $threads create $project "${files[@]}" --no-snp-bad-check --max-cover 10000 
    python3 ${scripts}/mixalime/limiter.py --threads $threads fit $project NB
    python3 ${scripts}/mixalime/limiter.py --threads $threads test $project
    python3 ${scripts}/mixalime/limiter.py --threads $threads combine $project
    python3 ${scripts}/mixalime/limiter.py --threads $threads export all $project $project
    python3 ${scripts}/mixalime/limiter.py --threads $threads plot all $project $project
fi
python3 ${scripts}/mixalime/lessthan500k_sclices.py --home ${home}

# Prepare list of files with different TFs and cell lines for MixALiMe combine 

filtered="${home}/mixalime/file_lists/filtered_indivs.txt"
filtered_names="${home}/mixalime/file_lists/filtered_names.txt"

awk -F/ '{print $NF}' "$filtered" | sort -u > "$filtered_names"
cut -f2 "${home}/clustering/metadata.clustered.tsv" | tail -n +2 | sort -u > "${home}/mixalime/groups/factors.list"
while read -r tf; do
    awk -F'\t' -v tf="$tf" 'NR > 1 && $2 == tf { print $7 }' "${home}/clustering/metadata.clustered.tsv" \
        | sort -u \
        | awk -v home="$home" '{ printf "%s/BEDs/%s.with_bad.bed\t%s.with_bad.bed\n", home, $1, $1 }' \
        | awk -F'\t' 'NR==FNR { ok[$1]=1; next } $2 in ok { print $1 }' "$filtered_names" - \
        > "${home}/mixalime/groups/factors_${tf}.list"
done < "${home}/mixalime/groups/factors.list"

cut -f3 "${home}/clustering/metadata.clustered.tsv" | tail -n +2 | sort -u > "${home}/mixalime/groups/cell.list"
while read -r cell; do
    awk -F'\t' -v cell="$cell" 'NR > 1 && $3 == cell { print $7 }' "${home}/clustering/metadata.clustered.tsv" \
        | sort -u \
        | awk -v home="$home" '{ printf "%s/BEDs/%s.with_bad.bed\t%s.with_bad.bed\n", home, $1, $1 }' \
        | awk -F'\t' 'NR==FNR { ok[$1]=1; next } $2 in ok { print $1 }' "$filtered_names" - \
        > "${home}/mixalime/groups/cell_${cell}.list"
done < "${home}/mixalime/groups/cell.list"

# MixALiMe combine

mkdir -p ${home}/mixalime/multiple_combine

echo [INFO] $(date '+%Y-%m-%d %H:%M:%S') START MIXALIME MULTIPLE_COMBINE FOR TFs > ${home}/logs/status_multiple_combine_factors.txt
while read -r tf; do
    echo [INFO] $(date '+%Y-%m-%d %H:%M:%S') START ${tf} >> ${home}/logs/status_multiple_combine_factors.txt

    projects=$(
        awk -F/ '{print $NF}' ${home}/mixalime/groups/factors_${tf}.list \
        | sed 's/\.with_bad\.bed$//' \
        | while read -r indiv_id; do
            if grep -Fxq "${indiv_id}" ${home}/mixalime/file_lists/halfmillions.filtered.txt; then
                echo ${home}/mixalime/${indiv_id}/${indiv_id}
            elif grep -Fxq "${indiv_id}" ${home}/mixalime/file_lists/not_halfmillions.filtered.txt; then
                echo ${home}/mixalime/less_than_500K_SNPs/less_than_500K_SNPs
            else
                echo [WARNING] $(date '+%Y-%m-%d %H:%M:%S') ${tf} UNKNOWN_PROJECT ${indiv_id} >> ${home}/logs/status_multiple_combine_factors.txt
            fi
        done | sort -u
    )

    if [[ -n "${projects}" ]]; then
        mixalime multiple_combine \
            --n-jobs $threads \
            --subname TF_${tf} \
            --group ${home}/mixalime/groups/factors_${tf}.list \
            ${home}/mixalime/multiple_combine/TF_${tf}/TF_${tf} \
            ${projects}
    else
        echo [WARNING] $(date '+%Y-%m-%d %H:%M:%S') SKIP_EMPTY ${tf} >> ${home}/logs/status_multiple_combine_factors.txt
    fi
done < ${home}/mixalime/groups/factors.list




















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
