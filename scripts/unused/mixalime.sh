TF=$1
model=$2
home='/home/subpolare/adastra-v7'

mixalime create ${home}/mixalime/${TF}_${model} ${home}/VCFs/${TF}_*.vcf.gz #--bad-maps ${home}/BABACHI/${TF}.badmap.bed
mixalime fit ${home}/mixalime/${TF}_${model} $model
mixalime test ${home}/mixalime/${TF}_${model}
mixalime combine ${home}/mixalime/${TF}_${model}

# mixalime export all ${home}/mixalime/${TF}_${model} ${home}/mixalime/results_${TF}_${model}
# mixalime plot all ${home}/mixalime/${TF}_${model} ${home}/mixalime/results_${TF}_${model}

mixalime export all ${home}/mixalime/${TF}_${model} ${home}/mixalime/results_${TF}_${model}_without_BAD
mixalime plot all ${home}/mixalime/${TF}_${model} ${home}/mixalime/results_${TF}_${model}_without_BAD
