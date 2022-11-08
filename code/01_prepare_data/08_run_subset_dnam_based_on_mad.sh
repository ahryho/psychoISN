# LD pruning on MAFed genotypes

cd /binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/input/dnam

# MAD filtering

treatment="dex"
mad_thr=75
in_gds_fn="gds/methyl_beta_mtrx_corrected_for_cov_"$treatment".gds"
out_gds_fn="mad_filtered/gds/methyl_beta_mtrx_corrected_for_cov_mad"$mad_thr"_filtered_"$treatment".gds"
out_gds_fn_prefix="mad_filtered/gds/"

Rscript --vanilla ~/kul/dex-stim-human-array-isns/code/01_prepare_data/08_subset_dnam_based_on_mad.R \
$treatment $mad_thr $in_gds_fn $out_gds_fn \
> mad_filtered/logs/subset_dnam_${treatment}_based_on_mad_${mad_thr}.log

# Split by chromosome

Rscript --vanilla ~/kul/dex-stim-human-array-isns/code/01_prepare_data/06b_split_dnam_by_chrom.R \
$treatment $out_gds_fn $out_gds_fn_prefix \
> mad_filtered/logs/split_dnam_${treatment}_mad_${mad_thr}_gds_by_chrom.log
