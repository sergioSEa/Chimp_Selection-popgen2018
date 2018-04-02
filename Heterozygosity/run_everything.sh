grep Schwein ../Data/pop.info  > separate_files/filter_Schwein ; python scripts/change_info_format.py separate_files/filter_Schwein > separate_files/filter_Schwein_good

grep Verus ../Data/pop.info  > separate_files/filter_Verus ; python scripts/change_info_format.py separate_files/filter_Verus > separate_files/filter_Verus_good

grep Trogl ../Data/pop.info  > separate_files/filter_Trogl ; python scripts/change_info_format.py separate_files/filter_Trogl > separate_files/filter_Trogl_good

plink --bfile ../Data/pruneddata --double-id --keep separate_files/filter_Schwein_good --out separate_files/filter_Schwein --recode
plink --bfile ../Data/pruneddata --double-id --keep separate_files/filter_Verus_good --out separate_files/filter_Verus --recode
plink --bfile ../Data/pruneddata --double-id --keep separate_files/filter_Trogl_good --out separate_files/filter_Trogl --recode

plink --noweb --file separate_files/filter_Schwein --freq --out heterozygosity/Pt_schwein
plink --noweb --file separate_files/filter_Trogl --freq --out heterozygosity/Pt_troglo
plink --noweb --file separate_files/filter_Verus --freq --out heterozygosity/Pt_verus


cat heterozygosity/Pt_schwein.frq |grep -v NA > heterozygosity/Pt_schwein_noNA.frq
cat heterozygosity/Pt_troglo.frq |grep -v NA > heterozygosity/Pt_troglo_noNA.frq
cat heterozygosity/Pt_verus.frq |grep -v NA > heterozygosity/Pt_verus_noNA.frq

awk '$1!=23 ' heterozygosity/Pt_troglo_noNA.frq > heterozygosity/temporal_troglo ; mv heterozygosity/temporal_troglo heterozygosity/Pt_troglo_noNA.frq
awk '$1!=23 ' heterozygosity/Pt_schwein_noNA.frq > heterozygosity/temporal_schwein ; mv heterozygosity/temporal_schwein heterozygosity/Pt_schwein_noNA.frq
awk '$1!=23 ' heterozygosity/Pt_verus_noNA.frq > heterozygosity/temporal_verus ; mv heterozygosity/temporal_verus  heterozygosity/Pt_verus_noNA.frq

Rscript heterozygosity/compute_heterozygosity_v2.R 

mv heterozygosity/*.pdf ./results/


