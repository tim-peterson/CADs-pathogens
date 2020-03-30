

select external_gene_name as gene /*, (avg_nortrip/avg_dmso_hi + avg_verap/avg_dmso_hi + avg_amiod/avg_dmso_hi + avg_sert/avg_dmso_hi + avg_fluox/avg_dmso_hi + avg_chloro/avg_dmso_hi)/6 avg_cads,
(avg_aln/avg_dmso_hi + avg_atorv/avg_dmso_hi + avg_phen/avg_dmso_hi)/3 avg_non_cads, 

(
((avg_nortrip/avg_dmso_hi + avg_verap/avg_dmso_hi + avg_amiod/avg_dmso_hi + avg_sert/avg_dmso_hi + avg_fluox/avg_dmso_hi + avg_chloro/avg_dmso_hi)/6)
/
((avg_aln/avg_dmso_hi + avg_atorv/avg_dmso_hi + avg_phen/avg_dmso_hi)/3)
) 
ratio*/

 from CADs_ge_v3_with_avgs where avg_dmso_hi > 1
and avg_chloro/avg_dmso_hi > 5
and avg_fluox/avg_dmso_hi > 5
and avg_sert/avg_dmso_hi > 5
and avg_amiod/avg_dmso_hi > 5
and avg_verap/avg_dmso_hi > 5
and avg_nortrip/avg_dmso_hi > 5
/*and avg_dmso_hi > 1
and avg_aln/avg_dmso_hi < 2
and avg_atorv/avg_dmso_hi < 2
and avg_phen/avg_dmso_hi < 2
and avg_celecox/avg_dmso_hi < 2*/
order by (
((avg_nortrip/avg_dmso_hi + avg_verap/avg_dmso_hi + avg_amiod/avg_dmso_hi + avg_sert/avg_dmso_hi + avg_fluox/avg_dmso_hi + avg_chloro/avg_dmso_hi)/6)
/
((avg_aln/avg_dmso_hi + avg_atorv/avg_dmso_hi + avg_phen/avg_dmso_hi)/3)
)  desc