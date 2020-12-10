

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


/* CADs 3X downreguled genes fatty acid elongation */
select gene from (select external_gene_name as gene, (avg_nortrip/avg_dmso_hi + avg_verap/avg_dmso_hi + avg_amiod/avg_dmso_hi + avg_sert/avg_dmso_hi + avg_fluox/avg_dmso_hi + avg_chloro/avg_dmso_hi)/6 avg_cads,
(avg_aln/avg_dmso_hi + avg_atorv/avg_dmso_hi + avg_phen/avg_dmso_hi)/3 avg_non_cads

 from CADs_ge_v3_with_avgs
/*where avg_chloro/avg_dmso_hi > 10
and avg_fluox/avg_dmso_hi > 10
and avg_sert/avg_dmso_hi > 10
and avg_amiod/avg_dmso_hi > 10
and avg_verap/avg_dmso_hi > 10
and avg_nortrip/avg_dmso_hi > 10
and avg_dmso_hi > 1
and avg_aln/avg_dmso_hi < 2
and avg_atorv/avg_dmso_hi < 2
and avg_phen/avg_dmso_hi < 2
and avg_celecox/avg_dmso_hi < 2*/
where avg_dmso_hi > 0.1
order by (avg_nortrip/avg_dmso_hi + avg_verap/avg_dmso_hi + avg_amiod/avg_dmso_hi + avg_sert/avg_dmso_hi + avg_fluox/avg_dmso_hi + avg_chloro/avg_dmso_hi)/6 asc) a where  avg_cads > 0 and avg_cads/avg_non_cads < 0.33


/* CADs all genes for paper scatter plot */
select external_gene_name as gene, avg_cads, avg_non_cads, avg_cads/avg_non_cads ratio from (select external_gene_name, (avg_nortrip/avg_dmso_hi + avg_verap/avg_dmso_hi + avg_amiod/avg_dmso_hi + avg_sert/avg_dmso_hi + avg_fluox/avg_dmso_hi + avg_chloro/avg_dmso_hi)/6 avg_cads,
(avg_aln/avg_dmso_hi + avg_atorv/avg_dmso_hi + avg_phen/avg_dmso_hi)/3 avg_non_cads

 from CADs_ge_v3_with_avgs
/*where avg_chloro/avg_dmso_hi > 10
and avg_fluox/avg_dmso_hi > 10
and avg_sert/avg_dmso_hi > 10
and avg_amiod/avg_dmso_hi > 10
and avg_verap/avg_dmso_hi > 10
and avg_nortrip/avg_dmso_hi > 10
and avg_dmso_hi > 1
and avg_aln/avg_dmso_hi < 2
and avg_atorv/avg_dmso_hi < 2
and avg_phen/avg_dmso_hi < 2
and avg_celecox/avg_dmso_hi < 2*/
where avg_dmso_hi > 0.1
order by (avg_nortrip/avg_dmso_hi + avg_verap/avg_dmso_hi + avg_amiod/avg_dmso_hi + avg_sert/avg_dmso_hi + avg_fluox/avg_dmso_hi + avg_chloro/avg_dmso_hi)/6 desc) a /* where avg_non_cads > 0.1 and avg_cads > 0.1*/ order by ratio desc






/* p53 analysis */
select * from p53_targets a
left join `CADs_ge_v3_gt0.1_with_ratio` b
on a.gene=b.gene where ratio is not null and ratio > 3


select * from `CADs_ge_v3_gt0.1_with_ratio` where ratio > 3
















