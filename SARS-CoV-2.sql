/* get all CADs from DrugBank */
select * from (select kind,value, parent_key from drug_experimental_properties where (kind = 'pKa' and value > 8)  ) a
join (select kind,value, parent_key from drug_experimental_properties where (kind = 'logP' and value > 2) group by parent_key ) b
on a.parent_key=b.parent_key
join drug c
on a.parent_key=c.primary_key


select * from (select pka, val, logP, value, c.name from (select kind as pka ,value as val, parent_key from drug_calculated_properties where (kind = 'pKa (strongest basic)' and value > 8)  ) a
join (select kind as logP, value, parent_key from drug_calculated_properties where (kind = 'logP' and value > 2) group by parent_key ) b
on a.parent_key=b.parent_key
join drug c
on a.parent_key=c.primary_key ) d
right join SARS_COV2.`SARS-CoV-2-GSE147507-gt1-results_drugs_only` e
on d.name=e.name order by summary desc


/* gene expression */

select Gene from PMID_32747830_SuppTable_3_4_DEGandPathways_preprint where SampleGroup = "all_samples" and Comparison = "Positive_vs_Negative" and log2FoldChange > 1

select * from SARS_COV2.`PMID_32747830_2xup_pos_vs_neg_all_samples (575 genes)` e
left join (select * from (select kind as pka ,value as val, parent_key as pk_a from drug_calculated_properties where (kind = 'pKa (strongest basic)' and value > 8)  ) a
join (select kind as logP, value, parent_key from drug_calculated_properties where (kind = 'logP' and value > 2) group by parent_key ) b
on a.pk_a=b.parent_key
join drug c
on a.pk_a=c.primary_key ) d
on e.Term like concat("", d.name, "%") order by `P-value` asc 


select Gene_name from `PMID_32416070_patients_DEGs-mmc4` where log2FoldChange > 2.5


select SUBSTRING_INDEX(aa.Term_e, ' ', 1) AS name_, val as pka_, value as logP_, pval_e, pval_x, odds_e, odds_x from  (select val, value, e.Term as Term_e, pka,logP, e.`P-value` pval_e, e.`Odds Ratio` odds_e, x.`P-value` pval_x, x.`Odds Ratio` odds_x from SARS_COV2.`PMID_32747830_2xup_pos_vs_neg_all_samples (575 genes)` e
join SARS_COV2.`PMID_32416070_patients_DEGs-mmc4_log(2.5)~5.65Xup (495 genes)` x
on e.Term=x.Term
left join (select * from (select kind as pka ,value as val, parent_key as pk_a from drug_calculated_properties where (kind = 'pKa (strongest basic)' and value > 8 and value < 12)  ) a
join (select kind as logP, value, parent_key from drug_calculated_properties where (kind = 'logP' and value > 3) group by parent_key ) b
on a.pk_a=b.parent_key
join drug c
on a.pk_a=c.primary_key ) d
on e.Term like concat("", d.name, "%") order by x.`P-value` asc ) aa group by name_ order by pval_e asc

/* count CADs out of 100 */
select SUM(case when pka_ is not null then 1 else 0 end) as countKindOne from (select SUBSTRING_INDEX(aa.Term_e, ' ', 1) AS name_, val as pka_, value as logP_, pval_e, pval_x, odds_e, odds_x from  (select val, value, e.Term as Term_e, pka,logP, e.`P-value` pval_e, e.`Odds Ratio` odds_e, x.`P-value` pval_x, x.`Odds Ratio` odds_x from SARS_COV2.`PMID_32747830_2xup_pos_vs_neg_all_samples (575 genes)` e
join SARS_COV2.`PMID_32416070_patients_DEGs-mmc4_log(2.5)~5.65Xup (495 genes)` x
on e.Term=x.Term
left join (select * from (select kind as pka ,value as val, parent_key as pk_a from drug_calculated_properties where (kind = 'pKa (strongest basic)' and value > 8 and value < 12)  ) a
join (select kind as logP, value, parent_key from drug_calculated_properties where (kind = 'logP' and value > 3) group by parent_key ) b
on a.pk_a=b.parent_key
join drug c
on a.pk_a=c.primary_key ) d
on e.Term like concat("", d.name, "%") order by x.`P-value` asc ) aa group by name_ order by pval_e asc limit 100) z 

select SUM(case when pka_ is not null then 1 else 0 end) as countKindOne from (select SUBSTRING_INDEX(aa.Term_e, ' ', 1) AS name_, val as pka_, value as logP_, pval_e, pval_x, odds_e, odds_x from  (select val, value, e.Term as Term_e, pka,logP, e.`Adjusted P-value` pval_e, e.`Odds Ratio` odds_e, x.`Adjusted P-value` pval_x, x.`Odds Ratio` odds_x from SARS_COV2.`PMID_32747830_2xup_pos_vs_neg_all_samples (575 genes)` e
join SARS_COV2.`PMID_32416070_patients_DEGs-mmc4_log(2.5)~5.65Xup (495 genes)` x
on e.Term=x.Term
left join (select * from (select kind as pka ,value as val, parent_key as pk_a from drug_calculated_properties where (kind = 'pKa (strongest basic)' and value > 7)  ) a
join (select kind as logP, value, parent_key from drug_calculated_properties where (kind = 'logP' and value > 0) group by parent_key ) b
on a.pk_a=b.parent_key
join drug c
on a.pk_a=c.primary_key ) d
on e.Term like concat("", d.name, "%") order by x.`P-value` asc ) aa group by name_ order by pval_e asc /*limit 100*/) z 




select SUM(case when pka_ is not null then 1 else 0 end) as countKindOne from (select SUBSTRING_INDEX(aa.Term_e, ' ', 1) AS name_, val as pka_, value as logP_, pval_e, odds_e from  (select val, value, e.Term as Term_e, pka,logP, e.`P-value` pval_e, e.`Odds Ratio` odds_e from ssri.`pathogens_tomics_gt2_pathogen_gt1_ (302 genes)` e

left join (select * from (select kind as pka ,value as val, parent_key as pk_a from drug_calculated_properties where (kind = 'pKa (strongest basic)' and value > 7.5 and value < 13)  ) a
join (select kind as logP, value, parent_key from drug_calculated_properties where (kind = 'logP' and value > 1) group by parent_key ) b
on a.pk_a=b.parent_key
join drug c
on a.pk_a=c.primary_key ) d
on e.Term like concat("", d.name, "%") order by e.`P-value` asc ) aa group by name_ order by pval_e asc  limit 100) z 


select SUM(case when pka_ is not null then 1 else 0 end) as countKindOne from (select SUBSTRING_INDEX(aa.Term_e, ' ', 1) AS name_, val as pka_, value as logP_, pval_e, odds_e from  (select val, value, e.Term as Term_e, pka,logP, e.`Adjusted P-value` pval_e, e.`Odds Ratio` odds_e from ssri.`pathogens_tomics_gt2_pathogen_gt1_ (302 genes)` e

left join (select * from (select kind as pka ,value as val, parent_key as pk_a from drug_calculated_properties where (kind = 'pKa (strongest basic)' and value > 8 and value < 12)  ) a
join (select kind as logP, value, parent_key from drug_calculated_properties where (kind = 'logP' and value > 3) group by parent_key ) b
on a.pk_a=b.parent_key
join drug c
on a.pk_a=c.primary_key ) d
on e.Term like concat("", d.name, "%") order by e.`P-value` asc ) aa /*where pval_e < 0.00005*/ group by name_ order by odds_e desc) z 



select SUM(case when pka_ is not null then 1 else 0 end) as countKindOne from (select SUBSTRING_INDEX(aa.Term_e, ' ', 1) AS name_, val as pka_, value as logP_, pval_e, odds_e from  (select val, value, e.Term as Term_e, pka,logP, e.`Adjusted P-value` pval_e, e.`Odds Ratio` odds_e from ssri.`pathogens_v2_upregulated_only_used_only_gt2 (475 genes)` e

left join (select * from (select kind as pka ,value as val, parent_key as pk_a from drug_calculated_properties where (kind = 'pKa (strongest basic)' and  value > 7.5 and value < 13)  ) a
join (select kind as logP, value, parent_key from drug_calculated_properties where (kind = 'logP' and value > 1) group by parent_key ) b
on a.pk_a=b.parent_key
join drug c
on a.pk_a=c.primary_key ) d
on e.Term like concat("", d.name, "%") order by e.`P-value` asc ) aa /*where pval_e < 0.05*/ group by name_ order by pval_e asc) z 



select SUM(case when pka_ is not null then 1 else 0 end) as countKindOne from (select SUBSTRING_INDEX(aa.Term_e, ' ', 1) AS name_, val as pka_, value as logP_, pval_e, pval_x, odds_e, odds_x from  (select val, value, e.Term as Term_e, pka,logP, e.`Adjusted P-value` pval_e, e.`Odds Ratio` odds_e, x.`Adjusted P-value` pval_x, x.`Odds Ratio` odds_x from SARS_COV2.`PMID_32747830_2xup_pos_vs_neg_all_samples (575 genes)` e
join SARS_COV2.`PMID_32416070_patients_DEGs-mmc4_log(2.5)~5.65Xup (495 genes)` x
on e.Term=x.Term
left join (select * from (select kind as pka ,value as val, parent_key as pk_a from drug_calculated_properties where (kind = 'pKa (strongest basic)' and value > 8 and value < 12)  ) a
join (select kind as logP, value, parent_key from drug_calculated_properties where (kind = 'logP' and value > 3) group by parent_key ) b
on a.pk_a=b.parent_key
join drug c
on a.pk_a=c.primary_key ) d
on e.Term like concat("", d.name, "%") order by x.`P-value` asc ) aa group by name_ order by pval_e asc limit 100) z 