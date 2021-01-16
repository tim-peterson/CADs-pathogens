/* CAD selection from Biogrid */
select * from biogrid_crispr_screens a
join DrugBank.drug b
on a.CONDITION_NAME=b.name
join (select * from DrugBank.drug_calculated_properties where kind = 'logP' and value > 2) c
on b.`primary_key`=c.`parent_key`

join (select * from DrugBank.drug_calculated_properties where kind = 'pKa (strongest basic)' and value > 6.5 ) d
on b.`primary_key`=d.`parent_key`
where FULL_SIZE > 10000

/* non-CAD selection from Biogrid */
select * from biogrid_crispr_screens a
join DrugBank.drug b
on a.CONDITION_NAME=b.name
join (select * from DrugBank.drug_calculated_properties where (kind = 'logP' and value < 2) or (kind = 'pKa (strongest basic)' and value < 6.5)) c
on b.`primary_key`=c.`parent_key`
where FULL_SIZE > 10000
order by name desc
