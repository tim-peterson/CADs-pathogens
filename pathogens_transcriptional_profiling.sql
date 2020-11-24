

/* no self made jitter using RAND() *, this is actually what went in paper */

select /*all_gene, */ cnt_all, 


CASE
        WHEN (cnt_down IS not NULL or cnt_up is not null)
        THEN 
        	CASE
        		WHEN (cnt_up is null)
        		THEN  (0 - cnt_down)
        		ELSE cnt_up
        	END
        ELSE 0 
       
        END AS cnt_all_0,
        
        
	cnt_up,  (0 - cnt_down) as cnt_down_
 from 

(select final_human as all_gene, cnt as cnt_all from (select *, count(final_human) cnt from (select  final_human, Pathogen as pathogen1, FC_or_Diff, count(final_human) count_per_pathogen from 
(select *, 

/* https://stackoverflow.com/a/11618798/702275*/
    CASE
        WHEN (Gene_symbol IS NULL or Gene_symbol = '')
        THEN
            CASE
            WHEN (`Human gene name` IS NULL or `Human gene name` = '')
                THEN HGNC_symbol
		        ELSE
		            CASE
		            WHEN ( HGNC_symbol IS NULL or  HGNC_symbol = '')
		                THEN `Human gene name`
		                ELSE 'no_gene_human'   
		            END
            
            END
		ELSE Gene_symbol
        END AS final_human
  
 from 
 
(select ensembl_Gene_ID, Gene_symbol, `Human gene name`, HGNC_symbol, FC_or_Diff, Pathogen from pathogens_transcriptomics_w_pathogen_v2 a
left join (select * from morpheome.ensembl_mouse_w_human where `Human gene name` is not null group by `Gene stable ID`)  b

on a.ensembl_Gene_ID = b.`Gene stable ID`
left join (select * from morpheome.aliases_ensembl where `HGNC_symbol` is not null group by `Gene stable ID`)  c

on a.ensembl_Gene_ID = c.`Gene stable ID`) d) e where final_human != "" /*and FC_or_Diff < 0  /*and FC_or_Diff > 1  /*and lower(final_human) = "chrna3"*/ group by final_human, Pathogen order by count(final_human) desc) f

left join (select Pathogen, count(Pathogen) meta_count from pathogens_transcriptomics_metadata_v2 where used = "used" group by Pathogen) g
on
f.Pathogen1=g.Pathogen /*where (count_per_pathogen > 1 and meta_count > 1) or (count_per_pathogen = 1 and meta_count = 1)*/ group by final_human order by cnt desc ) h) all_genes


left join 

(select final_human as up_gene, cnt as cnt_up from (select *, count(final_human) cnt from (select  final_human, Pathogen as pathogen1, FC_or_Diff, count(final_human) count_per_pathogen from 
(select *, 

/* https://stackoverflow.com/a/11618798/702275*/
    CASE
        WHEN (Gene_symbol IS NULL or Gene_symbol = '')
        THEN
            CASE
            WHEN (`Human gene name` IS NULL or `Human gene name` = '')
                THEN HGNC_symbol
		        ELSE
		            CASE
		            WHEN ( HGNC_symbol IS NULL or  HGNC_symbol = '')
		                THEN `Human gene name`
		                ELSE 'no_gene_human'   
		            END
            
            END
		ELSE Gene_symbol
        END AS final_human
  
 from 
 
(select ensembl_Gene_ID, Gene_symbol, `Human gene name`, HGNC_symbol, FC_or_Diff, Pathogen from pathogens_transcriptomics_w_pathogen_v2 a
left join (select * from morpheome.ensembl_mouse_w_human where `Human gene name` is not null group by `Gene stable ID`)  b

on a.ensembl_Gene_ID = b.`Gene stable ID`
left join (select * from morpheome.aliases_ensembl where `HGNC_symbol` is not null group by `Gene stable ID`)  c

on a.ensembl_Gene_ID = c.`Gene stable ID`) d) e where final_human != "" and FC_or_Diff > 1 /*and FC_or_Diff < 0  /*and FC_or_Diff > 1  /*and lower(final_human) = "chrna3"*/ group by final_human, Pathogen order by count(final_human) desc) f

left join (select Pathogen, count(Pathogen) meta_count from pathogens_transcriptomics_metadata_v2 where used = "used" group by Pathogen) g
on
f.Pathogen1=g.Pathogen where (count_per_pathogen > 1 and meta_count > 1) or (count_per_pathogen = 1 and meta_count = 1) group by final_human order by cnt desc ) h where cnt > 2) upreg
on all_genes.all_gene=upreg.up_gene


left join 

(select final_human as down_gene, cnt as cnt_down from (select *, count(final_human) cnt from (select  final_human, Pathogen as pathogen1, FC_or_Diff, count(final_human) count_per_pathogen from 
(select *, 

/* https://stackoverflow.com/a/11618798/702275*/
    CASE
        WHEN (Gene_symbol IS NULL or Gene_symbol = '')
        THEN
            CASE
            WHEN (`Human gene name` IS NULL or `Human gene name` = '')
                THEN HGNC_symbol
		        ELSE
		            CASE
		            WHEN ( HGNC_symbol IS NULL or  HGNC_symbol = '')
		                THEN `Human gene name`
		                ELSE 'no_gene_human'   
		            END
            
            END
		ELSE Gene_symbol
        END AS final_human
  
 from 
 
(select ensembl_Gene_ID, Gene_symbol, `Human gene name`, HGNC_symbol, FC_or_Diff, Pathogen from pathogens_transcriptomics_w_pathogen_v2 a
left join (select * from morpheome.ensembl_mouse_w_human where `Human gene name` is not null group by `Gene stable ID`)  b

on a.ensembl_Gene_ID = b.`Gene stable ID`
left join (select * from morpheome.aliases_ensembl where `HGNC_symbol` is not null group by `Gene stable ID`)  c

on a.ensembl_Gene_ID = c.`Gene stable ID`) d) e where final_human != "" and FC_or_Diff < 0 /*and FC_or_Diff < 0  /*and FC_or_Diff > 1  /*and lower(final_human) = "chrna3"*/ group by final_human, Pathogen order by count(final_human) desc) f

left join (select Pathogen, count(Pathogen) meta_count from pathogens_transcriptomics_metadata_v2 where used = "used" group by Pathogen) g
on
f.Pathogen1=g.Pathogen where (count_per_pathogen > 1 and meta_count > 1) or (count_per_pathogen = 1 and meta_count = 1) group by final_human order by cnt desc ) h where cnt > 1) downreg
on all_genes.all_gene=downreg.down_gene























