select *, (`crispri.fluoxetine_v2` +`crispri.nortriptyline` +
`crispri.amiodarone_static`+ `crispri.chloroquine` + `crispri.sertraline_v2_correct_gammarhotau` + `crispri.verapamil`)/6 sum from  top_crispria_screen_per_gene_v13_correct_nonCADs     /*top_crispria_screen_per_gene_v4_more_external*/
where abs(`crispri.fluoxetine_v2`) > 2.5 and abs(`crispri.nortriptyline`) > 2.5 and
abs(`crispri.amiodarone_static`) > 2.5 and abs(`crispri.chloroquine`) > 2.5 and
abs(`crispri.sertraline_v2_correct_gammarhotau`) > 2.5 and abs(`crispri.verapamil`) > 2.5
order by (`crispri.fluoxetine_v2` +`crispri.nortriptyline` +
`crispri.amiodarone_static`+ `crispri.chloroquine` + `crispri.sertraline_v2_correct_gammarhotau` + `crispri.verapamil`)/6 asc