select *, (`crispri.fluoxetine_v2` +`crispri.nortriptyline` +
`crispri.amiodarone_static`+ `crispri.chloroquine` + `crispri.sertraline_v2_correct_gammarhotau` + `crispri.verapamil`)/6 sum from top_crispria_screen_per_gene_v4_more_external
where `crispri.fluoxetine_v2` > 2.5 and `crispri.nortriptyline` > 2.5 and
`crispri.amiodarone_static` > 2.5 and `crispri.chloroquine` > 2.5 and
`crispri.sertraline_v2_correct_gammarhotau` > 2.5 and `crispri.verapamil` > 2.5
order by (`crispri.fluoxetine_v2` +`crispri.nortriptyline` +
`crispri.amiodarone_static`+ `crispri.chloroquine` + `crispri.sertraline_v2_correct_gammarhotau` + `crispri.verapamil`)/6 desc