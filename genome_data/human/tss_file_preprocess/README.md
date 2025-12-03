# Instructions used for generating tss.tsv file.

1) wget https://api.data.igvf.org/reference-files/IGVFFI9573KOZR/@@download/IGVFFI9573KOZR.gtf.gz (from https://data.igvf.org/reference-files/IGVFFI9573KOZR/)
2) run `gunzip IGVFFI9573KOZR.gtf.gz`
3) run `python prep\_tss\_data.py`
