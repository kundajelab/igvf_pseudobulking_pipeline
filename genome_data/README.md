# Instructions for generating chr\_sizes.tsv
1) wget https://api.data.igvf.org/reference-files/IGVFFI0653VCGH/@@download/IGVFFI0653VCGH.fasta.gz (GRCh38 genome reference from https://data.igvf.org/reference-files/IGVFFI0653VCGH/)
2) wget https://api.data.igvf.org/reference-files/IGVFFI9282QLXO/@@download/IGVFFI9282QLXO.fasta.gz (GRCm39 genome reference from https://data.igvf.org/reference-files/IGVFFI9282QLXO/)
3) gunzip IGVFFI0653VCGH.fasta.gz
4) gunzip IGVFFI9282QLXO.fasta.gz
5) samtools faidx IGVFFI0653VCGH.fasta
6) samtools faidx IGVFFI9282QLXO.fasta
7) cut -f 1,2 IGVFFI0653VCGH.fasta.fai > ./human/chr\_sizes.tsv
8) cut -f 1,2 IGVFFI9282QLXO.fasta.fai > ./mouse/chr\_sizes.tsv

# Instructions for generating gene\_info.csv
1) wget https://api.data.igvf.org/reference-files/IGVFFI9573KOZR/@@download/IGVFFI9573KOZR.gtf.gz (GRCh38 GENECODE 43 transcriptome reference from https://data.igvf.org/reference-files/IGVFFI9573KOZR/)
2) wget https://api.data.igvf.org/reference-files/IGVFFI4777RDZK/@@download/IGVFFI4777RDZK.gtf.gz (GRCm39 GENECODE M36 transcriptome reference from https://data.igvf.org/reference-files/IGVFFI4777RDZK/)
3) gunzip IGVFFI9573KOZR.gtf.gz
4) gunzip IGVFFI4777RDZK.gtf.gz
5) python prep\_gene\_info.py

# Instructions for generating blacklist.bed
1) wget https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz (GRCh38 blacklist from https://www.encodeproject.org/files/ENCFF356LFX/)
2) download mm39.excluderanges.bed from https://dozmorovlab.github.io/excluderanges/ (GRCm39 blacklist through the package (excluderanges) or google drive related to the package(https://drive.google.com/drive/folders/1sF9m8Y3eZouTZ3IEEywjs2kfHOWFBSJT))
3) gunzip ENCFF356LFX.bed.gz
4) mv ENCFF356LFX.bed ./human/blacklist.bed
5) cut -f 1,2,3 mm39.excluderanges.bed > ./mouse/blacklist.bed

# Instructions for generating tss.tsv file.
1) wget https://api.data.igvf.org/reference-files/IGVFFI9573KOZR/@@download/IGVFFI9573KOZR.gtf.gz (GRCh38 GENECODE 43 transcriptome reference from https://data.igvf.org/reference-files/IGVFFI9573KOZR/)
2) wget https://api.data.igvf.org/reference-files/IGVFFI4777RDZK/@@download/IGVFFI4777RDZK.gtf.gz (GRCm39 GENECODE M36 transcriptome reference from https://data.igvf.org/reference-files/IGVFFI4777RDZK/)
3) gunzip IGVFFI9573KOZR.gtf.gz
4) gunzip IGVFFI4777RDZK.gtf.gz
5) python prep\_tss\_data.py
