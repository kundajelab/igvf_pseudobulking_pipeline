import pandas as pd

# for human transcriptome reference GTF
transcript_rows = []
with open("IGVFFI9573KOZR.gtf", 'r') as f:
    for line in f:
        # if header line
        if line.startswith("#"):
            continue
        # split
        line_split = line.strip().split("\t")
        chro, feature, start, end, strand, info = line_split[0], line_split[2], int(line_split[3]), int(line_split[4]), line_split[6], line_split[8]
        info_split = info.split("; ")
        info_dict = {x.split(" ")[0]: (x.split(" ")[1][1:-1] if x.split(" ")[1].startswith('"') else x.split(" ")[1]) for x in info_split}
        # add to df
        if feature == "transcript":
            new_transcript_dict = dict()
            new_transcript_dict["gene"] = info_dict["gene_name"]
            new_transcript_dict["transcript"] = info_dict["transcript_name"]
            new_transcript_dict["chro"] = chro
            new_transcript_dict["TSS"] = start if strand == "+" else end
            new_transcript_dict["strand"] = strand
            transcript_rows.append(new_transcript_dict)

transcript_df = pd.DataFrame(transcript_rows)
print(transcript_df)
transcript_df.to_csv("./human/tss.tsv", sep="\t", index=False)

# for mouse transcriptome reference GTF
transcript_rows = []
with open("IGVFFI4777RDZK.gtf", 'r') as f:
    for line in f:
        # if header line
        if line.startswith("#"):
            continue
        # split
        line_split = line.strip().split("\t")
        chro, feature, start, end, strand, info = line_split[0], line_split[2], int(line_split[3]), int(line_split[4]), line_split[6], line_split[8]
        info_split = info.split("; ")
        info_dict = {x.split(" ")[0]: (x.split(" ")[1][1:-1] if x.split(" ")[1].startswith('"') else x.split(" ")[1]) for x in info_split}
        # add to df
        if feature == "transcript":
            new_transcript_dict = dict()
            new_transcript_dict["gene"] = info_dict["gene_name"]
            new_transcript_dict["transcript"] = info_dict["transcript_name"]
            new_transcript_dict["chro"] = chro
            new_transcript_dict["TSS"] = start if strand == "+" else end
            new_transcript_dict["strand"] = strand
            transcript_rows.append(new_transcript_dict)

transcript_df = pd.DataFrame(transcript_rows)
print(transcript_df)
transcript_df.to_csv("./mouse/tss.tsv", sep="\t", index=False)