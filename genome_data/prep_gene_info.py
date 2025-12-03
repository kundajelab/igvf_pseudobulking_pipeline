import pandas as pd
import re

def generate_enhanced_gene_map(gtf_path, output_csv):
    gene_map = {}

    print(f"Parsing {gtf_path}...")
    with open(gtf_path, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if parts[2] != 'gene': continue
            
            attributes = parts[8]
            id_match = re.search(r'gene_id "([^"]+)"', attributes)
            name_match = re.search(r'gene_name "([^"]+)"', attributes)
            
            if id_match:
                gene_id = id_match.group(1)
                # Strip version if your matrix doesn't have them
                # gene_id = gene_id.split('.')[0] 
                
                gene_name = name_match.group(1) if name_match else gene_id
                
                # Check Mito/Ribo status HERE
                is_mito = gene_name.startswith(('MT-', 'mt-'))
                is_ribo = gene_name.startswith(('RPS', 'RPL', 'Rps', 'Rpl'))
                
                gene_map[gene_id] = {
                    'gene_name': gene_name,
                    'mt': is_mito,
                    'ribo': is_ribo
                }

    # Create DataFrame
    df = pd.DataFrame.from_dict(gene_map, orient='index')
    df.index.name = 'gene_id'
    
    # Save
    df.to_csv(output_csv)
    print(f"Saved enriched map to {output_csv}")

# Execute
human_gtf = "IGVFFI9573KOZR.gtf"
mouse_gtf = "IGVFFI4777RDZK.gtf"
generate_enhanced_gene_map(human_gtf, "./human/gene_info.csv")
generate_enhanced_gene_map(mouse_gtf, "./mouse/gene_info.csv")