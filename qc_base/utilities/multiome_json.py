import json
import sys
 
# Opening JSON file
f = open(sys.argv[1])
 
# returns JSON object as
# a dictionary
data = json.load(f)

atac_reads  = data['ATAC Sequencing reads']
atac_status = data['ATAC QC status']
cells       = data['Estimated number of cells']
ds_name     = data['Dataset name']
exons       = data['Reads Mapped to Exons']
gex_frac    = data['Fraction Reads in Cells']
gex_leng    = data['GEX metrics']['Mean trimmed sequence length']
gex_status  = data['GEX QC status']
mito        = data['Mitochondrial reads %']
umi         = data['Median UMI counts per cell']

atac_leng   = data['ATAC metrics']['Mean trimmed sequence length']
frac_map    = data['frac_mapped_confidently']
frac_cut    = data['frac_cut_fragments_in_peaks']

print( ds_name, "\t", cells, "\t", umi, "\t", exons, "\t", gex_frac, "\t", mito, "\t", frac_map, "\t", frac_cut, "\t", gex_leng, "\t", atac_leng, "\t\t\t", atac_status, "\t", gex_status )

# Closing file
f.close()
