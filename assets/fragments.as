table fragments
"ATAC-seq fragments file"
(
string  chrom;		"Reference sequence chromosome or scaffold"
uint    chromStart;	"Start position of feature on chromosome"
uint    chromEnd;	"End position of feature on chromosome"
string  barcode_sample;		"Barcode and sample accession joined by underscore."
uint    num_reads;	"Number of observed copies of this fragment."
)
