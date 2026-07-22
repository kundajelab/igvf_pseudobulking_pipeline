table peaks
"MACS3 Narrow Peak file"
(
string  chrom;		"Reference sequence chromosome or scaffold"
uint    chromStart;	"Start position of feature on chromosome"
uint    chromEnd;	"End position of feature on chromosome"
string  name;		"Name of gene"
uint    score;		"Score"
char[1] strand;	"+ or - for strand"
float   signalValue;	"Overal fold-enrichment for the region"
float   pValue;	"-log10 measure of significance"
float  	qValue;	"Meausure of signifcance using false-discovery"
uint    peak;	"peak summit location, relative to chromStart"
)
