


# biofiles -- Dealing with GenBank files in R

### Installation

This package is under development and currently ony available via github.
It depends on functions provided in the [rmisc](https://github.com/gschofl/rmisc)
package. The package [Rentrez](https://github.com/gschofl/rentrez) is helpful
to interface with the NCBI EUtilities.



```r
install_github("rmisc", "gschofl")
install_github("rentrez", "gschofl")
install_github("biofiles", "gschofl")
```





### Basic functionality

Let's download a small bacterial genome, _Chlamydophila psittaci_ DC15
(GenBank acc. CP002806; GI 334693792) from NCBI.



```r
require(Rentrez)
gb_file <- efetch("CP002806", "nuccore", rettype = "gbwithparts", retmode = "text")
```




Next, we parse the `efetch` object into a `gbRecord` instance.



```r
gb_record <- gbRecord(gb_file)
gb_record

 ##  'gbRecord' instance with 2040 features
 ##  LOCUS       CP002806 1172182 bp DNA circular BCT 2011-08-16
 ##  DEFINITION  Chlamydophila psittaci 02DC15, complete genome.
 ##  ACCESSION   CP002806
 ##  VERSION     CP002806.1 GI:334693792
 ##  DBLINK      Project: PRJNA66481
 ##  KEYWORDS    .
 ##  SOURCE      Chlamydophila psittaci 02DC15
 ##    ORGANISM  Chlamydophila psittaci 02DC15
 ##              Bacteria; Chlamydiae; Chlamydiales; Chlamydiaceae; Chlamydia
 ##               /Chlamydophila group; Chlamydia.
 ##  REFERENCE   Not implemented yet
 ##  COMMENT     Source DNA available from Dr. Konrad Sachse (Konrad.Sachse@f
 ##               li.bund.de) at the Institute for Molecular Pathogenesis, Jena, Germany.
 ##  ORIGIN      CAAAGTTTTAAACATGTTAACCACGTTGTTTTCCCTCTATAAACCGAGTTTTAAACAGTTT
 ##              ...
 ##              TGATCTTGTGTTCACAAACATGCAAAATCTCGATTTATTTACAAACTGAACAATGTTTTAA
 ##  CONTIG
```




The `summary` function provides an overview over the object:



```r
summary(gb_record)

 ##  [[CP002806]]
 ##    1172182 bp: Chlamydophila psittaci 02DC15, complete genome.
 ##    N    Key    Location                     Product
 ##    1    source 1..1172182                   NA
 ##    2    gene   149..1144                    NA
 ##    3    CDS    149..1144                    delta-aminolevulinic acid d...
 ##    4    gene   complement(1160..2578)       NA
 ##    5    CDS    complement(1160..2578)       Na(+)-translocating NADH:ub...
 ##    6    gene   complement(2601..3035)       NA
 ##    7    CDS    complement(2601..3035)       conserved hypothetical protein
 ##    ...  ...    ...                          ...
 ##    2034 CDS    1169125..1170006             geranylgeranyl pyrophosphat...
 ##    2035 gene   complement(1170021..1170094) NA
 ##    2036 tRNA   complement(1170021..1170094) tRNA-Pro
 ##    2037 gene   complement(1170106..1170360) NA
 ##    2038 CDS    complement(1170106..1170360) conserved hypothetical protein
 ##    2039 gene   complement(1170443..1172035) NA
 ##    2040 CDS    complement(1170443..1172035) conserved hypothetical protein
```




Various getter methods provide access to the data contained in a GenBank record;
for instance:



```r
getAccession(gb_record)

 ##  [1] "CP002806"
```





```r
getGeneID(gb_record)

 ##  [1] "334693792"
```





```r
getDefinition(gb_record)

 ##  [1] "Chlamydophila psittaci 02DC15, complete genome."
```





```r
getOrganism(gb_record)

 ##  [1] "Chlamydophila psittaci 02DC15"
```





```r
getSequence(gb_record)

 ##    A DNAStringSet instance of length 1
 ##        width seq                                        names               
 ##  [1] 1172182 CAAAGTTTTAAACATGTTAA...AAACTGAACAATGTTTTAA CP002806
```





### Genbank Features

The important part of GenBank record will generally be the list of annotions
or features.

We can access the `gbFeatureList` of a `gbRecord` using `getFeatures()`:



```r
gb_features <- getFeatures(gb_record)
gb_features

 ##  'gbFeatureList' with 2040 features:
 ##  Feature:         Location/Qualifiers:
 ##   source          1..1172182
 ##                   /organism="Chlamydophila psittaci 02DC15"
 ##                   /mol_type="genomic DNA"
 ##                   /strain="02DC15"
 ##                   /db_xref="taxon:1027845"
 ##  ...
 ##  Feature:         Location/Qualifiers:
 ##   CDS             complement(1170443..1172035)
 ##                   /locus_tag="CPS0B_1080"
 ##                   /note="[F] COG1351 Predicted alternative thymidylate synthase"
 ##                   /codon_start="1"
 ##                   /transl_table="11"
 ##                   /product="conserved hypothetical protein"
 ##                   /protein_id="AEG87987.1"
 ##                   /db_xref="GI:334694770"
 ##                   /translation="MLSRDDEFSSEQRKSLSHFVTNLETNIFALKNLPEVVKGA
 ##                   LFSKYSRSTLGLRSLLLKEFLEGEGGDFLDSSGVDFEVGIHKAADFYRRVLDGF
 ##                   GDDSIGELGGAHLAIESVSMLAAKILEDARIGGSPLEKSSRYVYFDQKVKGEYL
 ##                   YYRDPILMTSAFKDVFLGTCDFLFDTYADLIPKVRTYFEKIYPKESEVSQSAYT
 ##                   ISLRAKVLDCLRGLLPAATLTNLGFFGNGRFWQTLLHKIQGHNLTEIRQIGESS
 ##                   LTELMKIIPSFVSRAESHHHHHQAMLSYRQTLREQLTSLAEKFSGGSHPSKQTG
 ##                   VRLVYGDPEGIYKVAAGFLFPYSEHTYEELINICKSMPREDLIRVLEAGSSSRE
 ##                   NRRHKSPRGLECLEFGFDITADFGAYRDLQRHRILTQERQLLTTNLGYHIPEQL
 ##                   LDTPMEKDFREAMEKAEEAYNQISLEFPEEAQYVVPLAYNIRWFFHINGRALQW
 ##                   LCELRSQVQGHENYRKIAIDMVKEVIRFDPVYESFFKFVDYSECDLGRIKQESRKRSF"
 ##  Seqinfo:
 ##  CP002806  1172182 DNA  Chlamydophila psittaci 02DC15, complete genome.
```




We can narrow down the list of features of interest by using `select()`.
For instance, we can select for all coding sequences (CDS):



```r
cds <- select(gb_features, key = "CDS")
summary(cds)

 ##  N    Key Location                     Product
 ##  3    CDS 149..1144                    delta-aminolevulinic acid dehydr...
 ##  5    CDS complement(1160..2578)       Na(+)-translocating NADH:ubiquin...
 ##  7    CDS complement(2601..3035)       conserved hypothetical protein
 ##  9    CDS 3148..5301                   transcription elongation factor ...
 ##  13   CDS 5452..6645                   aromatic amino acid aminotransfe...
 ##  15   CDS 6889..7956                   putative cell shape-determining ...
 ##  17   CDS complement(7916..11053)      exodeoxyribonuclease V subunit beta
 ##  19   CDS complement(11040..14108)     putative exodeoxyribonuclease V ...
 ##  ...  ... ...                          ...
 ##  2024 CDS 1163308..1164144             conserved hypothetical protein
 ##  2026 CDS complement(1164154..1165668) exodeoxyribonuclease V subunit a...
 ##  2028 CDS complement(1165683..1167527) conserved hypothetical protein
 ##  2030 CDS 1167722..1168435             transcriptional regulatory protein
 ##  2032 CDS 1168495..1169109             putative glucosamine-1-phosphate...
 ##  2034 CDS 1169125..1170006             geranylgeranyl pyrophosphate syn...
 ##  2038 CDS complement(1170106..1170360) conserved hypothetical protein
 ##  2040 CDS complement(1170443..1172035) conserved hypothetical protein
```




or we can select for all elongation factors,



```r
elong <- select(gb_features, key = "CDS", product = "elongation factor")
summary(elong)

 ##  N    Key Location                   Product
 ##  9    CDS 3148..5301                 transcription elongation factor Gr...
 ##  103  CDS complement(56119..56967)   translation elongation factor Ts
 ##  408  CDS 206213..208297             translation elongation factor G
 ##  950  CDS complement(526454..527758) transcription elongation factor
 ##  1126 CDS 623470..624027             elongation factor P
 ##  1406 CDS complement(805145..806329) translation elongation factor Tu
 ##  1750 CDS complement(997270..997842) translation elongation factor P
```




now let's extract the sequence for all elongation factors, and using the tools
from the `Biostrings` packages, translate them into protein sequences



```r
require(Biostrings)
dna <- getSequence(elong)
dna

 ##    A DNAStringSet instance of length 7
 ##      width seq                                          names               
 ##  [1]  2154 GTGGACTATCTAGAAAAGTTG...AATCTATTTGGGATGCATAA lcl|CDS.9|gb|CP00...
 ##  [2]   849 TTAAGCTCCTATTTTCCATAA...TCCATAGAAAAGTTGCTCAT lcl|CDS.103|gb|CP...
 ##  [3]  2085 ATGAGTGATCAAGAATTCGAT...AAGAGATTGTTAAGAAGTAA lcl|CDS.408|gb|CP...
 ##  [4]  1305 TTAATCTTCAATTTTAGGTTT...GCTACAAGATCTTTATTCAT lcl|CDS.950|gb|CP...
 ##  [5]   558 ATGGTATTAAGTAGCCAGCTC...AATATATTCAACGCGTCTAA lcl|CDS.1126|gb|C...
 ##  [6]  1185 TTAGGCAATAATTTTTGAAAT...TGAAAAGTTTCTTTTGACAT lcl|CDS.1406|gb|C...
 ##  [7]   573 TTACTTGGAAACACGCGATTC...CTAGTGCTTACACGAACCAT lcl|CDS.1750|gb|C...
```






```r
str <- strand(elong)
dna_revcomp <- c(reverseComplement(dna[str == -1]), dna[str == 1])
aa <- translate(dna_revcomp)
names(aa) <- names(dna_revcomp)
aa

 ##    A AAStringSet instance of length 7
 ##      width seq                                          names               
 ##  [1]   283 MSNFSMETLKLLRQQTGVGLT...KTSGNSVEVKEFILWKIGA* lcl|CDS.103|gb|CP...
 ##  [2]   435 MNKDLVAIFDYMEKEKGIQRP...EQVSKYGEGKVDEKPKIED* lcl|CDS.950|gb|CP...
 ##  [3]   395 MSKETFQRTKPHINIGTIGHV...AIREGGRTIGAGTISKIIA* lcl|CDS.1406|gb|C...
 ##  [4]   191 MVRVSTSEFRVGLRIEIDGQP...GEVVKVDTRTGSYESRVSK* lcl|CDS.1750|gb|C...
 ##  [5]   718 VDYLEKLQVLIDEEQPSSFFN...LQFQGKKYKISRIQSIWDA* lcl|CDS.9|gb|CP00...
 ##  [6]   695 MSDQEFDLSKIRNIGIMAHID...EPAFFAKVPQKIQEEIVKK* lcl|CDS.408|gb|CP...
 ##  [7]   186 MVLSSQLSVGMFISTKDGLYK...EIGDIIKIDTRTCEYIQRV* lcl|CDS.1126|gb|C...
```




We can use the `ranges()` method to extract `GRanges` objects defined in the
Bioconductor package `GenomicRanges`:



```r
cds_ranges <- biofiles::ranges(cds, include = c("locus_tag", "protein_id", "product"))

 ##  Error: incorrect number of dimensions
r
cds_ranges

 ##  Error: object 'cds_ranges' not found
```




