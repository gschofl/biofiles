# biofiles - an interface to GenBank/GenPept files in R

[![Travis-CI Build Status](https://travis-ci.org/gschofl/biofiles.svg?branch=master)](https://travis-ci.org/gschofl/biofiles)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/gschofl/biofiles?branch=master&svg=true)](https://ci.appveyor.com/project/gschofl/biofiles)

This is an R package for interfacing with GenBank and GenPept flat
file records. It includes utilities for reading and writing GenBank
files, and methods for interacting with annotation and sequence data.

### Installation

This package is currently ony available via github. 
It depends on the Boost regex library and I have so far only tested it on Ubuntu.

On Ubuntu run
```sh
sudo apt install libboost-regex-dev
```
before attempting to install `biofiles`.


```r
install.packages("devtools")  # if not already installed
devtools::install_github("gschofl/biofiles")
```





### Basic functionality

Let's download a small bacterial genome, _Chlamydophila psittaci_ DC15
(GenBank acc. CP002806; GI 334693792) from NCBI.



```r
library(reutils)
gb_file <- reutils::efetch("CP002806", "nuccore", rettype = "gbwithparts", retmode = "text")
```




Next, we parse the `efetch` object into a `gbRecord` instance.



```r
rec <- gbRecord(gb_file)
rec

 ##  An object of class 'gbRecord', with 2040 features
 ##  LOCUS       CP002806             1172182 bp    DNA     circular BCT 31 ...
 ##  DEFINITION  Chlamydophila psittaci 02DC15, complete genome.
 ##  ACCESSION   CP002806
 ##  VERSION     CP002806.1  GI:334693792
 ##  DBLINK      Project: PRJNA66481
 ##  KEYWORDS    .
 ##  SOURCE      Chlamydia psittaci 02DC15
 ##    ORGANISM  Chlamydia psittaci 02DC15
 ##              Bacteria; Chlamydiae; Chlamydiales; Chlamydiaceae;
 ##              Chlamydia/Chlamydophila group; Chlamydia.
 ##  REFERENCE   1  (bases 1 to 1172182)
 ##    AUTHORS   Schofl,G., Voigt,A., Litsche,K., Sachse,K. and Saluz,H.P.
 ##    TITLE     Complete Genome Sequences of Four Mammalian Isolates of
 ##              Chlamydophila psittaci
 ##    JOURNAL   J. Bacteriol. 193 (16), 4258 (2011)
 ##    PUBMED    21705611
 ##  REFERENCE   2  (bases 1 to 1172182)
 ##    AUTHORS   Schofl,G.
 ##    TITLE     Direct Submission
 ##    JOURNAL   Submitted (20-MAY-2011) Cell and Molecular Biology,
 ##              Leibniz Institute for Natural Product Research and
 ##              Infection Biology, Beutenbergstrasse 11a, Jena D-07745,
 ##              Germany
 ##  COMMENT     Source DNA available from Dr. Konrad Sachse
 ##              (Konrad.Sachse@fli.bund.de) at the Institute for Molecular
 ##              Pathogenesis, Jena, Germany.
 ##  ORIGIN      CAAAGTTTTAAACATGTTAACCACGTTGTTTTCCCTCTATAAACCGAGTTTTAAACAGT
 ##              ...
 ##              ATCTTGTGTTCACAAACATGCAAAATCTCGATTTATTTACAAACTGAACAATGTTTTAA
 ##  CONTIG
```




The `summary` function provides an overview over the object:



```r
summary(rec)

 ##  [[CP002806]]
 ##    1172182 bp: Chlamydophila psittaci 02DC15, complete genome.
 ##     Id Feature Location                     GeneId Product           ...
 ##      1 source  1..1172182                   NA     NA                ...
 ##      2 gene    149..1144                    hemB   NA                ...
 ##      3 CDS     149..1144                    hemB   delta-aminolevuli ...
 ##      4 gene    complement(1160..2578)       NA     NA                ...
 ##      5 CDS     complement(1160..2578)       NA     Na(+)-translocati ...
 ##      6 gene    complement(2601..3035)       NA     NA                ...
 ##      7 CDS     complement(2601..3035)       NA     conserved hypothe ...
 ##    ... ...     ...                          ...    ...               ...
 ##   2034 CDS     1169125..1170006             NA     geranylgeranyl py ...
 ##   2035 gene    complement(1170021..1170094) NA     NA                ...
 ##   2036 tRNA    complement(1170021..1170094) NA     tRNA-Pro          ...
 ##   2037 gene    complement(1170106..1170360) NA     NA                ...
 ##   2038 CDS     complement(1170106..1170360) NA     conserved hypothe ...
 ##   2039 gene    complement(1170443..1172035) NA     NA                ...
 ##   2040 CDS     complement(1170443..1172035) NA     conserved hypothe ...
```




Various getter methods provide access to the data contained in a GenBank record;
for instance:



```r
getAccession(rec)

 ##  [1] "CP002806"
```





```r
getGeneID(rec)

 ##  [1] "334693792"
```





```r
getDefinition(rec)

 ##  [1] "Chlamydophila psittaci 02DC15, complete genome."
```





```r
getOrganism(rec)

 ##  [1] "Chlamydia psittaci 02DC15"
```





```r
getSequence(rec)

 ##    A DNAStringSet instance of length 1
 ##        width seq                                        names               
 ##  [1] 1172182 CAAAGTTTTAAACATGTTAA...AAACTGAACAATGTTTTAA CP002806
```





```r
getReference(rec)

 ##  A 'gbReferenceList' instance:
 ##  REFERENCE   1  (bases 1 to 1172182)
 ##    AUTHORS   Schofl,G., Voigt,A., Litsche,K., Sachse,K. and Saluz,H.P.
 ##    TITLE     Complete Genome Sequences of Four Mammalian Isolates of
 ##              Chlamydophila psittaci
 ##    JOURNAL   J. Bacteriol. 193 (16), 4258 (2011)
 ##    PUBMED    21705611
 ##  REFERENCE   2  (bases 1 to 1172182)
 ##    AUTHORS   Schofl,G.
 ##    TITLE     Direct Submission
 ##    JOURNAL   Submitted (20-MAY-2011) Cell and Molecular Biology,
 ##              Leibniz Institute for Natural Product Research and
 ##              Infection Biology, Beutenbergstrasse 11a, Jena D-07745,
 ##              Germany
 ##  
```




`uniqueQualifs()` provides an overview over the feature qualifiers used
in a record:



```r
uniqueQualifs(rec)

 ##   [1] "organism"            "mol_type"            "strain"             
 ##   [4] "db_xref"             "gene"                "locus_tag"          
 ##   [7] "note"                "codon_start"         "transl_table"       
 ##  [10] "product"             "protein_id"          "translation"        
 ##  [13] "EC_number"           "artificial_location"
```




### Genbank Features

The important part of GenBank record will generally be the list of annotions
or features.

We can access the `gbFeatureList` of a `gbRecord` using `getFeatures()` or `ft()`:



```r
f <- ft(rec)
f

 ##  'gbFeatureTable' with 2040 features:
 ##  Feature:        Location/Qualifiers:
 ##  source          1..1172182
 ##                  /organism = "Chlamydia psittaci 02DC15"
 ##                  /mol_type = "genomic DNA"
 ##                  /strain = "02DC15"
 ##                  /db_xref = "taxon:1112254"
 ##  ...
 ##  Feature:        Location/Qualifiers:
 ##  CDS             complement(1170443..1172035)
 ##                  /locus_tag = "CPS0B_1080"
 ##                  /note = "[F] COG1351 Predicted alternative thymidylate
 ##                  synthase"
 ##                  /codon_start = "1"
 ##                  /transl_table = "11"
 ##                  /product = "conserved hypothetical protein"
 ##                  /protein_id = "AEG87987.1"
 ##                  /db_xref = "GI:334694770"
 ##                  /translation = "MLSRDDEFSSEQRKSLSHFVTNLETNIFALKNLPEVVKGA
 ##                  LFSKYSRSTLGLRSLLLKEFLEGEGGDFLDSSGVDFEVGIHKAADFYRRVLDGFG
 ##                  DDSIGELGGAHLAIESVSMLAAKILEDARIGGSPLEKSSRYVYFDQKVKGEYLYY
 ##                  RDPILMTSAFKDVFLGTCDFLFDTYADLIPKVRTYFEKIYPKESEVSQSAYTISL
 ##                  RAKVLDCLRGLLPAATLTNLGFFGNGRFWQTLLHKIQGHNLTEIRQIGESSLTEL
 ##                  MKIIPSFVSRAESHHHHHQAMLSYRQTLREQLTSLAEKFSGGSHPSKQTGVRLVY
 ##                  GDPEGIYKVAAGFLFPYSEHTYEELINICKSMPREDLIRVLEAGSSSRENRRHKS
 ##                  PRGLECLEFGFDITADFGAYRDLQRHRILTQERQLLTTNLGYHIPEQLLDTPMEK
 ##                  DFREAMEKAEEAYNQISLEFPEEAQYVVPLAYNIRWFFHINGRALQWLCELRSQV
 ##                  QGHENYRKIAIDMVKEVIRFDPVYESFFKFVDYSECDLGRIKQESRKRSF"
 ##  Seqinfo:
 ##  CP002806  1172182 DNA  Chlamydophila psittaci 02DC15, complete genome.
```




We can extract features either by numeric subsetting:


```r
f[[1]]

 ##  Feature:        Location/Qualifiers:
 ##  source          1..1172182
 ##                  /organism = "Chlamydia psittaci 02DC15"
 ##                  /mol_type = "genomic DNA"
 ##                  /strain = "02DC15"
 ##                  /db_xref = "taxon:1112254"
 ##  Seqinfo:
 ##  CP002806  1172182 DNA  Chlamydophila psittaci 02DC15, complete genome.
```




or we can subset by feature key:


```r
f["CDS"][[2]]

 ##  Feature:        Location/Qualifiers:
 ##  CDS             complement(1160..2578)
 ##                  /locus_tag = "CPS0B_0002"
 ##                  /EC_number = "1.6.5.-"
 ##                  /note = "[C] COG1726 Na+-transporting NADH:ubiquinone
 ##                  oxidoreductase, subunit NqrA"
 ##                  /codon_start = "1"
 ##                  /transl_table = "11"
 ##                  /product = "Na(+)-translocating NADH:ubiquinone
 ##                  oxidoreductase subunit A"
 ##                  /protein_id = "AEG87011.1"
 ##                  /db_xref = "GI:334693794"
 ##                  /translation = "MKIAITRGLDLSLQGSPKESGFLKRIDPALVSVDLRPYSA
 ##                  LTLKLKVEQGQAISSGSPVAEYKNFPGVFITSPVSGTVQEIRRGDKRSLLDVVIK
 ##                  KNPGQNLTEYSYDLSKLSRQELLEIFKKEGLFALFKQRPFDIPALPTHHPRDVFI
 ##                  NLADNRPFTPSTEKHLSVFSSREEGFYVFNVGVRAIAKLFGLCPHIISTDRLAIP
 ##                  EKDLKSIAHLHKITGPYPSGSPSTHIHYIAPITSEKDVVFTISFQEVLAIGHLFL
 ##                  KGRILNEQVVALAGSGLKPSLRRYVITTRGADFQSLLPLDEIASDQVSLISGDPL
 ##                  TGRLCDKEHLPCLGMRDATISVIPNPQKRQAFNFLRLGINKPTLTRTYLSGFLKR
 ##                  KHTYMDPDTNLHGETRPIIDTEIYDKVMAMKIPVVPLIKSVITKNFELACMLGFL
 ##                  EVCPEDFALPTFIDPSKTEMLTIIKEALKHYAKETGILNPENTADTE"
 ##  Seqinfo:
 ##  CP002806  1172182 DNA  Chlamydophila psittaci 02DC15, complete genome.
```




A more versatile method to narrow down the list of features of interest is the 
function `filter()`.
For instance, we can filter for all coding sequences (CDS) with the
annotation "hypothetical" in the product qualifiers:



```r
hypo <- filter(rec, key = "CDS", product = "hypothetical")
summary(hypo)

 ##  [[CP002806]]
 ##    1172182 bp: Chlamydophila psittaci 02DC15, complete genome.
 ##     Id Feature Location                     GeneId Product           ...
 ##      7 CDS     complement(2601..3035)       NA     conserved hypothe ...
 ##     31 CDS     20427..21713                 NA     conserved hypothe ...
 ##     33 CDS     21822..23741                 NA     conserved hypothe ...
 ##     35 CDS     complement(23931..26531)     NA     conserved hypothe ...
 ##     37 CDS     26689..29340                 NA     conserved hypothe ...
 ##     39 CDS     29390..29620                 NA     conserved hypothe ...
 ##     43 CDS     complement(30347..30832)     NA     conserved hypothe ...
 ##    ... ...     ...                          ...    ...               ...
 ##   2002 CDS     1152648..1153925             NA     conserved hypothe ...
 ##   2008 CDS     complement(1155986..1156711) NA     conserved hypothe ...
 ##   2010 CDS     complement(1156746..1157474) NA     conserved hypothe ...
 ##   2024 CDS     1163308..1164144             NA     conserved hypothe ...
 ##   2028 CDS     complement(1165683..1167527) NA     conserved hypothe ...
 ##   2038 CDS     complement(1170106..1170360) NA     conserved hypothe ...
 ##   2040 CDS     complement(1170443..1172035) NA     conserved hypothe ...
r
length(hypo)

 ##  [1] 249
```




or we can filter for all elongation factors,



```r
elong <- filter(rec, key = "CDS", product = "elongation factor")
summary(elong)

 ##  [[CP002806]]
 ##    1172182 bp: Chlamydophila psittaci 02DC15, complete genome.
 ##     Id Feature Location                   GeneId Product             ...
 ##      9 CDS     3148..5301                 greA   transcription elong ...
 ##    103 CDS     complement(56119..56967)   tsf    translation elongat ...
 ##    408 CDS     206213..208297             fusA   translation elongat ...
 ##    950 CDS     complement(526454..527758) nusA   transcription elong ...
 ##   1126 CDS     623470..624027             NA     elongation factor P ...
 ##   1406 CDS     complement(805145..806329) tuf    translation elongat ...
 ##   1750 CDS     complement(997270..997842) efp    translation elongat ...
```




now let's extract the sequence for all elongation factors, and using the tools
from the `Biostrings` packages, translate them into protein sequences. Note, that
in order to do so, we first get the `gbFeatureTable` from the `gbRecord`, as otherwise
we'd just extract the complete sequence associated with the GenBank record.



```r
require(Biostrings)
dna <- getSequence(ft(elong))
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
elong_ranges <- ranges(elong, include = c("locus_tag", "protein_id", "product"))
elong_ranges

 ##  GRanges with 7 ranges and 4 metadata columns:
 ##               seqnames           ranges strand |         key   locus_tag
 ##                  <Rle>        <IRanges>  <Rle> | <character> <character>
 ##    CPS0B_0004 CP002806 [  3148,   5301]      + |         CDS  CPS0B_0004
 ##    CPS0B_0056 CP002806 [ 56119,  56967]      - |         CDS  CPS0B_0056
 ##    CPS0B_0219 CP002806 [206213, 208297]      + |         CDS  CPS0B_0219
 ##    CPS0B_0509 CP002806 [526454, 527758]      - |         CDS  CPS0B_0509
 ##    CPS0B_0600 CP002806 [623470, 624027]      + |         CDS  CPS0B_0600
 ##    CPS0B_0751 CP002806 [805145, 806329]      - |         CDS  CPS0B_0751
 ##    CPS0B_0928 CP002806 [997270, 997842]      - |         CDS  CPS0B_0928
 ##                protein_id
 ##               <character>
 ##    CPS0B_0004  AEG87013.1
 ##    CPS0B_0056  AEG87055.1
 ##    CPS0B_0219  AEG87200.1
 ##    CPS0B_0509  AEG87466.1
 ##    CPS0B_0600  AEG87550.1
 ##    CPS0B_0751  AEG87684.1
 ##    CPS0B_0928  AEG87849.1
 ##                                                           product
 ##                                                       <character>
 ##    CPS0B_0004 transcription elongation factor GreA domain protein
 ##    CPS0B_0056                    translation elongation factor Ts
 ##    CPS0B_0219                     translation elongation factor G
 ##    CPS0B_0509                     transcription elongation factor
 ##    CPS0B_0600                                 elongation factor P
 ##    CPS0B_0751                    translation elongation factor Tu
 ##    CPS0B_0928                     translation elongation factor P
 ##    ---
 ##    seqlengths:
 ##     CP002806
 ##      1172182
```




