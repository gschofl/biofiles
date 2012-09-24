###
GBFIELDS <- paste0("@G@I|accession|comment|date|dblink|dbsource|",
                   "definition|division|features|keywords|length|",
                   "lineage|locus|organism|references|sequence|",
                   "source|topology|type|version")

.GBFIELDS <- c("@G@I","accession","comment","date","dblink",
               "dbsource","definition","division","features",
               "keywords","length","lineage","locus","organism",
               "references","sequence","source","topology",
               "type","version")

###
# Declare the regex patterns for use in .getLocation
# single location possibly fuzzy
.SIL <- "[<>]?\\d+"
# within location
.WL <- "\\d+\\.\\d+"
# between location
.BL <- "\\d+\\^\\d+"
# paired location possibly fuzzy
.PL <- sprintf("%s\\.\\.%s", .SIL, .SIL)
# remote accession
.RA <- "([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)"
# simple location possibly with remote accession
.SL <- sprintf("(%s\\:)?(%s|%s|%s|%s)", .RA, .SIL, .BL, .WL, .PL)
# complemented simple location
.CSL <- sprintf("complement\\(%s\\)", .SL)
# possibly complemented simple location
.PCSL <- sprintf("(%s|%s)", .SL, .CSL)
# compound location
.CL <- sprintf("(join|order)\\(%s(,%s)*\\)", .PCSL, .PCSL)
# complemented compound location
.CCL <- sprintf("complement\\(%s\\)", .CL)
