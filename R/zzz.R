###
GBFIELDS <- paste0("@G@I|accession|comment|date|dblink|dbsource|",
                   "definition|division|features|keywords|length|",
                   "lineage|locus|organism|references|sequence|",
                   "source|topology|type|version")

###
# Declare the regex patterns for use in .getLocation

# single location possibly fuzzy
sil <- "[<>]?\\d+"
# within location
wl <- "\\d+\\.\\d+"
# between location
bl <- "\\d+\\^\\d+"
# paired location possibly fuzzy
pl <- sprintf("%s\\.\\.%s", sil, sil)
# remote accession
ra <- "([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)"

# simple location possibly with remote accession
sl <- sprintf("(%s\\:)?(%s|%s|%s|%s)", ra, sil, bl, wl, pl)
# complemented simple location
csl <- sprintf("complement\\(%s\\)", sl)
# possibly complemented simple location
pcsl <- sprintf("(%s|%s)", sl, csl)

# compound location
cl <- sprintf("(join|order)\\(%s(,%s)*\\)", pcsl, pcsl)
# complemented compound location
ccl <- sprintf("complement\\(%s\\)", cl)
