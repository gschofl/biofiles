#ifndef __BIOFILES__ // make sure the header is included only once
#define __BIOFILES__

#include <Rcpp.h>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>
#include <regex>

struct BaseSpan {
    std::vector<unsigned int> range;
    std::vector<bool> fuzzy;
    int strand;
    std::string accession;
    bool remote;
    std::string type;
};

void parse_simple_span(
    BaseSpan& bs,
    std::string& simple_span, 
    std::string& accession );
    
void parse_gb_location(
    Rcpp::S4 &location,
    std::string &span,
    std::string &accession );

void parse_gb_feature_table(
    Rcpp::S4& gb_feature,
    const std::vector<std::string>& feature_string,
    std::string& accession );

Rcpp::CharacterVector parse_gb_qualifiers( 
    const std::vector<std::string>& qualifiers );
    
static inline unsigned int extractNumber(const std::string& str);

static inline std::string ltrim(std::string str, std::string ch);

static inline std::string rtrim(std::string str, std::string ch);

static inline std::string trim(std::string str, std::string ch);

/*
 * The notation '55^56' describes a site between two adjoining nucleotides,
 * such as endonucleolytic cleavage sites. The permitted formats for 
 * this descriptor are n^n+1, or, for circular molecules, n^1, where "n" 
 * is the full length of the molecule.
 *
 * As of 2006 the use of the notation '12.21' to indicate a single base taken
 * from between the indicated points is deprecated. We are not gonna support
 * it.
 *
 * Sequence spans are indicated by e.g. '34..456'. The '<' and '>' symbols are 
 * used with the starting and ending base numbers to indicate that an end 
 * point is beyond the specified base number.
 *
 * A location in a remote entry is specified by e.g. 'J12345.1:1..15'.
 *
 * Complemented and compound locations:
 *   complement(location)
 *   join(location,location,...)
 *   order(location,location,...)
 *   bond(location,location,...) only in GenPept files
 *   gap(), gap(X), gap(unkX)    in contigs
 */

// complement
static const std::regex COMPL("^complement");

// gap
static const std::regex GAP("^gap");

// gap length
static const std::regex GAPLEN("gap\\((unk)?(\\d+)?\\)");

// compound
static const std::regex CMPND("(join|order|bond)");

// BETWEEN_BASES type
static const std::regex BETWEEN_BASES("\\d+\\^\\d+");

/* 
 split bases
 IMGT/HLA uses fuzzy location strings like "<1..546>" although they should
 read "<1..>546". Aaargh!
 */
static const std::regex BASESPLIT(
    "([<>]?\\d+>?)(\\.\\.|\\^)?(>?\\d+>?)?", std::regex::optimize);

// FUZZY start and make sure we also catch the aberrant usage in IMGT/HLA
static const std::regex FUZZY_START("^<\\d+$");

// FUZZY end
static const std::regex FUZZY_END("^>\\d+$|^\\d+>$");

// split on comma
static const std::regex COMMA("\\s*,\\s*");

// remote accession
static const std::regex RA(
    "[a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?\\:", std::regex::optimize);

// simple location possibly with remote accession
static const std::regex RASL(
    "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"
    "(<?\\d+\\.\\.>?\\d+>?|\\d+\\^\\d+|[<>]?\\d+>?)", std::regex::optimize);

// Possibly complemented simple location
static const std::regex PCSL(
    "^("
        "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"  
        "(<?\\d+\\.\\.>?\\d+>?|\\d+\\^\\d+|[<>]?\\d+>?)"
    "|"
        "(complement|gap)\\("
            "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"
            "(<?\\d+\\.\\.>?\\d+>?|\\d+\\^\\d+|[<>]?\\d+>?|(unk)?(\\d+)?)"
        "\\)"
    ")$", std::regex::optimize);

// Simple locations in compound
static const std::regex SLC(
    "("
        "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"  
        "(<?\\d+\\.\\.>?\\d+>?|\\d+\\^\\d+|[<>]?\\d+>?)"
    "|"
        "(complement|gap)\\("
            "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"
            "(<?\\d+\\.\\.>?\\d+>?|\\d+\\^\\d+|[<>]?\\d+>?|(unk)?(\\d+)?)"
        "\\)"
    ")"
    "(,"
        "("
            "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"  
            "(<?\\d+\\.\\.>?\\d+>?|\\d+\\^\\d+|[<>]?\\d+>?)"
        "|"
            "(complement|gap)\\("
                "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"
                "(<?\\d+\\.\\.>?\\d+>?|\\d+\\^\\d+|[<>]?\\d+>?|(unk)?(\\d+)?)"
            "\\)"
        ")"
    ")*", std::regex::optimize);

// compound location
static const std::regex CL(
    "(join|order|bond)\\("
        "("
            "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"  
            "(<?\\d+\\.\\.>?\\d+>?|\\d+\\^\\d+|[<>]?\\d+>?)"
        "|"
            "(complement|gap)\\("
                "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"
                "(<?\\d+\\.\\.>?\\d+>?|\\d+\\^\\d+|[<>]?\\d+>?|(unk)?(\\d+)?)"
            "\\)"
        ")"
        "(,"
            "("
                "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"  
                "(<?\\d+\\.\\.>?\\d+>?|\\d+\\^\\d+|[<>]?\\d+>?)"
            "|"
                "(complement|gap)\\("
                    "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"
                    "(<?\\d+\\.\\.>?\\d+>?|\\d+\\^\\d+|[<>]?\\d+>?|(unk)?(\\d+)?)"
                "\\)"
            ")"
        ")*"
    "\\)", std::regex::optimize);

// possibly complemented compound location
static const std::regex PCCL(
    "^("
        "(join|order|bond)\\("
            "("
                "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"  
                "(<?\\d+\\.\\.>?\\d+>?|\\d+\\^\\d+|[<>]?\\d+>?)"
            "|"
                "(complement|gap)\\("
                    "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"
                    "(<?\\d+\\.\\.>?\\d+>?|\\d+\\^\\d+|[<>]?\\d+>?|(unk)?(\\d+)?)"
                "\\)"
            ")"
            "(,"
                "("
                      "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"  
                      "(<?\\d+\\.\\.>?\\d+>?|\\d+\\^\\d+|[<>]?\\d+>?)"
                "|"
                    "(complement|gap)\\("
                        "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"
                        "(<?\\d+\\.\\.>?\\d+>?|\\d+\\^\\d+|[<>]?\\d+>?|(unk)?(\\d+)?)"
                    "\\)"
                ")"
            ")*"
        "\\)"
    "|"
        "complement\\("
            "(join|order|bond)\\("
                "("
                    "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"  
                    "(<?\\d+\\.\\.>?\\d+>?|\\d+\\^\\d+|[<>]?\\d+>?)"
                "|"
                    "complement\\("
                        "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"
                        "(<?\\d+\\.\\.>?\\d+>?|\\d+\\^\\d+|[<>]?\\d+>?)"
                    "\\)"
                ")"
                "(,"
                    "("
                        "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"  
                        "(<?\\d+\\.\\.>?\\d+>?|\\d+\\^\\d+|[<>]?\\d+>?)"
                    "|"
                        "complement\\("
                            "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"
                            "(<?\\d+\\.\\.>?\\d+>?|\\d+\\^\\d+|[<>]?\\d+>?)"
                        "\\)"
                    ")"
                ")*"
            "\\)"
        "\\)"
    ")$", std::regex::optimize);

#endif // __BIOFILES__
