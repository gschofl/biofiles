#include <Rcpp.h>
#include <string>
#include <vector>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string_regex.hpp>

using namespace Rcpp;
using namespace std;

void print( std::vector<std::string> &v ) {
    for (size_t n = 0; n < v.size(); n++)
        std::cout << "\"" << v[ n ] << "\"\n";
    std::cout << std::endl;
}

// complement
static const boost::regex COMPL("^complement");

// compound
static const boost::regex CMPND("(join|order)");

// split bases 
static const boost::regex BASESPLIT(
    "([<>]?\\d+)(\\.\\.?|\\^)?([<>]?\\d+)?");

// remote accession
static const boost::regex RA(
    "([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)");
    
static const boost::regex RA2(
    "([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:");
  
// simple location possibly with remote accession
static const boost::regex SL(
    "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"
    "([<>]?\\d+\\.\\.[<>]?\\d+|\\d+\\^\\d+|\\d+\\.\\d+|[<>]?\\d+)");

// Possibly complemented simple location
static const boost::regex PCSL(
    "^("
        "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"  
        "([<>]?\\d+\\.\\.[<>]?\\d+|\\d+\\^\\d+|\\d+\\.\\d+|[<>]?\\d+)"
    "|"
        "complement\\("
            "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"
            "([<>]?\\d+\\.\\.[<>]?\\d+|\\d+\\^\\d+|\\d+\\.\\d+|[<>]?\\d+)"
        "\\)"
    ")$");

// Simple locations in compound
static const boost::regex SLC(
    "("
        "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"  
        "([<>]?\\d+\\.\\.[<>]?\\d+|\\d+\\^\\d+|\\d+\\.\\d+|[<>]?\\d+)"
    "|"
        "complement\\("
            "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"
            "([<>]?\\d+\\.\\.[<>]?\\d+|\\d+\\^\\d+|\\d+\\.\\d+|[<>]?\\d+)"
        "\\)"
    ")"
    "(,"
        "("
            "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"  
            "([<>]?\\d+\\.\\.[<>]?\\d+|\\d+\\^\\d+|\\d+\\.\\d+|[<>]?\\d+)"
        "|"
            "complement\\("
                "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"
                "([<>]?\\d+\\.\\.[<>]?\\d+|\\d+\\^\\d+|\\d+\\.\\d+|[<>]?\\d+)"
            "\\)"
        ")"
    ")*");

// compound location
static const boost::regex CL(
    "(join|order)\\("
        "("
            "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"  
            "([<>]?\\d+\\.\\.[<>]?\\d+|\\d+\\^\\d+|\\d+\\.\\d+|[<>]?\\d+)"
        "|"
            "complement\\("
                "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"
                "([<>]?\\d+\\.\\.[<>]?\\d+|\\d+\\^\\d+|\\d+\\.\\d+|[<>]?\\d+)"
            "\\)"
        ")"
        "(,"
            "("
                "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"  
                "([<>]?\\d+\\.\\.[<>]?\\d+|\\d+\\^\\d+|\\d+\\.\\d+|[<>]?\\d+)"
            "|"
                "complement\\("
                    "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"
                    "([<>]?\\d+\\.\\.[<>]?\\d+|\\d+\\^\\d+|\\d+\\.\\d+|[<>]?\\d+)"
                "\\)"
            ")"
        ")*"
    "\\)");

// possibly complemented compound location
static const boost::regex PCCL(
    "^("
        "(join|order)\\("
            "("
                "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"  
                "([<>]?\\d+\\.\\.[<>]?\\d+|\\d+\\^\\d+|\\d+\\.\\d+|[<>]?\\d+)"
            "|"
                "complement\\("
                    "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"
                    "([<>]?\\d+\\.\\.[<>]?\\d+|\\d+\\^\\d+|\\d+\\.\\d+|[<>]?\\d+)"
                "\\)"
            ")"
            "(,"
                "("
                      "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"  
                      "([<>]?\\d+\\.\\.[<>]?\\d+|\\d+\\^\\d+|\\d+\\.\\d+|[<>]?\\d+)"
                "|"
                    "complement\\("
                        "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"
                        "([<>]?\\d+\\.\\.[<>]?\\d+|\\d+\\^\\d+|\\d+\\.\\d+|[<>]?\\d+)"
                    "\\)"
                ")"
            ")*"
        "\\)"
    "|"
        "complement\\("
            "(join|order)\\("
                "("
                    "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"  
                    "([<>]?\\d+\\.\\.[<>]?\\d+|\\d+\\^\\d+|\\d+\\.\\d+|[<>]?\\d+)"
                "|"
                    "complement\\("
                        "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"
                        "([<>]?\\d+\\.\\.[<>]?\\d+|\\d+\\^\\d+|\\d+\\.\\d+|[<>]?\\d+)"
                    "\\)"
                ")"
                "(,"
                    "("
                        "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"  
                        "([<>]?\\d+\\.\\.[<>]?\\d+|\\d+\\^\\d+|\\d+\\.\\d+|[<>]?\\d+)"
                    "|"
                        "complement\\("
                            "(([a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z0-9]+)?)\\:)?"
                            "([<>]?\\d+\\.\\.[<>]?\\d+|\\d+\\^\\d+|\\d+\\.\\d+|[<>]?\\d+)"
                        "\\)"
                    ")"
                ")*"
            "\\)"
        "\\)"
    ")$");

class BaseSpan {
public:
    int strand;
    std::string accn;
    bool remote;
    std::vector<bool> closed;
    std::vector<bool> partial;
    std::vector<int> position;
};

BaseSpan parse_simple_span(std::string base_span) {

    BaseSpan res;
    std::string::const_iterator start, end;
    start = base_span.begin();
    end = base_span.end();
    boost::smatch m;
    const string empty("");
  
    // test for strand
    int strand(1);
    if ( boost::regex_search(start, end, m, COMPL) ) {
        strand = -1;
    }

    // get span string
    boost::regex_search(start, end, m, SL);
    string span = m[0];
    start = span.begin();
    end = span.end();

    // get remote accession number
    std::string accn("");
    bool remote(false);
    boost::regex_search(start, end, m, RA);
    if (m[0].matched) {
        accn = m[0];
        remote = true;
    }
  
    // strip remote accession number from span
    span = boost::regex_replace(span, RA2, empty);
    start = span.begin();
    end = span.end();
  
    // get closed
    std::vector<bool> closed(2);
    if( not boost::regex_search(start, end, m, boost::regex("\\d+\\.\\d+")) ) {
        closed[0] = true;
        closed[1] = true;
    }
  
    // split span
    string start_pos, end_pos;
    boost::regex_search(start, end, m, BASESPLIT);
    start_pos = m[1];
    end_pos = m[1];
    if ( m[3].matched ) {
        end_pos = m[3];
    }
  
    // get partial
    std::vector<bool> partial(2);
    boost::regex_match(start_pos, m, boost::regex("^(<|>)\\d+") );
    partial[0] = m[0].matched;
    boost::regex_match(end_pos, m, boost::regex("^(<|>)\\d+") );
    partial[1] = m[0].matched;
  
    // get span
    std::vector<int> pos(2);
    static const boost::regex P("<|>");
    pos[0] = atoi( boost::regex_replace(start_pos, P, empty).c_str() );
    pos[1] = atoi( boost::regex_replace(end_pos, P, empty).c_str() );
    
    res.strand = strand;
    res.accn = accn;
    res.remote = remote;
    res.closed = closed;
    res.partial = partial;
    res.position = pos;

    return res;
}

// [[Rcpp::export]]
SEXP parse_gb_location(std::string gb_base_span) {
    
    // clean up anny possible whitespace
    static const string empty("");
    gb_base_span = boost::regex_replace(gb_base_span, boost::regex("\\s+"), empty);
    // std::cout << gb_base_span << std::endl;

    // iterator over gb_base_span
    std::string::const_iterator start, end;
    start = gb_base_span.begin();
    end = gb_base_span.end();
    boost::smatch m;
  
    // use Rcpp::Language to create and assign a 'gbLocation' object
    Rcpp::S4 obj = Rcpp::S4("gbLocation");
  
    // test for a possibly complemented simple location
    if ( boost::regex_match(start, end, m, PCSL) ) {
        BaseSpan s = parse_simple_span(gb_base_span);
        obj.slot("strand") = s.strand;
        obj.slot("accession") = s.accn;
        obj.slot("remote") = s.remote;
        
        LogicalMatrix clo(1,2);
        clo(0,0) = s.closed[0];
        clo(0,1) = s.closed[1];
        obj.slot("closed") = clo;

        LogicalMatrix part(1,2);
        part(0,0) = s.partial[0];
        part(0,1) = s.partial[1];
        obj.slot("partial") = part;

        NumericMatrix pos(1,2);
        pos(0,0) = s.position[0];
        pos(0,1) = s.position[1];
        obj.slot(".Data") = pos;

    } else if ( boost::regex_match( start, end, m, PCCL ) ) {
        // test for complementary strand
        int strand(1);
        if ( boost::regex_search( start, end, m, COMPL ) ) {
            strand = -1;
        }  
        
        // get compound span
        boost::regex_search( start, end, m, CL );
        std::string cmpnd_span( m[0] );
        start = cmpnd_span.begin();
        end = cmpnd_span.end();  
        
        // get compound type
        boost::regex_search(start, end, m, CMPND);
        std::string compound( m[0] );
        
        // get span strings
        boost::regex_search( start, end, m, SLC );
        std::string span_str( m[0] );
        
        vector<string> spans;
        boost::split_regex( spans, span_str, boost::regex( "\\s*,\\s*" ) );
        
        int nrows = spans.size();
        Rcpp::NumericMatrix posMat( nrows, 2 );
        Rcpp::LogicalMatrix partialMat( nrows, 2 );
        Rcpp::LogicalMatrix closedMat( nrows, 2 );
        Rcpp::IntegerVector strandVec(nrows);
        Rcpp::CharacterVector accnVec(nrows);
        Rcpp::LogicalVector remoteVec(nrows);
        
        for (int i = 0; i < nrows; ++i) {
            std::string span =  spans[i];
            // std::cout << span << std::endl;

            BaseSpan s = parse_simple_span(span);
            
            // get positions
            posMat(i, 0) = s.position[0];
            posMat(i, 1) = s.position[1];
            
            // get strand
            strandVec[i] = s.strand;
            
            // get partial
            partialMat(i, 0) = s.partial[0];
            partialMat(i, 1) = s.partial[1];
            
            // get accession and remote
            accnVec[i] = s.accn;
            remoteVec[i] = s.remote;
            
            // get closed
            closedMat(i, 0) = s.closed[0];
            closedMat(i, 1) = s.closed[1];
        }
        
        // if the whole span is complementary reset strandVec
        if (strand == -1) {
            for (int i = 0; i < nrows; ++i)
                strandVec[i] = -1;
        }
        
        obj.slot("compound") = compound;
        obj.slot("strand") = strandVec;
        obj.slot("accession") = accnVec;
        obj.slot("remote") = remoteVec;
        obj.slot("closed") = closedMat;
        obj.slot("partial") = partialMat;
        obj.slot(".Data") = posMat;
    } else {
        std::cout << "Not a valid GenBank base span." << std::endl;
    }
    
    return obj;
}

// [[Rcpp::export]]
Rcpp::CharacterVector get_qual( std::vector<std::string> lines ) {
    Rcpp::CharacterVector quals;
    std::vector<std::string> tag_qual;
    std::string tag("");
    std::string qual("");
    int n = lines.size();
    for(int i = 0; i < n; ++i) {
        std::string line = lines[i];
        boost::split_regex( tag_qual, line, boost::regex( "=" ) );
        tag = tag_qual.at(0).erase(0,1);
        try {
            qual = tag_qual.at(1);
            qual.erase(remove(qual.begin(), qual.end(), '\"'), qual.end());
        }
        catch(...) {
            qual = "TRUE";
        }

        if (tag == "translation") {
            qual.erase(remove(qual.begin(), qual.end(), ' '), qual.end());
        }
    quals[tag] = qual;
    }

    return quals; 
}

// [[Rcpp::export]]
SEXP parse_feature_table(int id = 0,
                         Rcpp::CharacterVector lines = Rcpp::CharacterVector::create(""),
                         Rcpp::S4 seqinfo = Rcpp::S4("gbInfo") ) {
    
    Rcpp::S4 obj = Rcpp::S4("gbFeature");
    obj.slot(".Info") = seqinfo;
    obj.slot(".Id") = id;
    
    std::string::const_iterator s, e;
    boost::smatch m;
    Rcpp::CharacterVector thisline;
    
    // initialize idx
    vector<int> idx;
    idx.push_back(0);
    int end = lines.size();
    for (int i = 0; i < end; ++i) {
        thisline = lines[i];
        std::string line = as<std::string>(thisline);
        s = line.begin();
        e = line.end();
        if ( boost::regex_search(s, e, m, boost::regex("^\\s{21}/")) ) {
            idx.push_back(i);
            // cout << i << endl;
        } 
    }
    idx.push_back(end);
    
    // merge lines
    vector<std::string> merged_lines;
    for (int i = 0; i < idx.size() - 1; ++i) {
        int j(idx[i]);
        string line("");
        while (j < idx[i + 1]) {
            thisline = lines[j];
            string newline = as<string>(thisline);
            boost::trim( newline );
            line += newline + " ";
            j++;
        }
        boost::trim( line );
        merged_lines.push_back(line);
    }
    
    // split key and location
    std::vector<std::string> key_loc;
    boost::split_regex( key_loc, merged_lines[0], boost::regex( "\\s+" ) );
    merged_lines.erase( merged_lines.begin() );

    std::string key = key_loc[0];
    std::string loc("");
    for (int i = 1; i < key_loc.size(); ++i) {
        loc += key_loc[i];
    }
    // cout << key << "\n" << loc << endl;
    obj.slot("key") = key;
    obj.slot("location") = parse_gb_location( loc );
    obj.slot("qualifiers") = get_qual( merged_lines );
    
    return obj;
}

