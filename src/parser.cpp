#include <Rcpp.h>
#include <string>
#include <vector>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string_regex.hpp>

using namespace Rcpp;
using namespace std;

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

Rcpp::List parse_simple_span(std::string base_span) {

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
    CharacterVector accn = CharacterVector::create(NA_STRING);
    bool remote(false);
    boost::regex_search(start, end, m, RA);
    if (m[0].matched) {
        accn = CharacterVector::create(m[0]);
        remote = true;
    }
  
    // strip remote accession number from span
    span = boost::regex_replace(span, RA2, empty);
    start = span.begin();
    end = span.end();
  
    // get closed
    LogicalMatrix closed(1, 2);
    if( not boost::regex_search(start, end, m, boost::regex("\\d+\\.\\d+")) ) {
        closed(0, 0) = true;
        closed(0, 1) = true;
    } else {
        closed(0, 0) = false;
        closed(0, 1) = false;
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
    LogicalMatrix partial(1, 2);
    boost::regex_match(start_pos, m, boost::regex("^(<|>)\\d+") );
    partial(0, 0) = m[0].matched;
    boost::regex_match(end_pos, m, boost::regex("^(<|>)\\d+") );
    partial(0, 1) = m[0].matched;
  
    // get span
    IntegerMatrix pos(1, 2);
    static const boost::regex P("<|>");
    pos(0, 0) = atoi( boost::regex_replace(start_pos, P, empty).c_str() );
    pos(0, 1) = atoi( boost::regex_replace(end_pos, P, empty).c_str() );
    
    return List::create(
        _["position"] = pos, _["strand"] = strand, _["partial"] = partial,
        _["accn"] = accn, _["remote"] = remote, _["closed"] = closed
    );
}

// [[Rcpp::export]]
SEXP parse_gb_location(std::string gb_base_span) {
    
    // clean up anny possible whitespace
    static const string empty("");
    gb_base_span = boost::regex_replace(gb_base_span, boost::regex("\\s+"), empty);
    
    // iterator over gb_base_span
    std::string::const_iterator start, end;
    start = gb_base_span.begin();
    end = gb_base_span.end();
    boost::smatch m;
  
    // use Rcpp::Language to create and assign a 'gbLocation' object
    Rcpp::Language loc("new", "gbLocation");
    Rcpp::S4 obj( loc.eval() );
  
    // test for a possibly complemented simple location
    if ( boost::regex_match(start, end, m, PCSL) ) {
        Rcpp::List locList = parse_simple_span(gb_base_span);
        obj.slot("strand") = locList["strand"];
        obj.slot("accession") = locList["accn"];
        obj.slot("remote") = locList["remote"];
        obj.slot("closed") = locList["closed"];
        obj.slot("partial") = locList["partial"];
        obj.slot(".Data") = locList["position"];
        
    } else if ( boost::regex_match( start, end, m, PCCL ) ) {
        // test for complementary strand
        int strand(1);
        if ( boost::regex_search( start, end, m, COMPL ) ) {
            strand = -1;
        }  
        
        // get compound span
        boost::regex_search( start, end, m, CL );
        string cmpnd_span( m[0] );
        start = cmpnd_span.begin();
        end = cmpnd_span.end();  
        
        // get compound type
        boost::regex_search(start, end, m, CMPND);
        string compound( m[0] );
        
        // get span strings
        boost::regex_search( start, end, m, SLC );
        string span_str( m[0] );
        
        vector<string> spans;
        boost::split_regex( spans, span_str, boost::regex( "\\s*,\\s*" ) );
        
        int nrows = spans.size();
        Rcpp::IntegerMatrix posMat( nrows, 2 );
        Rcpp::LogicalMatrix partialMat( nrows, 2 );
        Rcpp::LogicalMatrix closedMat( nrows, 2 );
        Rcpp::IntegerVector strandVec(nrows);
        Rcpp::CharacterVector accnVec(nrows);
        Rcpp::LogicalVector remoteVec(nrows);
        
        for (int i = 0; i < nrows; ++i) {
            string span =  spans[i];
            Rcpp::List locList = parse_simple_span(span);
            
            // get positions
            Rcpp::IntegerMatrix pos = locList["position"];
            posMat(i, 0) = pos(0, 0);
            posMat(i, 1) = pos(0, 1);
            
            // get strand
            strandVec[i] = locList["strand"];
            
            // get partial
            Rcpp::LogicalMatrix part = locList["partial"];
            partialMat(i, 0) = part(0, 0);
            partialMat(i, 1) = part(0, 1);
            
            // get accession and remote
            CharacterVector accn = locList["accn"];
            accnVec[i] = accn[0];
            remoteVec[i] = locList["remote"];
            
            // get closed
            Rcpp::LogicalMatrix clo = locList["closed"];
            closedMat(i, 0) = clo(0, 0);
            closedMat(i, 1) = clo(0, 1);
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
    }    
    return obj;
}

// [[Rcpp::export]]
CharacterVector get_qual( std::vector<std::string> lines ) {
  int n = lines.size();
  Rcpp::CharacterVector quals = CharacterVector::create("");
  std::string qual;

  for(int i = 0; i < n; ++i) {
    std::string line = lines[i];
    int equalLoc = line.find('=');
    std::string tag = line.substr(1, equalLoc - 1);
    if (equalLoc < 0) {
      qual = "TRUE";
    } else {
      qual = line.substr(equalLoc + 1, line.length() - equalLoc - 1);
      qual.erase(remove(qual.begin(), qual.end(),'\"'), qual.end());
    }
    if (tag == "translation") {
      qual.erase(remove(qual.begin(), qual.end(),' '), qual.end());
    }
    quals[tag] = qual;
  }
  
  return quals; 
}

// [[Rcpp::export]]
SEXP parse_feature_table(int id = 0,
                         CharacterVector lines = CharacterVector::create(""),
                         std::string db_dir = "", std::string accession = "",
                         std::string definition = "") {
                           
    // use Rcpp::Language to create and assign a 'gbFeature' object
    Rcpp::Language ftr("new", "gbFeature");
    Rcpp::S4 obj( ftr.eval() );
    obj.slot(".Dir") = db_dir;
    obj.slot(".ACCN") = accession;
    obj.slot(".DEF") = definition;
    obj.slot(".ID") = id;
    
    string::const_iterator s, e;
    boost::smatch m;
    CharacterVector thisline;
    
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
        // cout << line << endl;
        merged_lines.push_back(line);
    }
    
    // split key and location
    std::vector<std::string> key_loc;
    boost::split_regex( key_loc, merged_lines[0], boost::regex( "\\s+" ) );
    merged_lines.erase( merged_lines.begin() );
    // cout << key_loc[0] << "\n" << key_loc[1] << endl;
    
    obj.slot("key") = key_loc[0];
    obj.slot("location") = parse_gb_location( key_loc[1] );
    obj.slot("qualifiers") = get_qual(  merged_lines );
    
    return obj;
}
