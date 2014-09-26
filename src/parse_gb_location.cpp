#include "biofiles.h"
#include <boost/algorithm/string_regex.hpp>

using namespace Rcpp;

// [[Rcpp::export]]
SEXP gbLocation(
    std::string gb_base_span,
    std::string accession = "" ) 
{
    // create and assign a 'gbLocation' object
    Rcpp::S4 gb_location = Rcpp::S4("gbLocation");
    parse_gb_location( gb_location, gb_base_span, accession );
    return gb_location; 
}

void parse_gb_location(
    Rcpp::S4 &location,
    std::string &span,
    std::string &accession )
{
    // clean up possible whitespace
    span.erase( remove_if(span.begin(), span.end(), ::isspace), span.end());

    // iterator over the complete gb_base_span
    std::string::const_iterator b_it, e_it; 
    b_it = span.begin();
    e_it = span.end();
    boost::smatch m;
  
    // initialise a BaseSpan
    BaseSpan bs;
  
    // test for a possibly complemented simple location
    try {
      if ( boost::regex_match(b_it, e_it, m, PCSL) ) {
        
          parse_simple_span(bs, span, accession);
          
          // fill gbLocation@range
          // guarantied integer Matrix/two columns 
          IntegerMatrix range(1,2);
          range(0,0) = bs.range[0];
          range(0,1) = bs.range[1];
          location.slot("range") = range;
          
          // fill gbLocation@fuzzy
          // guarantied logical Matrix/two columns 
          LogicalMatrix fuzzy(1,2);
          fuzzy(0,0) = bs.fuzzy[0];
          fuzzy(0,1) = bs.fuzzy[1];
          location.slot("fuzzy") = fuzzy;
          
          location.slot("strand") = bs.strand;
          location.slot("accession") = bs.accession;
          location.slot("remote") = bs.remote;
          location.slot("type") = bs.type;
          
      // test for a possibly complemented compound location
      } else if ( boost::regex_match( b_it, e_it, m, PCCL ) ) {
          // test for complementary strand
          int strand(1);
          if ( boost::regex_search( b_it, e_it, m, COMPL ) ) {
              strand = -1;
          }  
          
          // get compound span
          boost::regex_search( b_it, e_it, m, CL );
          b_it = m[0].first;
          e_it = m[0].second;  
          
          // get compound type
          boost::regex_search( b_it, e_it, m, CMPND );
          std::string compound( m[0] );
          
          // get span strings
          boost::regex_search( b_it, e_it, m, SLC );
          std::string span_str( m[0] );
          
          std::vector<std::string > spans;
          boost::split_regex( spans, span_str, boost::regex("\\s*,\\s*") );
          
          int nrows = spans.size();
          Rcpp::IntegerMatrix   rangeMat( nrows, 2 );
          Rcpp::LogicalMatrix   fuzzyMat( nrows, 2 );
          Rcpp::IntegerVector   strandVec(nrows);
          Rcpp::CharacterVector accessionVec(nrows);
          Rcpp::LogicalVector   remoteVec(nrows);
          Rcpp::CharacterVector typeVec(nrows);
            
          for (int i = 0; i < nrows; ++i) {
              std::string span = spans[i];
              parse_simple_span( bs, span, accession );
              
              // get positions
              rangeMat(i, 0) = bs.range[0];
              rangeMat(i, 1) = bs.range[1];
              // get fuzzy
              fuzzyMat(i, 0) = bs.fuzzy[0];
              fuzzyMat(i, 1) = bs.fuzzy[1];
              // get strand, accession, remote and type
              strandVec[i] = bs.strand;
              accessionVec[i] = bs.accession;
              remoteVec[i] = bs.remote;
              typeVec[i] = bs.type;
          }
          
          // if the whole span is complementary reset strandVec
          if (strand == -1) {
              for (int i = 0; i < nrows; ++i)
                  strandVec[i] = -1;
          }
          
          location.slot("range") = rangeMat;
          location.slot("fuzzy") = fuzzyMat;
          location.slot("strand") = strandVec;
          location.slot("compound") = compound;
          location.slot("accession") = accessionVec;
          location.slot("remote") = remoteVec;
          location.slot("type") = typeVec;
      } else {
        throw std::range_error("Cannot parse location descriptor.");
      }
    } catch (std::exception &ex) {
      forward_exception_to_r(ex);
    }
}

void parse_simple_span(
    BaseSpan &bs,
    std::string &simple_span,
    std::string &accession )
{
    std::string::const_iterator b_it, e_it;
    b_it = simple_span.begin();
    e_it = simple_span.end();
    boost::smatch m;
    //Rcpp::Rcout << "'" << simple_span << "' => ";
    //Rcpp::Rcout << *b_it << " ... " << *e_it << std::endl;
    
    // initialize and set defaults
    std::string start, end;
    std::vector<unsigned int> range(2);
    std::vector<bool> fuzzy(2, false);
    int strand(1);
    bool remote(false);
    std::string type("R");
  
    // test for complement
    if ( boost::regex_search(b_it, e_it, m, COMPL) ) {
        strand = -1;
    }
    
    // test for gap
    if ( boost::regex_search(b_it, e_it, m, GAP) ) {
        // we have matched a gap -> gap(), gap(X), gap(unkX)
        type = "G";
        accession = "";
        start = "1";
        // gap length
        boost::regex_search(b_it, e_it, m, GAPLEN);
        if (m[1].matched /* gap(unkX) */ | !m[2].matched /* gap() */) {
            fuzzy[0] = true;
        }
        if (m[2].matched /* gap(34) | gap(unk34) */) {
            end = m[2];
        } else { /* gap() */
            end = "1";
        }
    } else {
        // match remote accession (capture group 1) and
        // genomic span (capture group 4)
        boost::regex_search(b_it, e_it, m, RASL);
        if (m[1].matched) {
            remote = true;
            accession = m[2];  // remote accession without colon
            b_it = m[4].first;
            e_it = m[4].second;
        } else {
            b_it = m[0].first;
            e_it = m[0].second;
        }
        //Rcpp::Rcout << accn << std::endl;
        //Rcpp::Rcout << "'" << simple_span << "' => ";
        //Rcpp::Rcout << *b_it << " ... " << *e_it << std::endl;

        // get type
        if ( boost::regex_search(b_it, e_it, m, BETWEEN_BASES) ) {
            type = "B";
        }
  
        // split span
        boost::regex_search(b_it, e_it, m, BASESPLIT);
        start = m[1];
        if ( m[3].matched ) {
            end = m[3];
        } else {
            end = start;
        }
        
        //Rcpp::Rcout << "'" << simple_span << "' => ";
        //Rcpp::Rcout << start << " ... " << end << std::endl;
        // get fuzzy and make sure we also catch the usage in IMGT/HLA
        static const boost::regex FUZZY_START("^<\\d+$");
        static const boost::regex FUZZY_END("^>\\d+$|^\\d+>$");
        fuzzy[0] = boost::regex_match(start, m, FUZZY_START );
        fuzzy[1] = boost::regex_match(end, m, FUZZY_END );
    }

    // get range
    range[0] = extractNumber(start);
    range[1] = extractNumber(end);

    try {
      // throw error if the end point comes before the start point.
      if (range[0] > range[1])
        throw std::range_error("Inadmissible range: start point is larger than end point.");
      // throw error if type = 'B' and end point is not adjacent to start point.
      if ((type == "B") & (range[1] - range[0] != 1))
        throw std::range_error("Inadmissible range: For span of type '36^37', start and end positions must be adjacent");
      
      bs.range = range;
      bs.fuzzy = fuzzy;
      bs.strand = strand;
      bs.accession = accession;
      bs.remote = remote;
      bs.type = type;
      
    } catch (std::exception &ex) {
      forward_exception_to_r(ex);
    }
}


unsigned int extractNumber(const std::string& str) {
    std::string temp;
    std::string::const_iterator it;
    for (it = str.begin(); it != str.end(); it++)
    {
        if (isdigit(*it)) {
            temp += *it;
        }
    }
    return( atoi(temp.c_str()) );
}

