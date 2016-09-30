#include "biofiles.h"

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
    span.erase(remove_if(begin(span), end(span), ::isspace), end(span));

    // initialise a BaseSpan
    BaseSpan bs;
    
    // initialise match results
    std::smatch m;
    
    // iterator over span
    std::string::const_iterator first = span.begin();
    std::string::const_iterator last = span.end();
  
    // test for a possibly complemented simple location
    try {
      if (std::regex_match(first, last, m, PCSL)) {
        
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
      } else if (std::regex_match(first, last, m, PCCL)) {
          // test for complementary strand
          int strand(1);
          if (std::regex_search(first, last, m, COMPL)) {
              strand = -1;
          }  
          
          // get compound span
          std::regex_search(first, last, m, CL);
          first = m[0].first;
          last = m[0].second;  
          
          // get compound type
          std::regex_search(first, last, m, CMPND);
          std::string compound( m[0] );
          
          // get span strings
          std::regex_search(first, last, m, SLC);
          std::string span_str( m[0] );
          
          std::vector<std::string > spans;
          std::sregex_token_iterator
            first_it{begin(span_str), end(span_str), COMMA, -1},
            last_it;
          std::copy(first_it, last_it, std::back_inserter(spans));
          
          int nrows = spans.size();
          Rcpp::IntegerMatrix   rangeMat( nrows, 2 );
          Rcpp::LogicalMatrix   fuzzyMat( nrows, 2 );
          Rcpp::IntegerVector   strandVec(nrows);
          Rcpp::CharacterVector accessionVec(nrows);
          Rcpp::LogicalVector   remoteVec(nrows);
          Rcpp::CharacterVector typeVec(nrows);
            
          for (int i = 0; i < nrows; ++i) {
              std::string span = spans[i];
              parse_simple_span(bs, span, accession);
              
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
    std::smatch m;
    std::string::const_iterator first = simple_span.begin();
    std::string::const_iterator last = simple_span.end();
    //Rcpp::Rcout << "'" << simple_span << "' => ";
    //Rcpp::Rcout << *first << " ... " << *last << std::endl;
    
    // initialize and set defaults
    std::string start, end;
    std::vector<unsigned int> range(2);
    std::vector<bool> fuzzy(2, false);
    int strand(1);
    bool remote(false);
    std::string type("R");
  
    // test for complement
    if (std::regex_search(first, last, m, COMPL)) {
        // we have matched "complement(m..n)"
        strand = -1;
    }
    
    // test for gap
    if (std::regex_search(first, last, m, GAP)) {
        // we have matched "gap(), gap(X), gap(unkX)"
        type = "G";
        accession = "";
        start = "1";
        // gap length
        std::regex_search(first, last, m, GAPLEN);
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
        std::regex_search(first, last, m, RASL);
        if (m[1].matched) {
            remote = true;
            accession = m[2];  // remote accession without colon
            first = m[4].first;
            last = m[4].second;
        } else {
            first = m[0].first;
            last = m[0].second;
        }
        //Rcpp::Rcout << accn << std::endl;
        //Rcpp::Rcout << "'" << simple_span << "' => ";
        //Rcpp::Rcout << *first << " ... " << *last << std::endl;

        // get type
        if (std::regex_search(first, last, m, BETWEEN_BASES)) {
            type = "B";
        }
  
        // split span
        std::regex_search(first, last, m, BASESPLIT);
        start = m[1];
        if (m[3].matched) {
            end = m[3];
        } else {
            end = start;
        }
        //Rcpp::Rcout << "'" << simple_span << "' => ";
        //Rcpp::Rcout << start << " ... " << end << std::endl;
        fuzzy[0] = std::regex_match(start, m, FUZZY_START);
        fuzzy[1] = std::regex_match(end, m, FUZZY_END);
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

static inline unsigned int extractNumber(const std::string& str) {
    std::string temp;
    std::string::const_iterator it;
    for (it = str.begin(); it != str.end(); it++)
    {
        if (isdigit(*it)) {
            temp += *it;
        }
    }
    return(atoi(temp.c_str()));
}
