#include "gb_regex.h"
#include <boost/algorithm/string_regex.hpp>

using namespace Rcpp;

//' @internal
// [[Rcpp::export]]
SEXP gbFeature(
    std::vector<std::string> feature,
    Rcpp::Environment seqenv,
    int id = 0 )
{
    // create and assign a 'gbFeature' object
    Rcpp::S4 gb_feature("gbFeature");
    std::string accession;
    
    if ( not seqenv.exists("seqinfo") ) {
        Rcpp::S4 seqinfo("Seqinfo");
        seqenv.assign("seqinfo", seqinfo);
        accession = "";
    } else {
        Rcpp::S4 seqinfo = seqenv.get("seqinfo");
        std::string accession = seqinfo.slot("seqnames");
    }
    
    if ( not seqenv.exists("sequence") ) {
      Rcpp::S4 sequence("DNAStringSet");
      seqenv.assign("sequence", sequence);
    }
    
    gb_feature.slot(".seqinfo") = seqenv;
    gb_feature.slot(".id") = id;
    parse_gb_feature_table( gb_feature, feature, accession );
    return gb_feature; 
} 

void parse_gb_feature_table(
    Rcpp::S4& gb_feature,
    const std::vector<std::string>& feature_string,
    std::string& accession )
{   
    boost::smatch m;
    std::string this_line("");
    std::vector<std::string> merged;
    static const std::string WS(" ");
    // merge qualifier values that occupy more than one line
    // Define qualifier position: Optional blanks followed by forward slash
    static const boost::regex QUAL_POS("^(FT)?\\s*/");
    static const boost::regex TRANS("^translation=");
    std::vector<std::string>::const_iterator s_it = feature_string.begin();
    for ( s_it; s_it != feature_string.end(); s_it++ ) {
        boost::regex_search( s_it->begin(), s_it->end(), m, QUAL_POS ); 
        if ( m[0].matched ) {
            merged.push_back( this_line );
            this_line = boost::trim_left_copy_if(boost::trim_copy( *s_it ), boost::is_any_of("/"));
        } else {
            if ( boost::regex_search( this_line, m, TRANS) ) {
                this_line += boost::trim_copy( *s_it ); // don't add whitespace when joining translations
            } else {
                this_line += WS + boost::trim_copy( *s_it );
            }
        }
        // Push the last qualifier
        if ( s_it == feature_string.end() - 1 ) 
            merged.push_back( this_line );
    }
    
    // Get key
    int begin_key = merged[0].find_first_not_of(" ");
    int end_key = merged[0].find_first_of(" ", begin_key);
    std::string key = merged[0].substr(begin_key, end_key - 1);
    gb_feature.slot("key") = key;
    
    // Get location
    Rcpp::S4 gb_location = Rcpp::S4("gbLocation");
    std::string gb_base_span = merged[0].substr(end_key, merged[0].length() - end_key);
    parse_gb_location( gb_location, gb_base_span, accession );
    gb_feature.slot("location") = gb_location;
    
    merged.erase( merged.begin() );
    gb_feature.slot("qualifiers") = parse_gb_qualifiers( merged );
}


Rcpp::CharacterVector parse_gb_qualifiers( 
    const std::vector<std::string>& qualifiers )
{
    int begin_, end_;
    std::vector<std::string> names, values;
    names.reserve( qualifiers.size() );
    values.reserve( qualifiers.size() );
    std::vector<std::string>::const_iterator s_it = qualifiers.begin();
    for ( s_it; s_it != qualifiers.end(); s_it++ ) {
        begin_ = s_it->find_first_not_of("=");
        end_ = s_it->find_first_of("=", begin_);
        names.push_back( s_it->substr(begin_, end_) );
        values.push_back( boost::trim_copy_if( s_it->substr(end_ + 1, s_it->length() - end_), boost::is_any_of("\"") ) );
    }
    Rcpp::CharacterVector Values = Rcpp::CharacterVector( values.begin(), values.end() );
    Values.names() =  names;
    return Values; 
}

