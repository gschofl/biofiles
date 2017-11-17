#include "biofiles.h"

using namespace Rcpp;

// [[Rcpp::export]]
SEXP gbFeature(
    std::vector<std::string> feature,
    std::string accession = "",
    int id = 0 )
{
    // create and assign a 'gbFeature' object
    Rcpp::S4 gb_feature("gbFeature");    
    //gb_feature.slot(".seqinfo") = seqinfo;
    gb_feature.slot(".id") = id;
    parse_gb_feature_table( gb_feature, feature, accession );
    return gb_feature; 
} 

void parse_gb_feature_table(
    Rcpp::S4& gb_feature,
    const std::vector<std::string>& feature_string,
    std::string& accession )
{   
    std::smatch m;
    std::string this_line("");
    std::vector<std::string> merged_lines;
    // merge qualifier values that occupy more than one line into "merged_lines"
    // Define qualifier position: Optional FT + blanks followed by forward slash
    static const std::regex QUAL_POS("^(FT)?\\s+/");
    static const std::regex TRANS("^translation=");
    std::vector<std::string>::const_iterator s_it = std::begin(feature_string);
    for ( ; s_it != std::end(feature_string); s_it++ ) {
        std::regex_search( std::begin(*s_it), std::end(*s_it), m, QUAL_POS ); 
        if ( m[0].matched ) {
            merged_lines.push_back(this_line);
            // If we have an embl file trim FT, whitespace and /
            // If we have a gbk file trim whitespace and /
            this_line = *s_it;
            this_line = ltrim(this_line, "FT");
            this_line = ltrim(this_line, " /");
            this_line = rtrim(this_line, " \n\r\t");
        } else {
            std::string tmp = *s_it;
            tmp = ltrim(tmp, "FT");
            tmp = ltrim(tmp, " /");
            tmp = rtrim(tmp, " \n\r\t");
            if ( std::regex_search( this_line, m, TRANS ) ) {
                // don't add whitespace when joining translations
                this_line += tmp;
            } else {
                this_line += " " + tmp;
            }
        }
        // Push the last qualifier
        if ( s_it == std::end(feature_string) - 1 ) {
            merged_lines.push_back(this_line);
        }
    }
    
    // Get key
    int begin_key = merged_lines[0].find_first_not_of(" ");
    int end_key = merged_lines[0].find_first_of( " ", begin_key );
    std::string key = merged_lines[0].substr( begin_key, end_key - 1 );
    gb_feature.slot("key") = key;
    
    // Get location
    Rcpp::S4 gb_location = Rcpp::S4("gbLocation");
    std::string gb_base_span = merged_lines[0].substr( end_key, merged_lines[0].length() - end_key );
    parse_gb_location( gb_location, gb_base_span, accession );
    gb_feature.slot("location") = gb_location;
    
    merged_lines.erase( std::begin(merged_lines) );
    gb_feature.slot("qualifiers") = parse_gb_qualifiers(merged_lines);
}


Rcpp::CharacterVector parse_gb_qualifiers( 
    const std::vector<std::string>& qualifiers )
{
    int begin_, end_;
    std::vector<std::string> names, values;
    names.reserve( qualifiers.size() );
    values.reserve( qualifiers.size() );
    std::vector<std::string>::const_iterator s_it = std::begin(qualifiers);
    for ( ; s_it != std::end( qualifiers ); s_it++ ) {
        begin_ = s_it->find_first_not_of("=");
        end_ = s_it->find_first_of( "=", begin_ );
        names.push_back( s_it->substr( begin_, end_ ) );
        values.push_back( trim( s_it->substr(end_ + 1, s_it->length() - end_), "\"") );
    }
    Rcpp::CharacterVector Values = Rcpp::CharacterVector( std::begin(values), std::end(values) );
    Values.names() = names;
    return Values; 
}

static inline std::string ltrim(std::string str, std::string ch = " \n\r\t") {
  str.erase( 0, str.find_first_not_of(ch) );
  return str;
}

static inline std::string rtrim(std::string str, std::string ch = " \n\r\t") {
  str.erase( str.find_last_not_of(ch) + 1 );
  return str;
}

static inline std::string trim(std::string str, std::string ch = " \n\r\t") {
  str.erase( 0, str.find_first_not_of(ch) );
  str.erase( str.find_last_not_of(ch) + 1 );
  return str;
}
