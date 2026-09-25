#ifndef LINEFILEUTILITIES_H
#define LINEFILEUTILITIES_H

#include <vector>
#include <string>
#include <cstring>
#include <sstream>
#include <cstdlib>

using namespace std;

// templated function to convert objects to strings
template <typename T>
inline
std::string ToString(const T & value) {
    std::stringstream ss;
    ss << value;
    return ss.str();
}

#ifndef _MSVC_LANG
// tokenize into a list of strings.
inline
void Tokenize(const string &str, vector<string> &elems, const string &delimiter = "\t")
{
    char* tok;
    char cchars [str.size()+1];
    char* cstr = &cchars[0];
    strcpy(cstr, str.c_str());
    tok = strtok(cstr, delimiter.c_str());
    while (tok != NULL) {
        elems.push_back(tok);
        tok = strtok(NULL, delimiter.c_str());
    }
}

// tokenize into a list of integers
inline
void Tokenize(const string &str, vector<int> &elems, const string &delimiter = "\t") {
    char* tok;
    char cchars [str.size()+1];
    char* cstr = &cchars[0];
    strcpy(cstr, str.c_str());
    tok = strtok(cstr, delimiter.c_str());
    while (tok != NULL) {
        elems.push_back(atoi(tok));
        tok = strtok(NULL, delimiter.c_str());
    }
}
#else
inline void Tokenize(const string &str, vector<string> &tokens, const string &delimiter = "\t")
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiter, 0);
    // Find first "non-delimiter".
    string::size_type pos = str.find_first_of(delimiter, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiter, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiter, lastPos);
    }
}
inline void Tokenize(const string &str, vector<int> &tokens, const string &delimiter = "\t")
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiter, 0);
    // Find first "non-delimiter".
    string::size_type pos = str.find_first_of(delimiter, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(atoi(str.substr(lastPos, pos - lastPos).c_str()));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiter, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiter, lastPos);
    }
}
#endif

#endif /* LINEFILEUTILITIES_H */

