#ifndef ALM_BASE_SRC_STRINGUTILS_H_
#define ALM_BASE_SRC_STRINGUTILS_H_

#include <string>
#include <vector>

namespace OptiMass {
std::vector<std::string> &split(const std::string &s, char delim,
                                std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);
std::vector<std::string> splitKeepParenthesis(const std::string &s);

// trim from start
std::string &ltrim(std::string &s, std::string &delim);
// trim from end
std::string &rtrim(std::string &s, std::string &delim);
// trim from both ends
std::string &trim(std::string &s, std::string &delim);
}  // end namespace OptiMass

#endif  // ALM_BASE_SRC_STRINGUTILS_H_
