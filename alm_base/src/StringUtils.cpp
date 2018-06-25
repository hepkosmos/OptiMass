#include "StringUtils.h"
#include <string>
#include <sstream>
#include <vector>

std::vector<std::string>& OptiMass::split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> OptiMass::split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}
std::vector<std::string> OptiMass::splitKeepParenthesis(const std::string &s){
    std::string delim("() \t");
    std::string whitespace(" \t\n");
    std::string strbuf(s);
    trim(strbuf,delim);
    /*
    bool chk = true;
    while(chk)
        if(s.at(0) == '('){
            s.erase(0,1);
            s.erase(s.length()-1,1);
        }else{
            chk = false;
        }
    }*/

    std::vector<std::string> elems;
    int level = 0;
    std::string buf;
    for(auto it=strbuf.begin();it != strbuf.end(); ++it){
        if( (*it != '(' && *it != ')' && *it != ',' ) || (*it == ',' && level > 0)){
            buf.push_back(*it);
        }else if(*it == ',' && level == 0){
            elems.push_back(trim(buf, whitespace));
            buf.clear();
        }else if(*it == '(' && level == 0){
            level += 1;
        }else if(*it == '(' && level > 0){
            buf.push_back(*it);
            level += 1;
        }else if(*it == ')' && level == 1){
            level -= 1;
        }else if(*it == ')' && level > 1){
            buf += *it;
            level -= 1;
        }
    }
    elems.push_back(trim(buf, whitespace));
    return elems;
}

std::string& OptiMass::ltrim(std::string &s, std::string& delim) {
    std::size_t found = s.find_first_not_of(delim);
    if(found == std::string::npos){
        s.clear();
    }else{
        s.erase(0,found);
    }
    return s;
}
// trim from end
std::string& OptiMass::rtrim(std::string &s, std::string& delim) {
    std::size_t found = s.find_last_not_of(delim);
    if(found != std::string::npos){
        s.erase(found+1);
    }else{
        s.clear();
    }
    return s;
}
// trim from both ends
std::string& OptiMass::trim(std::string &s, std::string& delim) {
    return ltrim(rtrim(s, delim),delim);
}
