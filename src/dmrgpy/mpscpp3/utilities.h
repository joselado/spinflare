#include <string>
bool compare_string(auto a, auto b) {
   std::string a2 ;
   std::string b2 ;
   a2 += a ;
   b2 += b ;
   return a2.compare(b2)==0 ;
};
