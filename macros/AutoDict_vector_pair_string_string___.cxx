#include "vector"
#include "string"
#ifdef __CINT__ 
#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedefs;
#pragma link C++ class vector<pair<string,string> >+;
#pragma link C++ class vector<pair<string,string> >::*;
#ifdef G__VECTOR_HAS_CLASS_ITERATOR
#pragma link C++ operators vector<pair<string,string> >::iterator;
#pragma link C++ operators vector<pair<string,string> >::const_iterator;
#pragma link C++ operators vector<pair<string,string> >::reverse_iterator;
#endif
#endif
