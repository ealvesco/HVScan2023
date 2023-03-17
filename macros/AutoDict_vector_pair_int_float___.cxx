#include "vector"
#include "utility"
#ifdef __CINT__ 
#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedefs;
#pragma link C++ class vector<pair<int,float> >+;
#pragma link C++ class vector<pair<int,float> >::*;
#ifdef G__VECTOR_HAS_CLASS_ITERATOR
#pragma link C++ operators vector<pair<int,float> >::iterator;
#pragma link C++ operators vector<pair<int,float> >::const_iterator;
#pragma link C++ operators vector<pair<int,float> >::reverse_iterator;
#endif
#endif
