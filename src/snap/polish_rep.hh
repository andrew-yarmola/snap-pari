#ifndef _polish_rep_hh_
#define _polish_rep_hh_

#include "snappea/SnapPea.h"
#include <vector> 

bool polish(GroupPresentation* G, bool filled, bool report, int limit=20);
bool normalize(GroupPresentation* G, MoebiusTransformation& conj, bool report=false, bool try_nice=false, std::vector<FGWord> const& wl = std::vector<FGWord>()); 

double accuracy(GroupPresentation* G);

void print_linear_errors(GroupPresentation* G);
void testmats(); 

#endif
