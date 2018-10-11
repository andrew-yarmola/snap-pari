#ifndef _fundamental_group_
#define _fundamental_group_

#include "Moebius_transformations.h"
#include "FGWord.h"

/* 
 * Added by oag. This is a semi-private interface
 * to the GroupPresentation, allowing a few extra
 * operations, beyond those provided in SnapPea.h.
 */

struct CyclicWord;
struct GroupPresentation;

void simplify_word(CyclicWord *word);
void init_name(CyclicWord **name);
void set_name(CyclicWord  **name,int i);
void print_word(CyclicWord *name, FILE* out);

void free_cyclic_word(CyclicWord *aCyclicWord);
void append_word(CyclicWord *source, CyclicWord *dest);
void append_inverse(CyclicWord *source, CyclicWord *dest);
void cancel_inverses_word(CyclicWord *word);
std::string string_form(const CyclicWord& cw);
CyclicWord* to_cyclic_word(FGWord const& word);

O31_matrix* fg_get_generators(GroupPresentation* G);
MoebiusTransformation* fg_Moebius_generators(GroupPresentation* G); 
GroupPresentation* read_presentation(FILE* fp);
void save_presentation(FILE* fp, GroupPresentation* G, bool closed=false);
void set_filling_flag(GroupPresentation* G, Boolean flag);

int  add_generator(GroupPresentation* group, const O31_matrix& m);
void add_relation(GroupPresentation* group, CyclicWord* word);
void simplify(GroupPresentation *group);

#endif
