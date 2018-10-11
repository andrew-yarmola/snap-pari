#include "alg_group.hh"
#include "kernel_extras.hh"
#include "pari_code.hh"
#include "polish_rep.hh"
#include "snappea/fundamental_group.h"

using std::cout;
using std::cerr;
using std::endl; 

const pari_printform pp_short = pari_printform(1,1,'g',16,0); 

void alg_group::clear() 
{
  exact_gens = pZERO; 
  hilbert[0] = pZERO;
  hilbert[1] = pZERO; 
  _integer_traces = -1; 

  trace.clear(); it.clear(); groupf.clear(); 

  itgens[0] = pZERO;
  itgens[1] = pZERO; 
  tracegens[0] = pZERO;
  tracegens[1] = pZERO; 
  G = 0;
}

void alg_group::set_group(GroupPresentation* g)
{
  clear();
  G = g; 
}

static int count_bits(unsigned long bitvec)
{
  int n = 0; 
  while (bitvec) {
    if (bitvec & 1L) n++;
    bitvec >>= 1;
  }
  return n;
}

static void next_subset(unsigned long& bitvec, int max)
{
  bitvec += (count_bits(bitvec) < max) ? 1 : (bitvec & ~(bitvec - 1));
}

static int num_itgens(int n)
{
  /* check we do not have too many generators */ 
  int n_itgens = 0; 
  unsigned long stop_bit = 1L << n; 
  if (!stop_bit) {
    warn("Too many generators to compute all traces\n");
    return 0;
  }

  /* count how many products of generators we have to consider */ 
  unsigned long bitvec = 1L; 
  while (bitvec < stop_bit) {
    n_itgens++; 
    next_subset(bitvec, 3);
  }

  return n_itgens;
}

static void get_trace_names(int n, vector<string>& tracegen_names, vector<string>& itgen_names)
{
  int n_itgens = num_itgens(n); 
  if (!n_itgens) return; 

  itgen_names.resize(n_itgens); 
  tracegen_names.resize(n_itgens); 

  /* get names of all products */ 
  char element[100], element2[100], s[10]; 
  unsigned long stop_bit = 1L << n; 
  unsigned long copy, bitvec;
  int i, j; 

  for (i = 0, bitvec = 1L; bitvec < stop_bit; i++, next_subset(bitvec,3)) {
    strcpy(element, "tr(");
    strcpy(element2, "tr(");
    for (j = 0, copy = bitvec; copy; j++, copy >>= 1) {
      if (copy & 1) {
	sprintf(s, "%c^2", j+'a');
	strcat(element, s);
	sprintf(s, "%c", j+'a');
	strcat(element2, s);
      }
    }
    strcat(element, ")"); 
    strcat(element2, ")"); 
    itgen_names[i] = element;
    tracegen_names[i] = element2; 
  }
}

static void print_list(vector<string> const& names, const pari& values)
{
  int j, num = length(values);
  for(j = 0; j < num; j++) {
    printf("%s = ", names[j].c_str()); 
    values[j].print(); 
  }
}

void alg_group::print() const
{
  printf("Fundamental group generators\n"); 
  print_group_generators(); 

  vector<string> itgen_names, tracegen_names; 
  get_trace_names(fg_get_num_generators(G), itgen_names, tracegen_names);

  int i = (tracegens[1].type()==1) ? 0 : 1; 
  if (tracegens[i].type()!=1) {
    printf("\nTraces\n"); 
    print_list(tracegen_names,tracegens[i]); 
  }

  i = (itgens[1].type()==1) ? 0 : 1; 
  if (itgens[i].type()!=1) {
    printf("\nInvariant trace field generators\n"); 
    print_list(itgen_names,itgens[i]); 
  }
}

void alg_group::print_fields(int print_fields_long) const
{
  if (it.is_set()) { 
    printf("Invariant trace field\n"); 
    it.print(print_fields_long); 
    printf("\n"); 
  }

  if (trace.is_set()) { 
    printf("Trace field\n"); 
    trace.print(print_fields_long); 
    printf("\n"); 
  }

  if (groupf.is_set()) {
    printf("Group coefficient field\n"); 
    groupf.print(print_fields_long); 
    printf("\n"); 
    if (trace.is_set() && 
	trace.degree() != groupf.degree()) {
      printf("  Trace field generator: "); 
      lift(groupf.exact_value(trace.root())).print(); 
    }
    if (it.is_set() &&
	(!trace.is_set() || 
	 trace.degree() != it.degree())) {
      printf("  Invariant trace field generator: "); 
      lift(groupf.exact_value(it.root())).print(); 
    }
  }
}
  
void alg_group::print_group_generators() const
{
  if (exact_gens.type()==1) {
    int i, n_gens = fg_get_num_generators(G);
    for (i=0; i<n_gens; ++i) 
      cout << fg_Moebius_gen(G,i) << endl; 
    return;
  }

  int i, n = length(exact_gens); 
  for (i = 0; i < n; i++) {
    printf("%c = ", i+'a'); 
    lift(exact_gens[i]).print();
  }
}

// Needs accurate generators for the fundamental group, generators[..]. 
// Computes the numeric field generators, itgens[0] and tracegens[0]. 
// Also calls set_trace_names() to set up n_itgens, tracegen_names and itgen_names. 

void alg_group::set_tracefield_gens()
{
  int ng = fg_get_num_generators(G);
  int n_itgens = num_itgens(ng); 
  if (!n_itgens) return; 

  itgens[0] = rvector(n_itgens); 
  tracegens[0] = rvector(n_itgens); 

  /* compute products of generators */ 
  pari acm; 
  pari m = matrix(2,2), mm = matrix(2,2);
  unsigned long stop_bit = 1L << ng; 
  unsigned long bitvec, copy; 
  int i, j, first; 

  for (i = 0, bitvec = 1L; bitvec < stop_bit; i++, next_subset(bitvec,3)) {
    first = 1; 
    for (j = 0, copy = bitvec; copy; j++, copy >>= 1) {
      if (copy & 1) {
	acm = fg_Moebius_gen(G,j).acc_matrix; 
	if (!first) {
	  m *= acm;
	  mm *= gsqr(acm); 
	} else {
	  first = 0; 
	  m = acm; 
	  mm = gsqr(acm);
	}
      }
    }
    itgens[0][i] = mm[0][0] + mm[1][1];
    tracegens[0][i] = m[0][0] + m[1][1];
  }
}

void alg_group::update_precision()
{
  if (!G) return; 
  polish(G, false, false);
  if (itgens[0]!=pZERO) set_tracefield_gens(); 

  trace.update_precision(); 
  it.update_precision();
  groupf.update_precision();
}

void alg_group::set_exact_gens()
{
  if (!groupf.is_set()) return; 

  int i, j, k, n_gens;

  n_gens = fg_get_num_generators(G); 
  exact_gens = rvector(n_gens); 
  for (i=0; i < n_gens; i++) {
    exact_gens[i] = matrix(2,2); 
      
    for (j=0; j<2; j++)
      for (k=0; k<2; k++)
	exact_gens[i][j][k] = 
	  groupf.exact_value(fg_Moebius_gen(G,i).acc_matrix[j][k]);
  }
}

static int get_groupfield_gens(GroupPresentation* G, vector<string>& names, pari& gens)
{
  MoebiusTransformation* gen = fg_Moebius_generators(G); 
  int n = fg_get_num_generators(G); 
  char s[50]; 
  names.resize(3*n); 
  gens = rvector(3*n); 
  pari acm; 
  int i; 
  for (i = 0; i < n; i++) {
    acm = gen[i].acc_matrix;
    gens[3*i] = acm[0][0];
    gens[3*i+1] = acm[0][1];
    gens[3*i+2] = acm[1][1];
    sprintf(s, "%c[0][0]", i+'a');
    names[3*i] = string(s);
    sprintf(s, "%c[0][1]", i+'a');
    names[3*i+1] = string(s);
    sprintf(s, "%c[1][1]", i+'a');
    names[3*i+2] = string(s); 
  }
  return 3*n; 
}

int alg_group::integer_traces() const
{
  if (_integer_traces==-1) 
    return has_integer_traces(); 
  return _integer_traces; 
}


int alg_group::has_integer_traces() const
{
  // Check for integer traces.
  int i; 
  int int_tr=1;
  int n_itgens = length(tracegens[0]);
  if (trace.is_set()) { 
    for (i=0; int_tr && i<n_itgens; i++)
      if (!is_integral(tracegens[1][i])) int_tr = 0; 
  } else if (it.is_set()) { 
    pari x_trsq; 
    for (i=0; int_tr && i<n_itgens; i++) {
      if (!it.contains(gsqr(tracegens[0][i]), x_trsq)) {
	warn("unable to find the exact value of a square of a trace\n"); 
	int_tr = -1; break; }
      if (!is_integral(x_trsq)) int_tr = 0;
    }
  } else return -1; 
  return int_tr; 
}


/* "which" is bitwise AND of 
   0x1 for invariant trace field, 
   0x2 for trace field
   0x8 for a group coefficient field

   "report" is bitwise AND of 
   0x1 for trace of the field computation
   0x2 for printout of the answers 
 */ 

void alg_group::compute_fields(int cc_poly, int report, unsigned int which)
{
  vector<string> names;
  pari gens, exact_gf_gens, exact_traces, exact_itgens; 
  int n; 
  int talk = report & 1; 

  /* get generators for the fields */ 
  set_tracefield_gens();

  if (which & 0x1) {

    if (talk) printf("Computing invariant trace field\n"); 
    if (trace.is_set())
      it.generated_by(itgens[0], exact_itgens, cc_poly, &trace, 1, talk);
    else if (groupf.is_set())
      it.generated_by(itgens[0], exact_itgens, cc_poly, &groupf, 1, talk);
    else
      it.generated_by(itgens[0], exact_itgens, cc_poly, talk);

    if (it.is_set()) {
      itgens[1] = exact_itgens;
      if (report & 2) { printf("Invariant trace field: "); it.print(0); printf("\n"); }
    } else {
      if (report & 2) { printf("Invariant trace field not found.\n"); }
    }
  }

  if (which & 0x2) {

    if (talk) printf("Computing trace field\n"); 

    if (it.is_set())
      trace.generated_by(tracegens[0], exact_traces, cc_poly, &it, 2, talk);
    else if (groupf.is_set())
      trace.generated_by(tracegens[0], exact_traces, cc_poly, &groupf, 1, talk);    else
      trace.generated_by(tracegens[0], exact_traces, cc_poly, talk);

    if (trace.is_set()) {
      tracegens[1] = exact_traces;
      if (report & 2) { printf("Trace field: "); trace.print(0); printf("\n"); }
    } else {
      if (report & 2) { printf("Trace field not found.\n"); }
    }
  }


  if (which & 0x8) {

    if (talk) printf("Computing group coeff field\n"); 
    n = get_groupfield_gens(G, names, gens);

    if (trace.is_set())
      groupf.generated_by(gens, exact_gf_gens, cc_poly>1, &trace, 2, talk);
    else if (it.is_set())
      groupf.generated_by(gens, exact_gf_gens, cc_poly>1, &it, 2, talk);
    else 
      groupf.generated_by(gens, exact_gf_gens, cc_poly>1, talk);

    if (groupf.is_set()) {
      set_exact_gens(); 
      if (report & 2) { printf("Group coeff field: "); groupf.print(0); printf("\n"); }
    } else {
      if (report & 2) printf("Group coeff field not found.\n"); 
    }
  }

  // Check for integer traces if still unknown. 
  if (_integer_traces==-1 && (it.is_set() || trace.is_set())) {
    _integer_traces = has_integer_traces(); 
  }

}

/* must have invariant trace field and hilbert symbol before calling this */ 

pari alg_group::real_ramification() const 
{
  int i,j, r1 = (it.p_nf())[1][0].int_value();
  int nrr_places = 0; 
  int *which_places = new int[r1]; 

  /* see which real places ramify and count them */ 
  for (i = 1; i <= r1; i++) {
    if (it.numeric_value(hilbert[0], i) < pZERO && 
	it.numeric_value(hilbert[1], i) < pZERO) {
      nrr_places++;
      which_places[i-1] = 1; 
    } else 
      which_places[i-1] = 0; 
  }

  /* copy them out into a pari row vector */ 
  j = 0; 
  pari real_ram = rvector(nrr_places);
  for (i = 1; i <= r1; i++)
    if (which_places[i-1]) real_ram[j++] = integer(i); 

  delete [] which_places; 

  return real_ram; 
}

int alg_group::arithmetic() const
{
  // Need the hilbert symbol for this function. 
  if (hilbert[0].type() == t_INT) return -1; 

  /* more than one complex place? */ 
  if ((it.p_nf())[1][1] > pONE) {
    return 0; 
  }

  /* any unramified real places? */ 
  int i, r1 = (it.p_nf())[1][0].int_value();
  pari tmp = rvector(2); 
  for (i = 1; i <= r1; i++) {
    tmp[0] = it.numeric_value(hilbert[0], i);
    tmp[1] = it.numeric_value(hilbert[1], i);

    if (tmp[0] > pZERO || tmp[1] > pZERO) {
      return 0;
    }
  }

  // Now it just depends on whether we have integer traces. 
  return integer_traces(); 
}

// exact_gens is expected to be a vector of 2x2 complex matrices 
// in SL(2,C). 

static pari word_to_Moebius(pari exact_gens, FGWord const& word)
{
  pari mt = idmat(2); 
  int i, n = word.length();

  int g; 
  for (i=0; i<n; i++) {
    g = word[i] > 0 ? word[i] : -word[i]; 
    if (g > length(exact_gens) || g==0) {
      warn("non existent generator in word passed to word_to_Moebius\n"); 
      continue; 
    }
    if (word[i] > 0) 
      mt *= exact_gens[g-1];
    else 
      mt *= ginv(exact_gens[g-1]);
  }
  return mt; 
}

static pari word_to_Moebius(pari exact_gens, const int* word)
{
  pari mt = idmat(2); 
  const int* it = word; 
  int g; 
  for (; *it != 0; it++) {
    g = *it > 0 ? *it : -*it; 
    if (g > length(exact_gens) || g==0) {
      warn("non existent generator in word passed to word_to_Moebius\n"); 
      continue; 
    }
    if (*it > 0) 
      mt *= exact_gens[g-1];
    else 
      mt *= ginv(exact_gens[g-1]);
  }
  return mt; 
}

static bool in_interval_minus2_2(pari const& x)
{
  return gabs(greal(x)) < real(2.01) && 
    gabs(gimag(x)) < real(0.01); 
}

bool alg_group::compute_hilbert_symbol()
{
  FGWord a, b; 
  int n_gens = fg_get_num_generators(G);

  a.next_word(n_gens); 
  b.next_word(n_gens); 
  b.next_word(n_gens); 

  MoebiusTransformation m0, m1;
  pari sq0, sq1, comm;

  // Check accurate values are set. 
  m0 = word_to_Moebius(G, a);
  if (m0.acc_matrix.type()==1) {
    cout << "Problem computing accurate group elements for hilbert symbol.\n";
    return false;
  }

  bool found=false; 
  while (b.length() < 4) {
    while (!(a == b)) {

      m0 = word_to_Moebius(G, a);
      m1 = word_to_Moebius(G, b);

      sq0 = gsqr(m0.acc_matrix);
      sq1 = gsqr(m1.acc_matrix);

      comm = sq0 * sq1 * ginv(sq0) * ginv(sq1); 

      if (!in_interval_minus2_2(gtrace(sq0)) && !in_interval_minus2_2(gtrace(comm))) {
	found = true; 
	break;
      }
      a.next_word(n_gens); 
    }
    if (found) break; 

    b.next_word(n_gens); 
    a = FGWord();
    a.next_word(n_gens); 
  }

  if (!found) return false; 

  pari h0 = gsqr(gtrace(sq0)) - pFOUR;
  pari h1 = gtrace(comm) - pTWO; 
  return it.contains(h0, hilbert[0]) && it.contains(h1, hilbert[1]);
}

void alg_group::print_hilbert_symbol(int show) const
{
  printf("("); 
  lift(hilbert[0]).print(no_newline);
  printf(", "); 
  lift(hilbert[1]).print(no_newline);
  printf(")\n"); 

  if (show & 0x1) { // Show factorizations. 
    printf("Factorization of a:\n"); 
    idealfactor(it.p_nf(), hilbert[0]).print(); 
    if (hilbert[1]!=pZERO) {
      printf("Factorization of b:\n"); 
      idealfactor(it.p_nf(), hilbert[1]).print();
    }
    printf("Factorization of 2:\n"); 
    idealfactor(it.p_nf(), pTWO).print();
  }
  if (show & 0x2) { // Show real places. 
    int r1 = (it.p_nf())[1][0].int_value(), i;
    pari tmp = rvector(2); 
    printf("Real embeddings of [a,b]:\n"); 
    for (i = 1; i <= r1; i++) {
      tmp[0] = it.numeric_value(hilbert[0], i);
      tmp[1] = it.numeric_value(hilbert[1], i);
      printf("%d: ",i); tmp.print(pp_short);
    }
  }
}

void print_relation(const int* l, FILE *fp)
{
  if (!fp) fp = stdout;

  const int* it = l;
  int c; 
  for (; *it; it++) {
    c = (*it > 0) ? 'a' + *it - 1 : 'A' - *it - 1;
    fprintf(fp,"%c", c);
  }
}

pari alg_group::alg_group_element(const FGWord& word) const
{
  if (exact_gens.type()==1) return pZERO; 
  return word_to_Moebius(exact_gens, word); 
}

int alg_group::verify() const
{
  if (exact_gens.type()==1) return 0; 

  int i, n_rels = fg_get_num_relations(G);
  int *word;
  int OK = 1; 
  pari mx, idmx = idmat(2); 

  printf("Fundamental Group Relations\n");
  for (i=0; i<n_rels; i++) {

    word = fg_get_relation(G, i); 
    print_relation(word, stdout); 

    printf(" -> ");

    mx = lift(word_to_Moebius(exact_gens, word)); 
    mx.print(); 

    if (mx != idmx && mx != -idmx) OK = 0; 
    fg_free_relation(word); 
  }

  return OK; 
}


bool alg_group::imquad_integer_itgens() const
{
  if (itgens[0]==pZERO) {
    cerr << "this test requires accurate invariant trace field generators to be set\n";
    return true; // this is the less conclusive of the two possible results. 
  }

  int i, n_itgens = length(itgens[0]); 
  for (i=0; i<n_itgens; i++) {
    if (!is_imag_quad_integer(itgens[0][i])) return false; 
  }
  return true; 
}

pari alg_group::exact_trace(FGWord const& w) const
{
  pari tr = real(0);

  if (!trace.is_set()) {
    cout << "This function requires the trace field to be set.\n";
    return tr;
  }

  pari mx = word_to_Moebius(G, w).acc_matrix;
  if (mx.type()==t_INT) {
    cout << "Problem evaluating this element.\n";
    return tr; 
  }

  tr = mx[0][0] + mx[1][1]; 
  return trace.exact_value(tr);
}

pari alg_group::finite_ramification(int max_time) const
{
  return ::ramification(hilbert_symbol(), itfield().p_nf(), 
			integer_traces(), 0, max_time);
}

pari alg_group::ramification(int max_time) const
{
  pari fr = finite_ramification(max_time);
  pari ram = concat(real_ramification(), fr);

  // There should be an even number of places.
  // If there is at most one place where the test was
  // aborted we can determine whether the place should
  // be included or not. 

  int i, n = length(ram), num_unknown=0;
  for (i=0; i<n; ++i) {
    if (ram[i].type()!=t_INT && ram[i][0] < pZERO) ++num_unknown;
  }

  if (n%2==1 && num_unknown==0) {
    cout << "Ramification returning an odd number of places!!\n";
  }
  if (num_unknown==0) return ram;
  if (num_unknown >1) {
    cout << "Can't fix up ramification, too many places unknown\n";
    return ram;
  }

  // Fix the unknown place.
  if (n%2==1) { // n is odd so get rid of the unknown place
    pari rfix = rvector(n-1);
    int j=0;
    for (i=0; i<n; ++i) {
      if (ram[i].type()!=t_INT && ram[i][0] < pZERO) continue;
      rfix[j] = ram[i]; ++j;
    }
    ram = rfix;
  } else { // n is even so keep the unknown place
    for (i=0; i<n; ++i)
      if (ram[i].type()!=t_INT && ram[i][0] < pZERO) break;
    ram[i][0] = -ram[i][0];
  }
  return ram; 
}

