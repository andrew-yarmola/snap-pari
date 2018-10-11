#include "polish_rep.hh"
#include "printable.hh"
#include "snappea/fundamental_group.h"
#include "snappea/kernel_typedefs.h"
#include <ctype.h>
#include "p_VS.hh"
#include "to_pari.hh"

#undef normalize

using std::vector; 
using std::cout;
using std::endl; 

inline static pari& Gen(GroupPresentation* G, int i)
{
  return fg_Moebius_generators(G)[i].acc_matrix;
}

static MoebiusTransformation evaluate(GroupPresentation* G, int* rel, int len=10000)
{
  int i, g;
  bool ginv; 
 
  MoebiusTransformation mt(One,1); 
  mt.acc_matrix = p_lisexpr("[1,0;0,1]"); 
  for (i=0; i<len; ++i) {
    if (!rel[i]) break; 
    ginv = rel[i] < 0; 
    g  = ginv ? (-1-rel[i]) : (rel[i]-1);
    mt*= ginv ? 
      inverse(fg_Moebius_generators(G)[g]) : 
      fg_Moebius_generators(G)[g];
  }

  return mt; 
}

static int get_length(int* rel)
{
  int i=0; 
  while (rel[i]) ++i; 
  return i; 
}

// This is dm/da for k=0,..., dm/dc for k=2
// where m = [a,b;c,d]   with d = (1+bc)/a. 
// or    m = [d,-b;-c,a] if minv is true. 

static MoebiusTransformation dmdk(MoebiusTransformation const& mt, bool minv, int k, bool use_a, bool& holo)
{
  pari const& m(mt.acc_matrix); 
  pari dm = matrix(2,2);
  int i = minv?1:0;
  if (!use_a) i = 1-i; 
  int a = use_a?0:1; 
  if (!k) {
    dm[i][i] = pONE; 
    dm[1-i][1-i] = -m[1-a][1-a]/m[a][a];
  } else {
    dm[k%2][k/2] = (minv?-pONE:pONE);
    dm[1-i][1-i] = m[1-k%2][1-k/2]/m[a][a];
  }

  holo = (!minv || mt.parity==orientation_preserving);
  if (!holo) dm[1-i][1-i] = gconj(dm[1-i][1-i]);

  MoebiusTransformation mdm;
  mdm.acc_matrix = dm; 
  mdm.parity = mt.parity; 
  return mdm; 
}

static void set_R_DR_from_rel(GroupPresentation* G, int* rel, int row, pari& R, pari& DR, vector<bool> const& use_a, bool trace_only) 
{
  MoebiusTransformation mt, ml, mr;
  pari m; 
  int g, rlen, j, k, o;
  int dzbar_cols=3*fg_get_num_generators(G); 
  bool ginv, holo; 

  // evaluate R
  m = evaluate(G, rel).acc_matrix;

  bool pos_trace = greal(m[0][0] + m[1][1]) > pZERO;
  
  if (trace_only) {
    R[row]   = m[0][0] + m[1][1] - (pos_trace ? pTWO : -pTWO); 
  } else {
    R[row]   = m[0][0] - (pos_trace ? pONE : -pONE);
    R[row+1] = m[1][0]; // col 1, row 0 = b entry. 
    R[row+2] = m[0][1]; // c entry
  }

  // evaluate DR
  rlen = get_length(rel);
  for (j=0; j<rlen; ++j) {
    ginv = rel[j] < 0; 
    g = (ginv) ? (-1-rel[j]) : (rel[j]-1);
    
    ml = evaluate(G, rel, j );
    mr = evaluate(G, rel+j+1);

    // find dependency on gen g: a, b, c
    for (k=0; k<3; k++) {
      mt = ml * dmdk(fg_Moebius_generators(G)[g],ginv,k,use_a[g],holo) * mr; 
      m = mt.acc_matrix; 

      o = (ml.parity == holo) ? 0 : dzbar_cols; 

      if (trace_only) {
	DR[3*g+k+o][row  ] += (m[0][0] + m[1][1]); 
      } else {
	DR[3*g+k+o][row  ] += m[0][0];
	DR[3*g+k+o][row+1] += m[1][0];
	DR[3*g+k+o][row+2] += m[0][1];
      }
    }
  }
}

static bool or_pres(GroupPresentation* G)
{
  int i, n=fg_get_num_generators(G);
  for (i=0; i<n; ++i)
    if (fg_Moebius_generators(G)[i].parity==orientation_reversing) return false; 
  return true; 
}

bool compute_R_DR(GroupPresentation* G, pari& R, pari& DR, vector<bool>& use_a)
{
  int n_gens = fg_get_num_generators(G);
  int n_rels = fg_get_num_relations(G); 
  int n_cusps= fg_get_num_cusps(G);
  bool opres = or_pres(G); 

  if (fg_integer_fillings(G)) n_cusps = 0; 

  // set sizes and initialize to zero
  int n_rows = 3*n_rels + 2*n_cusps; 
  R  = cvector(n_rows);
  DR = matrix((opres?3:6)*n_gens, n_rows); 

  int i; 

  // figure out which coordinates to use
  use_a.resize(n_gens); 
  for (i=0; i<n_gens; ++i) 
    use_a[i] = (abs(Gen(G,i)[0][0]) > abs(Gen(G,i)[1][1])); 

  int *rel; 
  for (i=0; i<n_rels; ++i) {
    rel = fg_get_relation(G, i);
    set_R_DR_from_rel(G, rel, 3*i, R, DR, use_a, false);
    fg_free_relation(rel); 
  }

  int rlen; 
  int row0 = 3*n_rels; 
  for (i=0; i<n_cusps; ++i) {
    rel = fg_get_meridian(G, i);
    set_R_DR_from_rel(G, rel, row0+2*i, R, DR, use_a, true);
    fg_free_relation(rel); 

    rel = fg_get_longitude(G, i);
    set_R_DR_from_rel(G, rel, row0+2*i+1, R, DR, use_a, true);
    fg_free_relation(rel); 
  }

  return true; 
}

static inline int real_prec(pari const& r)
{
  return int(double(length(abs(r))-2)* pariK);
}

static void fix_det(pari m, int a=-1)
{
  if (a < 0) a = (abs(m[0][0]) > abs(m[1][1])) ? 0:1; 
  m[1-a][1-a] = (pONE + m[1][0]*m[0][1])/m[a][a];
}

static void set_prec(MoebiusTransformation& m, int ndig)
{
  m.acc_matrix = gprec(to_pari(m)[0],ndig);
  fix_det(m.acc_matrix); 
}

static void set_approx(MoebiusTransformation& m)
{
  m.matrix[0][0] = to_complex(m.acc_matrix[0][0]); 
  m.matrix[1][0] = to_complex(m.acc_matrix[0][1]); 
  m.matrix[0][1] = to_complex(m.acc_matrix[1][0]); 
  m.matrix[1][1] = to_complex(m.acc_matrix[1][1]); 
}

static void adjust_matrices(GroupPresentation* G, pari const& step, vector<bool> const& use_a)
{
  int i, n_gens = fg_get_num_generators(G);
  int a=0; 
  for (i=0; i<n_gens; ++i) {
    a = use_a[i] ? 0:1;
    Gen(G,i)[a][a] += step[3*i];
    Gen(G,i)[1][0] += step[3*i+1]; // db
    Gen(G,i)[0][1] += step[3*i+2]; // dc
    fix_det(Gen(G,i),a); 
  }
}

// We're assuming matrix DR = [A|B] 
// which is to be evaluated on z as
// [A|B] * [     z ]
//         [conj(z)]
// To get the equivalent real matrix we
// replace entries a = x+i*y, b = u+i*v
// with block submatrices:
// [ x  -y]   [ u   v]
// [ y   x] + [ v  -u]

static void make_DR_real(pari& DR)
{
  int nc = length(DR)/2;
  int nr = length(DR[0]);
  int r, c; 

  pari a,b,RDR = matrix(2*nc, 2*nr);
  for (c=0; c<nc; ++c) {
    for (r=0; r<nr; ++r) {
      a = DR[c   ][r];
      b = DR[c+nc][r];
      RDR[2*c  ][2*r  ] = greal(a)+greal(b);
      RDR[2*c+1][2*r  ] =-gimag(a)+gimag(b);
      RDR[2*c  ][2*r+1] = gimag(a)+gimag(b);
      RDR[2*c+1][2*r+1] = greal(a)-greal(b);
    }
  }
  DR = RDR;
}

static void make_real(pari& V)
{
  int i, n = length(V);
  pari RV = (V.type()==t_VEC) ? rvector(2*n) : cvector(2*n);
  for (i=0; i<n; ++i) {
    RV[2*i  ] = greal(V[i]);
    RV[2*i+1] = gimag(V[i]); 
  }
  V = RV; 
}

static void make_complex(pari& V)
{
  int i, n = length(V)/2;
  pari CV = (V.type()==t_VEC) ? rvector(n) : cvector(n);
  for (i=0; i<n; ++i) {
    CV[i] = complex(V[2*i],V[2*i+1]); 
  }
  V = CV; 
}


// newton_step returns true while the iteration should
// continue. 
   
static bool newton_step(GroupPresentation* G, bool& found, bool report)
{
  found = false; 

  vector<bool> use_a; 

  /* compute error, derivative wrt. generators */
  pari R, DR;
  compute_R_DR(G, R, DR, use_a);

  // how far from solution? 
  pari ersz = vecmax(abs(R)); 
  if (report) ersz.print();

  bool op = or_pres(G); 

  if (!op) {
    make_real(R); 
    make_DR_real(DR); 
  }

  // discard rows which are (approximately) lin dep on others. 
  pari eps = pow(real(10),-6); 
  pari mask = li_col_mask(gtrans(DR), eps); // get row mask
  pari all_cols = pow(pTWO,integer(length(DR))) - pONE; 
  if (mask == pZERO) { 
    printf("vanishing jacobian in polish_rep\n");
    return false; 
  }
  DR = matextract(DR, mask, all_cols); 
  R  = extract(R, mask); 

  // newton step for under-determined system
  pari DRT = gconj(gtrans(DR)); 
  pari step = -DRT*gauss(DR*DRT, R);

  if (!op) 
    make_complex(step); 

  adjust_matrices(G, step, use_a);

  found = (ersz < pow(real(10),-digits)); 
  return !found; 
}

bool polish(GroupPresentation* G, bool filled, bool report, int limit)
{
  // Use this flag to tell compute_R_DR to use traces too. 
  if (!filled) set_filling_flag(G,FALSE);

  // Boost the precision first. 
  int i, n=fg_get_num_generators(G);
  for (i=0; i<n; ++i) set_prec(fg_Moebius_generators(G)[i],digits+18);

  // Polish the accurate matrices. 
  bool found; 
  i=0; 
  while (newton_step(G,found,report) && ++i < limit);
  if (!filled) set_filling_flag(G,TRUE);
  if (!found) {
    if (report) cout << "Newton's method failed to converge\n"; 
    return false;
  }

  // Set the machine precision parts too. 
  for (i=0; i<n; ++i) {
    set_approx(fg_Moebius_generators(G)[i]); 
    fg_get_generators(G)[i] = O31_matrix(fg_Moebius_generators(G)[i]); 
    // cout << fg_Moebius_generators(G)[i] << endl;
  }
  return found; 
}


void testmats()
{
#if 0
  pari f = p_lisexpr("[1,0,.5,0]~");
  pari mink = diagonal(p_lisexpr("[-1,1,1,1]"));
  pari ref = reflection(f);
  cout << "Reflection\n";
  ref.print(); 
  cout << "Reflection^2\n";
  (ref*ref).print(); 
  pari dev = ref * mink * gtrans(ref) - mink; 
  cout << "Deviation\n";
  dev.print(); 
  return; 
#endif

  // Test of above functions.
  pari T = p_lisexpr("[1+2*I,2-I,0,0;2+3*I,0,0,0]");
  cout << "Matrix   "; T.print(); 
  pari RT=T;
  make_DR_real(RT);
  cout << "R-matrix "; RT.print(); 

  pari Z = p_lisexpr("[-2+3*I,-I]~");
  pari RZ=Z;
  make_real(RZ); 
  cout << "Vector   "; Z.print(); 
  cout << "R-vector "; RZ.print(); 

  pari P,RP,CP;
  P  = T*concat(Z,gconj(Z));
  RP = RT*RZ;
  CP = RP; make_complex(CP); 
  cout << "Product   "; P.print();
  cout << "R-product "; RP.print(); 
  cout << "C-product "; CP.print(); 

  pari RS = gauss(RT,RZ);
  cout << "R-solution  "; RS.print();
  cout << "Check "; (RT*RS-RZ).print(); 
  pari CS=RS; make_complex(CS); 
  cout << "C-solution  "; CS.print(); 
  cout << "Check "; (T*concat(CS,gconj(CS)) - Z).print(); 
}

double accuracy(GroupPresentation* G)
{
  pari R, DR;
  vector<bool> use_a; 
  compute_R_DR(G, R, DR, use_a);
  return vecmax(abs(R)).double_value(); 
}

void print_linear_errors(GroupPresentation* G)
{

  // Set the accurate matrices. 
  int i, n=fg_get_num_generators(G);
  for (i=0; i<n; ++i) set_prec(fg_Moebius_generators(G)[i],digits);

  pari R0, R, DR, DR2, LR;
  vector<bool> use_a, tmp_use_a; 
  compute_R_DR(G, R0, DR, use_a);

  cout << "R0\n";
  R0.print(); 

  bool op = or_pres(G); 
  pari dm = cvector(3*n);
  pari eps = pow(real(10),integer(-7));
  pari evec = rvector(op?1:2); 

  for (i=0; i<3*n; ++i) {

    dm[i] = eps; 
    adjust_matrices(G, dm, use_a);

    compute_R_DR(G, R, DR2, tmp_use_a); // just need R. 
    
    LR = (op ? DR*dm : DR*concat(dm, gconj(dm)));
    evec[0] = vecmax(abs(R - R0 - LR))/eps;
    
    adjust_matrices(G,-dm, use_a);

    if (!op) {
      dm[i] = eps*pI;

      adjust_matrices(G, dm, use_a);

      compute_R_DR(G, R, DR2, tmp_use_a); 
    
      LR = DR*concat(dm, gconj(dm));
      evec[1] = vecmax(abs(R - R0 - LR))/eps;
    
      adjust_matrices(G,-dm, use_a);
    }

    evec.print(); 
    dm[i] = pZERO;
  }
}

static bool small(pari const& z)
{
  return gabs(z) < pEPSILON;
}

static bool big(pari const& z)
{
  return gabs(z) > pow(real(10),int(0.75*digits));
}

static bool same_point(pari const& a, pari const& b)
{
  if (gabs(a) < pHALF) {
    if (gabs(b) > pTWO) return false; 
    return small(b - a); 
  }
  return small(ginv(a) - ginv(b));
}

static pari fixed_points(pari const& m)
{
  pari a,b,c,d, diff, discr, tr;

  a = m[0][0]; b = m[1][0]; 
  c = m[0][1]; d = m[1][1]; 

  diff = a - d; 

  pari fix = rvector(2);

  if (small(c)) {      /* c = 0 */
    if (small(diff)) { /* a - d = 0 */
      if (small(b))    /* b = 0, implies identity */
	return p_lisexpr("[0,1,I]");
      fix = rvector(1); /* +/-[1,b;0,1] */ 
      fix[0] = complex(pHUGE,pZERO); 
      return fix; 
    } 

    /* [a,b;0,d] */
    fix[0] = b/(d-a); 
    fix[1] = complex(pHUGE,pZERO);
    return fix; 
  }

  tr = a + d; 
  discr = tr * tr - integer(4); 
  if (small(discr)) {
    fix = rvector(1); 
    fix[0] = (a-d)/(pTWO * c); 
    return fix; 
  }

  discr = sqrt(discr);
  fix[0] = (diff + discr)/(pTWO * c); 
  fix[1] = (diff - discr)/(pTWO * c); 

  /* put the fixed points into canonical order, repelling, attracting.
     derivative is 1/(cz+d)^2. we consider (c fix.end[0] + d) */ 
  pari inv_sqrt_deriv = gabs(c * fix[0] + d); 

  if (small(inv_sqrt_deriv - pONE))
    return fix; /* pure rotation */

  if (inv_sqrt_deriv < pONE) 
    return fix; /* fix[0] repelling */
 
  /* if we got here fix.end[0] was attracting */
  pari tmp = fix[0]; 
  fix[0] = fix[1];
  fix[1] = tmp;
  return fix; 
}   

static pari apply_mt(pari const& m, pari const& z)
{
  pari zz = rvector(2);
  zz[0] = z; 
  zz[1] = pONE; 
  zz = m * zz; 
  if (small(zz[1])) 
    return complex(pHUGE,pZERO);
  return zz[0]/zz[1]; 
}

static pari map01inf(const pari& uvw)
{
  pari u=uvw[0], v=uvw[1], w=uvw[2];
  pari a, b, c, d;

  if (big(u)) {
    a = w;     b = v-w; 
    c = pONE;  d = pZERO; 
  } else if (big(v)) {
    a = w;     b = -u;
    c = pONE;  d = -pONE; 
  } else if (big(w)) {
    a = v-u;   b = u;
    c = pZERO; d = pONE;
  } else {
    a = w*(u-v); b = u*(v-w);
    c = u-v;     d = v-w;
  }

  pari mx = matrix(2,2);
  mx[0][0] = a; mx[1][0] = b; 
  mx[0][1] = c; mx[1][1] = d; 

  pari det = a*d - b*c; 
  if (small(det)) {
    printf("degenerate triple in map01inf\n");
    return mx; 
  }

  return mx/gsqrt(det); 
}

static bool nearly_nice(pari const& p, int& l)
{
  if (abs(p) > real(10000)) {
    l = 2; 
    return true; 
  }
  if (abs(p).double_value() < 1e-4) {
    l = 0; 
    return true; 
  }
  if (abs(p-pONE).double_value() < 1e-4) {
    l = 1; 
    return true; 
  }
  return false; 
}

bool normalize(GroupPresentation* G, MoebiusTransformation& c, bool report, bool try_nice, vector<FGWord> const& wl)
{
  int n_gens = fg_get_num_generators(G);
  FGWord w, u;
  MoebiusTransformation mt;

  if (fg_Moebius_generators(G)[0].acc_matrix.type()==t_INT) {
    printf("normalize requires accurate matrices\n");
  }

  // find 3 fixed points
  pari fp, fixed = rvector(3); 
  int i, j, k, m, nf=0, wn=0, l;
  for (i=0; i<30 && nf<3; ++i) {

    if (wn < wl.size()) { 
      w = wl[wn]; ++wn; 
    } else {
      u.next_word(n_gens);
      w = u; 
      if (inverse(w) < w) continue; // skip inverses
    }
    mt = word_to_Moebius(G, w);
    if (!mt.parity) continue; // skip or. rev. 

    fp = fixed_points(mt.acc_matrix);
    m = length(fp); 
    if (m==3) continue; // ignore identity elements

    // add fixed points to our list
    for (j=0; j<m; ++j) {
      for (k=0; k<nf; ++k) // fp[j] in fixed? 
	if (same_point(fixed[k],fp[j])) break;
      if (k < nf) continue; 

      if (try_nice) {
	if (!nearly_nice(fp[j],l) || 
	    fixed[l].type()!=t_INT) // Already have it. 
	  continue;
	fixed[l] = fp[j];
	++nf; 
      } else {
	fixed[nf] = fp[j];
	++nf; 
      }
      if (nf==3) break; 
    }
  }

  if (nf < 3) return false; 

  // want first two to be 0,inf. 
  if (!try_nice) {
    pari tmp = fixed[2];
    fixed[2] = fixed[1];
    fixed[1] = tmp; 
  }

  if (report) {
    cout << "Fixed points\n";
    cout << fixed << endl; 
    cout << "Normalized matrices\n";
  }

  // conjugate the matrices. 
  c = MoebiusTransformation(One); 
  c.acc_matrix  = map01inf(fixed);
  set_approx(c); 
  MoebiusTransformation ci = inverse(c); 

  for (i=0; i<n_gens; ++i) {
    fg_Moebius_generators(G)[i] = 
      ci * fg_Moebius_generators(G)[i] * c;
    fg_get_generators(G)[i] = O31_matrix(fg_Moebius_generators(G)[i]); 
    if (report)
      cout << fg_Moebius_generators(G)[i] << endl;
  }
  return true; 
}

