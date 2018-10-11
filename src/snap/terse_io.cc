#include "terse_io.hh"
#include "snappea/terse_triangulation.h"
#include <set>
#include <algorithm>
#include "file_array.hh"
#include "printable.hh"

using std::string;
using std::multiset;
using std::pair;
using std::ostream;
using std::cout;
using std::endl;

class bitstream {
  size_t n;
  unsigned char* B;
  int pos;

public:
  bitstream(size_t n_bytes) : n(n_bytes), pos(0)
  { B=new unsigned char[n]; clear(); }
  ~bitstream() { delete[] B; }

  void put(int bits, unsigned int value);
  unsigned int get(int bits);

  void set_pos(int p) { pos = p; }
  void clear() { int i; for (i=0; i<n; ++i) B[i] = '\0'; }; 

  bool read(FILE* fp); 
  bool read(FILE* fp, int k); 
  bool write(FILE* fp) const; 

  void print() const; 
};

static int file_size(FILE* fp)
{
  if (fseek(fp, 0, SEEK_END)) return -1; 
  long int p = ftell(fp); 
  if (p<0) return -1; 
  return p; 
}

bool bitstream::read(FILE* fp)
{
  pos = 0; 
  return fread((void*)B, sizeof(unsigned char), n, fp)==n; 
}

// get k'th bitstream of this size from a file 
// assumed to contain at least k+1 bitstreams. 

bool bitstream::read(FILE* fp, int k)
{
  if (fseek(fp, n*k, SEEK_SET)) return false; 
  return read(fp); 
}

bool bitstream::write(FILE* fp) const
{
  // if (fseek(fp, 0, SEEK_END)) return false; 
  return fwrite((void*)B, sizeof(unsigned char), n, fp)==n;
}

void bitstream::put(int bits, unsigned int value)
{
  int i, bit_offset, bits_this_byte;
  unsigned char mask; 

  while (bits > 0) {
    // get position in B. 
    i = pos/8;
    bit_offset = pos%8;

    if (i>=n) { printf("bitstream overflow!\n"); return; }

    // work out how many bits we'll write this time. 
    bits_this_byte = 8-bit_offset;
    if (bits < bits_this_byte) 
      bits_this_byte = bits; 

    // write the data. 
    mask = 0xff >> (8-bits_this_byte); 
    B[i] |= (value & mask) << bit_offset;  
    value >>= bits_this_byte; 
    
    // update everything. 
    pos    += bits_this_byte; 
    bits   -= bits_this_byte; 
  }
}

unsigned int bitstream::get(int bits)
{
  unsigned int value = 0; 

  int i, bit_offset, bits_this_byte, data_offset = 0;
  unsigned char mask; 
  unsigned int data; 

  while (bits > 0) {
    // get position in B. 
    i = pos/8;
    bit_offset = pos%8;

    if (i>=n) { printf("bitstream overread!\n"); return value; }

    // work out how many bits we'll write this time. 
    bits_this_byte = 8-bit_offset;
    if (bits < bits_this_byte) 
      bits_this_byte = bits; 

    // read the data. 
    mask = (0xff >> (8-bits_this_byte)) << bit_offset; 
    data = (B[i] & mask) >> bit_offset; 
    value |= data << data_offset; 
    
    // update everything. 
    pos         += bits_this_byte; 
    bits        -= bits_this_byte; 
    data_offset += bits_this_byte; 
  }

  return value; 
}




void bitstream::print() const
{
  int i,j;
  unsigned char u; 
  for (i=n-1; i>=0; --i) {
    u = B[i]; 
    for (j=7; j>=0; --j) {
      printf("%d", (u>>j)&1);
    }
  }
}

static int num_bits(unsigned long n)
{
  int nbits = 0;
  while (n) { n >>= 1; ++nbits; }
  return nbits; 
}


int terse_triangulation_bytes(int n_tetrahedra)
{
  if (n_tetrahedra < 1) return 0; // avoid possible trouble.  

  int tet_num_bits = num_bits(n_tetrahedra - 1); 
  int num_bits = 2*n_tetrahedra + (5 + tet_num_bits) * (n_tetrahedra + 1); 

  int num_bytes = (num_bits + 7)/8; // divide by 8, rounding up. 

  return num_bytes; 
}

int num_terse(string const& dir, int n_tet)
{
  char fname[25];
  sprintf(fname, "/trsT%d", n_tet);
  string name = dir + fname; 

  FILE* fp = fopen(name.c_str(),"rb");
  if (!fp) return 0; 

  int num = file_size(fp)/terse_triangulation_bytes(n_tet);
  fclose(fp);

  return num;
}
  

bool write_terse(FILE* fp, TerseTriangulation const& T, int& index)
{
  int n_bytes = terse_triangulation_bytes(T.num_tetrahedra);
  bitstream B(n_bytes); 

  index = file_size(fp)/n_bytes;

  int i; 
  for (i=0; i < 2*T.num_tetrahedra; i++)
    B.put(1, T.glues_to_old_tet[i] ? 1:0);

  int tet_num_bits = num_bits(T.num_tetrahedra - 1); 
  for (i=0; i < T.num_tetrahedra+1; i++) {
    B.put(tet_num_bits, T.which_old_tet[i]);
    B.put(5, index_by_permutation[T.which_gluing[i]]);
  }

  return B.write(fp);
}

void print_terse(TerseTriangulation* tt)
{
  int i, n=tt->num_tetrahedra;

  printf("Num tetrahedra = %d\n", n); 
  for (i=0; i<2*n; ++i) {
    printf("%d", (int)tt->glues_to_old_tet[i]);
  }
  printf("\n");
  for (i=0; i < n+1; i++) {
    printf("%d:%d ", tt->which_old_tet[i], 
	   index_by_permutation[tt->which_gluing[i]]); 
  }
  printf("\n");
}

bool read_terse(FILE* fp, TerseTriangulation& T, int index)
{
  bitstream B(terse_triangulation_bytes(T.num_tetrahedra)); 

  if (!B.read(fp,index)) return false; 

  int i; 
  for (i=0; i < 2*T.num_tetrahedra; i++)
    T.glues_to_old_tet[i] = B.get(1);

  int tet_num_bits = num_bits(T.num_tetrahedra - 1); 
  for (i=0; i < T.num_tetrahedra+1; i++) {
    T.which_old_tet[i] = B.get(tet_num_bits);
    T.which_gluing[i] = permutation_by_index[B.get(5)];
  }
  return true;
}

Triangulation* read_terse(string const& dir, int n_tet, int index)
{
  char fname[25];
  sprintf(fname, "/trsT%d", n_tet);
  string name = dir + fname; 

  FILE* fp = fopen(name.c_str(),"rb");
  if (!fp) return 0; 

  TerseTriangulation* tt = alloc_terse(n_tet);
  tt->num_tetrahedra = n_tet;

  bool res = read_terse(fp, *tt, index);
  fclose(fp);
  if (!res) { 
    free_terse_triangulation(tt); 
    return 0; 
  }

  Triangulation* T = NEW_STRUCT(Triangulation);
  initialize_triangulation(T);
  terse_to_tri(tt, T);

  free_terse_triangulation(tt);

  char manifold_name[25];
  sprintf(manifold_name, "T%d.%d", n_tet, index);
  set_triangulation_name(T, manifold_name);

  return T;
}

bool save_terse(string const& dir, Triangulation* T, int& index)
{
  char fname[25];
  sprintf(fname, "/trsT%d", get_num_tetrahedra(T));
  string name = dir + fname; 

  FILE* fp = fopen(name.c_str(),"a+b");
  if (!fp) return false; 

  TerseTriangulation* tt = tri_to_terse(T);

  bool res = write_terse(fp, *tt, index);
  fclose(fp);
  free_terse_triangulation(tt);

  return res; 
}

struct terse_key {
  double vol;
  unsigned char n_tet;
  unsigned int index;

  terse_key(double v=0.) : vol(v), n_tet(0), index(0) {}

  friend bool operator < (terse_key const& a, terse_key const& b)
  { return a.vol < b.vol - 1e-9; }
  friend ostream& operator << (ostream& out, terse_key const& a)
  { return out << a.vol << ' ' << int(a.n_tet) << '.' << a.index; }
};


class manifold_db {
  string dir;
  multiset<terse_key> pending_keys;
  FILE* pending_log; 

  typedef file_array<terse_key>::const_iterator key_iter;
  bool isometric(terse_key const& K, Triangulation* T) const;

public:
  manifold_db(string const& dir);
  ~manifold_db() { if (pending_log) fclose(pending_log); }

  bool locate(Triangulation* T, int& index, int& nt) const;
  bool insert(Triangulation* T, int& index, int& nt);

  void save_pending_keys();
  void print() const; 
};

manifold_database::manifold_database(string const& dir)
{ rep = new manifold_db(dir); }
manifold_database::~manifold_database()
{ delete rep; }
bool manifold_database::locate(Triangulation* T, int& index, int& nt) const
{ return rep->locate(T,index, nt); }
bool manifold_database::insert(Triangulation* T, int& index, int& nt)
{ return rep->insert(T,index, nt); }
void manifold_database::save_pending_keys()
{ rep->save_pending_keys(); }
void manifold_database::print() const
{ rep->print(); }

void manifold_db::print() const
{
  file_array<terse_key> keys((dir+"/terse_by_volume").c_str());

  cout << "keys \n" << PSeq(keys, "\n") << endl; 

  if (pending_keys.size()) 
    cout << "pending keys\n" << PSeq(pending_keys,"\n") << endl;
}

manifold_db::manifold_db(string const& d)
 : dir(d)
{
  pending_log = fopen((dir + "/pending_log").c_str(), "r+b");

  if (!pending_log) {
    return;
  }

  terse_key K;
  while (fread((void*)&K, sizeof(terse_key), 1, pending_log)==1)
    pending_keys.insert(K);

  fseek(pending_log, 0, SEEK_END);
}

bool manifold_db::isometric(terse_key const& K, Triangulation* T) const
{
  Triangulation* U = read_terse(dir, K.n_tet, K.index);
  if (!U) {
    printf("problem finding triangulation for key in database\n");
    return false;
  }
  Boolean isometric;
  if (compute_isometries(T,U,&isometric,0,0)!=func_OK) {
    printf("problem checking an isometry\n");
    return false;
  }
  return isometric;
}


bool manifold_db::locate(Triangulation* T, int& index, int& nt) const
{
  file_array<terse_key> keys((dir+"/terse_by_volume").c_str());
  int first;
  double vol = volume(T,0);
  int n = keys.find(vol,first);

  int i;
  for (i=first; i<first+n; i++)
    if (isometric(keys[i],T)) break;
  if (i<first+n) {
    index = keys[i].index;
    nt = keys[i].n_tet; 
    return true;
  }

  typedef multiset<terse_key>::const_iterator pkey_iter;


  pair<pkey_iter,pkey_iter> range = 
    equal_range(pending_keys.begin(), pending_keys.end(), terse_key(vol));

#if 0
  printf("checking %d pending keys\n", 
	 std::distance(range.first, range.second)); 
#endif

  pkey_iter it;
  for (it=range.first; it!=range.second; ++it)
    if (isometric(*it,T)) break;
  if (it!=range.second) {
    index = it->index;
    nt = it->n_tet; 
    return true;
  }

  return false;
}

bool manifold_db::insert(Triangulation* T, int& index, int& nt)
{
  if (locate(T,index,nt)) return false;

  if (!save_terse(dir, T, index)) {
    index = -1;
    printf("problem saving a terse triangulation\n");
    return false;
  }

  double vol = volume(T,0);
  terse_key K(vol);
  K.n_tet = get_num_tetrahedra(T);
  K.index = index;
  nt = K.n_tet; 

  pending_keys.insert(K);

  string logname = dir+"/pending_log";
  pending_log = fopen(logname.c_str(),"ab");
  if (!pending_log) {
    printf("problem opening pending_log for append\n");
    return true; 
  }
  if (fwrite((void*)&K, sizeof(terse_key), 1, pending_log)!=1) {
    printf("problem updating log of keys\n");
  }
  return true;
}

class terse_key_inserter : 
  public std::iterator<std::output_iterator_tag, int>
{
  FILE* fp;
 public:
  terse_key_inserter(FILE* f) : fp(f) {}

  terse_key_inserter& operator * ()    { return *this; }
  terse_key_inserter& operator ++()    { return *this; }
  terse_key_inserter& operator ++(int) { return *this; }
  terse_key_inserter& operator = (terse_key const& K) 
  { fwrite((void*)&K, sizeof(terse_key), 1, fp); return *this; }
};

void manifold_db::save_pending_keys()
{

#if 0
  {
    file_array<terse_key> keys((dir+"/terse_by_volume").c_str());
    printf("have %d key(s)\n", keys.size()); 
  }

  printf("have %d pending key(s)\n", pending_keys.size()); 
#endif

  if (!pending_keys.size()) return;

  string kf_name = dir + "/terse_by_volume";
  string kf_old  = kf_name + "~";
  rename(kf_name.c_str(), kf_old.c_str());

  FILE* k_out = fopen(kf_name.c_str(), "wb");
  if (!k_out) {
    printf("problem opening new key file for writing\n");
    rename(kf_old.c_str(), kf_name.c_str());
    return;
  }

  terse_key_inserter tki(k_out);
  file_array<terse_key> okeys(kf_old.c_str());

  if (!okeys.size()) {
    copy(pending_keys.begin(), pending_keys.end(), tki);
  } else {
    merge(okeys.begin(), okeys.end(), 
	  pending_keys.begin(), pending_keys.end(), tki);
  }

  fclose(k_out);
  pending_keys.clear();

  fclose(pending_log);
  string logname = dir+"/pending_log";
  rename(logname.c_str(), (logname+"~").c_str());
  pending_log = fopen(logname.c_str(),"ab");
}


#if 0
int main()
{
  int i; 
  
#if 0
  printf("n-tet\tn-bytes\n");
  for (i=1; i<25; i++) {
    printf("%d\t%d\n", i, terse_triangulation_bytes(i));
  }
#endif
  bitstream B(4);

  for (i=0; i<6; i++) B.put((i%2)?5:4,i);
  B.print();
  printf("\n"); 

  B.set_pos(0); 

  unsigned int j; 
  for (i=0; i<6; i++) { j = B.get((i%2)?5:4); printf("%u ",j); }
  printf("\n"); 

  // use mode: "rb" "wb" "ab" "r+b" "w+b" "a+b"
  bool open(const char* name, const char* mode)
  { fp = fopen(name, mode); return fp!=0; }
  void close() { if (fp) fclose(fp); }

  FILE* fp = fopen("my_bitstream", "ab");
  B.write(fp); 
  fclose(fp); 

  B.clear(); 

  fp = fopen("my_bitstream", "rb");

  int s = B.file_size(fp); 
  printf("my_bitstream contains %d bitstreams\n", s); 


  if (!B.read(fp, 0)) {
    printf("error reading a bitstream\n"); 
  } else {
    printf("bitstream 0 in file\n"); 
    B.print(); 
    printf("\n");
  }

  fclose(fp); 

  for (i=0; i<6; i++) { j = B.get((i%2)?5:4); printf("%u ",j); }
  printf("\n"); 

  return 0;
}
#endif
