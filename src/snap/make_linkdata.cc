#include <cstdlib>
#include <cstdio>

using std::sscanf;
using std::fopen;
using std::fclose;
using std::printf;
using std::fgets;

#define BSIZE 256

struct link {
  double volume;
  char dowker[26];
};

int main()
{
  char buf[BSIZE];
  char codebuf[BSIZE];

  size_t linksize = sizeof(link);
  printf("sizeof link is %d\n", linksize); 

  FILE* out = fopen("links_by_volume", "wb");

  if (!out) {
    printf("problem opening file links_by_volume for append\n");
    return 0;
  }

  link L;

  int line=1, i, nargs; 

  while (fgets(buf, BSIZE, stdin)) {

    nargs = sscanf(buf, "%s %lf ", codebuf, &L.volume);
    if (nargs != 2) {
      printf("invalid data on input line %d\n", line); 
      fclose(out);
      return 0;
    }

    for (i=0; buf[i] && i<26; ++i)
      L.dowker[i] = codebuf[i];
    for (; i<26; ++i) 
      L.dowker[i] = '\0';

    if (fwrite((void*)&L, linksize, 1, out) != 1) {
      printf("trouble writing link number %d\n", line);
      fclose(out); 
      return 0; 
    }

    ++line;
  }

  printf("wrote %d links\n", line-1);
  return 0;
}





