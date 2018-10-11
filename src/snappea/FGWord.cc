#include "FGWord.h"

#include <algorithm>
#include <ctype.h>
#include <stdio.h>

using namespace std; 

char invchar(char c)
{ return (c >= 'a') ? (c + 'A' - 'a') : (c + 'a' - 'A'); }

int toint(char a)
{ return (a >= 'a') ? a - 'a' + 1 : -(int)a + 'A' - 1; }

char tochar(int i)
{ return (i > 0) ? char(i + 'a' - 1) : char(-i+ 'A' - 1); }

int FGWord::operator [] (int i) const
{ return toint(w[i]); }

FGWord::FGWord(int g)
{
  if (g) w += tochar(g); 
}

FGWord::FGWord(const char* c)
{
  int n; 
  while (*c && isspace(*c)) c++; // eat space
  if (*c != '(') { // read letter form
    w = c; 
  } else { // read numeric form
    c++; // eat open paren
    while (*c && isspace(*c)) c++; // eat space
    if (*c == '\0') return; // premature end of line
    while (*c != ')') {
      if (!sscanf(c, "%d", &n)) return;
      w += tochar(n); 
      while (*c && *c != ')' && !isspace(*c)) c++; // eat to end of word or closing paren. 
      while (*c && isspace(*c)) c++; // eat space
      if (*c == '\0') return; // premature end of line
    }
  }
}

FGWord::operator string () const
{
  return w; 
}

string FGWord::numeric_form() const
{
  char buf[20];
  string s = "("; 
  int i, n = length(); 
  for (i=0; i<n; i++) { 
    if (i!=0) s += ' '; 
    sprintf(buf, "%d", toint(w[i])); 
    s += buf;
  }
  s += ")";
  return s; 
}

void FGWord::invert()
{
  int i, n = w.length(); 
  char c; 

  // Reverse and negate simultaneously. 
  for (i=0; i<(n+1)/2; i++) {
    c = w[i];
    w[i] = invchar(w[n-1-i]);
    w[n-1-i] = invchar(c); 
  }
}

FGWord inverse(FGWord a)
{
  a.invert(); 
  return a; 
}

FGWord& FGWord::operator *= (FGWord b)
{
  int i, j, n = min(length(), b.length());

  // Find first point where words are not mutually inverse. 
  for (i=0, j=length()-1; i<n; i++,j--)
    if (invchar(w[j]) != b.w[i]) break; 

  // w.erase(j+1);
  string empty; 
  w.replace(j+1, length()-j-1, empty);
  w += b.w.substr(i);
  return *this; 
}

FGWord& FGWord::operator *= (int c)
{
  char bi = tochar(-c); 
  if (length() && w[length()-1] == bi) {
    string empty; 
    // w.erase(length()-1);
    w.replace(length()-1, 1, empty);
  } else {
    w += tochar(c);
  }
  return *this;
}

FGWord operator * (FGWord a, FGWord b)
{
  int i, j, n = min(a.length(), b.length());

  // Find first point where words are not mutually inverse. 
  for (i=0, j=a.length()-1; i<n; i++,j--)
    if (invchar(a.w[j]) != b.w[i]) break; 

  FGWord c; 
  c.w = a.w.substr(0,j+1) + b.w.substr(i); 
  return c; 
}

void FGWord::left_multiply_by(FGWord a)
{
  int i, j, n = min(a.length(), length());

  // Find first point where words are not mutually inverse. 
  for (i=0, j=a.length()-1; i<n; i++,j--)
    if (invchar(a.w[j]) != w[i]) break; 

  w.replace(0, i, a.w.substr(0,j+1));
}
  
void FGWord::left_multiply_by(int c)
{
  if (length() && w[0] == tochar(-c)) {
    // w.erase(0,1);
    string empty;
    w.replace(0,1,empty);
  } else {
    w.insert(w.begin(),tochar(c));
  }
}

void FGWord::next_word(int n_gens)
{
  int l = w.length(); 
  if (l==0) {
    w = string(1,'a');
    return; 
  }
  if (toint(w[l-1])==-n_gens) {
    w = w.substr(0,l-1); /* get rid of last letter */
    next_word(n_gens); 
    w += 'a';
    l = w.length(); 
  } else {
    w[l-1] = (w[l-1] < 'a') ? invchar(w[l-1])+1 : invchar(w[l-1]); 
  }
  if (l > 1 && (w[l-2] == invchar(w[l-1])))
    next_word(n_gens); 
  return;
}

ostream& operator << (ostream& out, const FGWord& a)
{
  return out << string(a);
}

bool operator == (FGWord const& a, FGWord const& b)
{
  return a.w==b.w; 
}

void FGWord::print() const
{
  cout << *this << endl;
}

