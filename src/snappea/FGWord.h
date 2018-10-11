#ifndef _FGWord_
#define _FGWord_

#include <string>
#include <iostream>

/* An FGWord is an element of the free group on the symbols 1,2,3,...
   with the inverses being represented as -1,-2,-3,... .
   Multiplication keeps words freely reduced.  */

char invchar(char c);

int toint(char a);

char tochar(int i);

class FGWord 
{ 
  std::string w; 
public: 
  FGWord() {} 
  FGWord(int g); 
  FGWord(const char* c); // Read in form abAB or "(1 2 -1 -2)" 

  operator std::string() const; // Convert to form abAB 

  void invert(); 

  FGWord& operator *= (int c);
  FGWord& operator *= (FGWord b);

  // char& operator [] (int i) 
  //  { return w[i]; }
  int operator [] (int i) const;
  int length() const 
    { return w.size(); }

  void left_multiply_by(int c);
  void left_multiply_by(FGWord a);

  void next_word(int n_gens); 

  friend FGWord operator * (FGWord a, FGWord b);
  friend bool operator == (FGWord const& a, FGWord const& b); 
  friend bool operator < (FGWord const& a, FGWord const& b) // short-lex. 
    { return a.length() < b.length() || (a.length()==b.length() && a.w < b.w); }

  std::string numeric_form() const; // Convert to form "(1 2 -1 -2)".
  void print() const;
};

FGWord inverse(FGWord a);
std::ostream& operator << (std::ostream& out, const FGWord& a);

#endif
