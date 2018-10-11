#include "SL2RMatrix.h"

void SL2R_to_o21(SL2RMatrix source,O21Matrix target)
{
  SL2RMatrix m[3] = {{1,0,0,1},{0,1,1,0},{1,0,0,-1}};
  SL2RMatrix tmp[3];
  int i;
  
 
  SL2RMatrix_conjugate(m[0],source,tmp[0]);
  SL2RMatrix_conjugate(m[1],source,tmp[1]);
  SL2RMatrix_conjugate(m[2],source,tmp[2]);

  for(i=0;i<3;i++)
    {
      target[0][i] = (tmp[i][0][0]+tmp[i][1][1])/2;
      target[1][i] = tmp[i][0][1];
      target[2][i] = (tmp[i][0][0]-tmp[i][1][1])/2;

      
    }
}

void o21_to_o31(O21Matrix source, O31_matrix& target)
{
  int i,j;
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
  target(i,j) = source[i][j];

  for(i=0;i<3;i++)
    target(i,3) = 0;
   for(j=0;j<3;j++)
     target(3,j) = 0;
  target(3,3) = 1;

}

void SL2RMatrix_product(SL2RMatrix a,SL2RMatrix b, SL2RMatrix product)
{
  register int i,j,k;
  register double sum;
  SL2RMatrix tmp;

  for(i=0;i<2;i++)
    for(j=0;j<2;j++)
      {
    sum = 0.0;
    for(k=0; k<2; k++)
      sum += a[i][k]*b[k][j];
    tmp[i][j] = sum;
      }
  for(i = 0; i<2; i++)
    for(j=0; j<2 ; j++)
      product[i][j] = tmp[i][j];
}

void SL2RMatrix_transpose(SL2RMatrix t, SL2RMatrix t_transpose)
{
  register int i,j;
  SL2RMatrix tmp;

  for(i=0;i<2;i++)
    for(j=0;j<2;j++)
      tmp[i][j]= t[j][i];

  for(i=0;i<2;i++)
    for(j=0;j<2;j++)
      t_transpose[i][j] = tmp[i][j];
}

void SL2RMatrix_conjugate(SL2RMatrix m,SL2RMatrix t, SL2RMatrix tmT)
{
  
  SL2RMatrix t_transpose,temp;

  SL2RMatrix_transpose(t,t_transpose);
  SL2RMatrix_product(t,m,temp);
  SL2RMatrix_product(temp,t_transpose,tmT);
}















