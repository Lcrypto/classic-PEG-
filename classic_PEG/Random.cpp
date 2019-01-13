#include "Random.h"

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <algorithm>


Random::Random()
  : m_seed{ 987'654'321 }, m_seed_u{ 123'456'789 }
{

}

Random::~Random() = default;

void Random::Sort(int a[], int size)
{
  return std::sort( a, a + size );
}

double Random::gauss(double sdev, double mean)
{ 
  double sum=0.0;
  for (int i=1;i<=12;i++)
    { 
      m_seed_u = 1664525lu * m_seed_u + 123456789lu; 
      sum=sum+double(m_seed_u); 
    }
  return (double(sum)/4.29497e9-6.0)*sdev+mean;
}

double Random::uniform(double a, double b)
{
  double t;
  for(int i=0; i<10;i++){
    m_seed=2045*m_seed+1;
    //seed=seed -(seed/1048576)*1048576;
    m_seed%=1048576;
  }
  t=m_seed/1048576.0;
  t=a+(b-a)*t;
  return(t);
}

int Random::uniform(int  a, int b) // [a, b-1]
{
  double t;
  int i, tt;
  if(b==a+1) return(a);
  for(i=0; i<10;i++){
    m_seed=2045*m_seed+1;
    //seed=seed -(seed/1048576)*1048576;
    m_seed%=1048576;
  }
  t=m_seed/1048576.0;
  t=a+(b-a)*t;
  tt=(int)t;
  if(tt<a) tt=a;
  else if(tt>=b) tt=b-1;
  return(tt);
}

int Random::nonUniform(int  a, int b) // [a, b-1]
{
  double t;
  int i, tt;
  if(b==a+1) return(a);
  for(i=0; i<10;i++){
    m_seed=2045*m_seed+1;
    //seed=seed -(seed/1048576)*1048576;
    m_seed%=1048576;
  }
  t=m_seed/1048576.0;
  t=a+(b-a)*pow(t, 0.6667); //t^1.5 
  tt=(int)t;
  if(tt<a) tt=a;
  else if(tt>=b) tt=b-1;
  return(tt);
}

/*
void m_uniform(int a , int b, long int *seed, int *itmp, int num)
{
  int index, imed, i, k;

  index=0;
  itmp[0]=uniform(a, b, seed);
 Loop1:
  imed=uniform(a, b, seed);
  k=0;
  for(i=0;i<=index;i++)
    if(imed==itmp[i]) {k=1; break;}
  if(k==0) {index++; itmp[index]=imed;}
  if(index<num-1) goto Loop1;
  bubble(itmp, num); //bubble sorting
}
*/









