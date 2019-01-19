#ifndef RANDOM
#define RANDOM

class Random {
public:
  Random();
  ~Random();

  static void Sort(int a[], int size);
  
  double gauss(double sdev, double mean);
  double uniform(double a, double b);
  int uniform(int  a, int b); // [a, b)
  int nonUniform(int  a, int b);

private:
  unsigned long int m_seed;
  unsigned long int m_seed_u;

}; 

#endif
