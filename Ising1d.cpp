#include<iostream>
#include<random>
#include<omp.h>
#include<cmath>
#include<vector>
const double J=1;
const double k=1;
const int N=500;
void fill(std::vector<double> &s);
/*void change();
double energy();
double total_energy();
double magnetization();
double fluctuation();
double specific_heat();*/
int main(int argc, char **argv)
{
  double T=10000;
  std::vector<double> s(N);
  fill(s);

  

}
void fill(std::vector<double> &s)
{
  int seed=5;
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> dis(0,1);
  
  for (int i=0; i<N; i++)
    {
      double n=dis(gen);
      if(n<=0.5)
	s[i]=-1;
      else if(n>0.5)
	s[i]=1;
    }
}

void change(std::vector<double> &s)
{
 int seed=2;
  std::mt19937 gen(seed);
  std::uniform_int_distribution<int> dis(0,N-1);
  int i=dis(gen);
  s[i]*=-1;
}
