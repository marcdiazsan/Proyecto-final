#include<iostream>
#include<random>
#include<omp.h>
#include<cmath>
#include<vector>
const double J=1;
const double k=1;
const int N=500;
void fill(std::vector<double> &s);
void change(std::vector<double> &s);
double initial_energy(std::vector<double>s);
double change_energy(std::vector<double>s, int i);
/*double total_energy();
double magnetization();
double fluctuation();
double specific_heat();*/
int main(int argc, char **argv)
{
  double T=10000, energy;
  std::vector<double> s(N);
  int index=0;
  fill(s);
  energia=initial_energy(s)
  change(s,index);
  

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

void change(std::vector<double> &s, int &i)
{
 int seed=2;
  std::mt19937 gen(seed);
  std::uniform_int_distribution<int> dis(0,N-1);
  i=dis(gen);
  s[i]*=-1;
}

double initial_energy(std::vector<double>s )
{
  double sum=0;
  for(int i=0; i<N;i++)
    {
      if(i==N-1)
	{
	  sum+=s[i]*s[0];
	}
      else{

	sum+=s[i]*s[i*1];

    }
    }
  return sum*J;
}

double change_energy(std::vector<double>s, int i)
{


}
