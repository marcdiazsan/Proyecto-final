#include<iostream>
#include<random>
#include<omp.h>
#include<cmath>
#include<vector>


const double J=1;
const double k=1;
const int N=10000;
const double T=10000;
const double epsilon=0.1;


void fill(std::vector<double> &s);
void change(std::vector<double> &s,int &i);
double initial_energy(std::vector<double>s);
double change_energy(std::vector<double>s, int i, double energy);
void energy_comparision(double &energy, double &tmp_energy,int &k);
double magnetization(std::vector<double>s);
/*double fluctuation();
double specific_heat();*/


int main(int argc, char **argv)
{
  double  energy, energy_tmp, mgt,e2=0,u2,e=0,u,c;
  std::vector<double> s(N);
  int index=0, n=0;//n, contador de pasos
  fill(s);
  
  /* for(int i=0; i<N; i++)
    {
      std::cout<<s[i];
    }
  
    std::cout<<endl;*/

  energy=initial_energy(s);
  std::cout<<energy<<std::endl;

    /*for(int i=0; i<N; i++)
    {
      std::cout<<s[i];
    }
    std::cout<<std::endl;*/
    energy=initial_energy(s);
    std::cout<<energy; 
    change(s,index);
    energy_tmp=change_energy(s,index,energy);
    std::cout<<energy_tmp;
    energy_comparision(energy,energy_tmp,n);  //al llegar al equilibrio, la energia total es la energy que quede al final(no se necesita una función mas)
    e+=energy;
    u=e/n; // Se calcula el promedio de la energia
    e2+=(energy*energy);
    u2=e2/n; // Se calcula el promedio de la energia al cuadrado
    std::cout<<" "<<u2<<std::endl;
  
  
  mgt=magnetization(s);
  c=(u2-u*u)/(n*k*T*T);
}
void fill(std::vector<double> &s)
{
  std::random_device rd;
  std::mt19937 gen(rd());
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
  std::random_device rd;
  std::mt19937 gen(rd());
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
	  sum+=(s[i]*s[0]);
	}
      else{

	sum+=(s[i]*s[i+1]);

    }
    }
  return -sum*J;
}

double change_energy(std::vector<double> s, int i,double energy)
{
  energy-=-J*(-s[i]*s[i-1]-s[i]*s[i+1]);//se restan las 2 de la energía previa
  energy+=-J*(s[i]*s[i-1]+s[i]*s[i+1]);// se suman las 2 de la energía al cambiar un unico estado
  return energy;

}
void energy_comparision(double &energy, double &tmp_energy, int &n)
{ std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dis(0,1);
  if(tmp_energy<=energy)
    {
      energy=tmp_energy;  //accept the new configuration
    }
  else
    {
      double deltaE=tmp_energy-energy;
      double P=std::exp((-deltaE)/(k*T));
      double rj=dis(gen);
      if(P>=rj)
	{
	  energy=tmp_energy;
	}
      else
	{
	  energy=energy; //energy keeps the same
	}
      
    }
  n+=1;
  
}
double magnetization(std::vector<double>s)
{
  double m;
  for(int k=0; k<N;k++)
    {
      m+=s[k];
    }
  return m;
}

