#include <iostream>
#include <stdlib.h>
#include <vector>
#include <math.h>

//Se definen Constantes globales
const int N = 2000;
const int rep = 50;//repeticiones de cada temperatura
//const double kT = 1.5;
const int PASOS = N*100;

//Funciones
double Energia(std::vector<double> & estado);
void print(const std::vector<double> & estado);
void Metropolis(std::vector<double> & estado, int i, double kT);
double vecinos(std::vector<double> & estado, int i);
double deltaE(std::vector<double> & estado, int i);
double Magnetico(std::vector<double> & estado);
void inicio(std::vector<double> & estado,double  M, double E);

int main(void){

  
  //Inicializa vector y variables e imprime
  double M=0,E=0, Mprom=0, Eprom=0;
  std::vector<double> estado(N);
  srand (time(NULL));
  //print(estado);
  
  //Corridas de metrópolis
  for (double kT =0.1; kT<5.0;kT+=0.05){
    Mprom=0;
    Eprom=0;
    for(int ii = 0; ii<rep; ii++){
      inicio(estado, M, E);
      for (int n=0;n < PASOS;n++){
	int i = rand() % N;//Elige posición al azar a cambiar spin
	Metropolis(estado, i, kT);
      }
      Mprom+=Magnetico(estado);
      Eprom+=Energia(estado);
      // std::cout << Mprom << " " <<Eprom<<std::endl;
    }
    std::cout << kT << " \t" << Mprom/(N*rep) << "\t "<< Eprom/(N*rep)<<std::endl;
  }
  return 0;
  
}

void inicio(std::vector<double> & estado, double M, double E){

  for (int i = 0; i < N; i++){
    estado[i]=1.0;
  }
}


double Magnetico(std::vector<double> & estado){
  double M=0;
  for(int i=0;i<N;i++){
    M+=estado[i];
  }
  return abs(M);
}


void print(const std::vector<double> & estado){
  for(int i = 0;i < N; i++){
    std::cout << estado[i] << " ";
  }
  std::cout << "\n";
}

void Metropolis(std::vector<double> & estado, int i, double kT){
    double p = (rand()/(double)(RAND_MAX));//Número al azar entre 0 y 1
    
    //Calcula el cambio en energía si se cambia el spin
    double cambioE = deltaE(estado, i);

    if (cambioE < 0){
      estado[i]*=-1;
    }
    else if (exp(-cambioE/kT) >= p ){
      estado[i]*=-1;
    }
    else{
      estado[i] *=1;
    }
}

double vecinos(std::vector<double> & estado, int i){//Da suma de los valores de spines de vecinos
  return estado[(i-1+N) % N]+estado[(i+1) % N];
}

double deltaE(std::vector<double> & estado, int i){//Calcula el cambio de energía al cambiar el spin de i (J=1)
  return 2.0 *estado[i] * vecinos(estado, i);
}

double Energia(std::vector<double> & estado){
  double E=0;
  for(int i=0;i<N;i++){
    E+=-estado[i]*estado[(i+1)% N];
  }
  return E;
}


