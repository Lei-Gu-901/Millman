#include <iostream>
#include <vector>
#include <fstream>
#include <iterator>
#include <assert.h>
#include <random>
#include <chrono>
#include <stdio.h>
#include <cmath>
#include "neuron.h"
// #include "parameter.h"
using namespace std;

//unit of the rates is 1/ms

double taus = 5.0;                           //current decay time 5ms
double we = 95.0*0.1/3.0;                    //pA  has /C
double win = 50.0*0.1/3.0;                   //pA  has /C
double threshold = 20.0;                    //mV  
int rsN = 6;                                // number of relase site
double taur = 100.0;                      //synapse recover time 100ms
double pr = 0.3;                            //rate for release
double iRC = 0.05;

double Ve = we*taus;
double pVin = win*pr*rsN*taus;
double rec = pr*taur;
int numV = 10000;
double Vol[10000];
double dV = threshold/(numV-1);

double Norm(double &D, double &v)
{
    double norm=0.0;
    double con =exp((0.25*threshold*threshold*iRC-0.5*v*threshold)/D);
    double x;
    for(int i=0;i<numV;++i){
        x = Vol[i];
        norm += 1.0/con-exp(x/D*(v-0.5*iRC*x))*con;
    }
    return norm*dV;
}

int main()
{   
    for (int i=0;i<numV;++i){  // generate voltage grid of [0, 20] meV
        Vol[i]=i*dV;
    }
    
    double utillp, rate, norm;
    double Dsum, vsum, factor;
    for (int i=0; i<200; i++){
            double ratep = 0.001+i*0.001;
            utillp = 1.0/(1.0+pr*taur*ratep);
            // factor = 1-(ratep*taus*(1-exp(-1/(ratep*taus))));    //the rescaling factor proposed in Hidalgo et al. 2012
            factor = 1;
            double Vin = pVin*factor;
            Dsum = 7.5*utillp*utillp*Vin*Vin*ratep;
            Dsum *= 0.5;
            vsum = 7.5*utillp*Vin*ratep;
            norm = Norm(Dsum, vsum);
            rate = (vsum-iRC*threshold)*exp(-(0.25*threshold*threshold*iRC-0.5*vsum*threshold)/Dsum);
            rate /= norm;
            //print the presynaptic and postsynaptic firing rate, separated by "..."
            cout<< ratep<<"..."<<rate <<endl;
    }
}