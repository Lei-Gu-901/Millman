#include <iostream>
#include <vector>
#include <fstream>
#include <iterator>
#include <assert.h>
#include <algorithm>
#include <random>
#include <chrono>
#include <stdio.h>
#include <cmath>
#include "neuron.h"
#include "parameter.h"
using namespace std;
// the simulation is based on the event driven method

unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);
std::uniform_real_distribution<double> uniform(0.0,1.0);
std::exponential_distribution<double> expdis(edt);

Neuron *neurons[num_neuron];                         //array of pointers to neurons
double spiking_time[num_neuron] = {0};                  //time of spikings   
double external_time[num_neuron] = {0};                 //time of external stimulations
double event_time[num_neuron] = {0};                    //time of the latest event (external stimulation or voltage resetting)
double voltage[num_neuron];                           //voltages at the event 
double current[num_neuron];                              //currents at the event
// double spike_curr[num_neuron];                        //currents at the spiking

double next_spiking[num_neuron] = {0};                  //time of the next spiking
double next_exteral[num_neuron] = {0};                  //time of the next extermal spiking
double next_setting[num_neuron] = {0};

int num_firing[num_neuron] = {0};                    //number of spikings for all neurons 
int num_branching[num_neuron] = {0};                 //number of branching from all neurons
double fact = 20.0/3.0;                              // a numberical factor
double inner_time[num_neuron];
const int s_neuron = t_neuron-81;                   //number of neurons counted for spikes
int avalanche_num[num_neuron-s_neuron] = {0};        //the avalance num of the last spiking        
int avalanche_t[t_neuron] = {0};                                        

// function for calculating the spiking time based on voltage and current
double next_time(double &cur, double &vol)
{
    double a = fact*cur+1.0e-15; //
    double d = -(vol+a);
    double delta = a*a*(256.0*a*pow(threshold,3)-27.0*pow(d,4));
    if (delta>0.0){
        return 0.0;
    }else{
        double cq = cbrt(0.5*(27.0*a*d*d-sqrt(-27.0*delta)));
        double cs = 0.5*sqrt(1.0/(3.0*a)*(cq+12.0*a*threshold/cq));
        double ct = 0.5*sqrt(-4.0*cs*cs-d/(a*cs));
        if(cs+ct<1.0){
            return -20.0*log(cs+ct);
        }else{
            return 0.0;
        }
    }
}

void intialization()
{
    // read the files
    int csize = 9117;               //copy the output number in the network generation step
    int esize = 81;                                  //set the first layer as receipt layer
    int connect[csize];
    ifstream con("connection.dat", ios::binary);
    con.read((char*)&connect, csize*sizeof(int));
    int extneuron[esize];
    ifstream ext("extneuron.dat", ios::binary);
    ext.read((char*)&extneuron, esize*sizeof(int));
    ifstream inf("con_info.dat", ios::binary);
    int con_inf[num_neuron];
    inf.read((char*)&con_inf, num_neuron*sizeof(int));
    vector<int> in_con[num_neuron];                  //array of vectors saving the connections outward
    int nind = 0;                                    //indices of current neuron
    int cnum = con_inf[0];
    for (int i=0;i<csize;i++){
        if (i==cnum){
            ++nind;
            cnum+=con_inf[nind];
        }
        in_con[connect[i]].push_back(nind);
    }

    cnum = 0;

    //initialization of neuron status
    for (int i=0; i<num_neuron; i++){
        int start = cnum;
        int osize = con_inf[i];
        cnum += osize;
        int out_con[osize];
        for (int j=start; j<cnum; ++j) {
            out_con[j-start] = connect[j];
        }
        current[i] = 3.0*win*uniform(generator)+0.001*we*uniform(generator);
        voltage[i] = 3.0*uniform(generator);

        double util[osize*rsN];
        for (int j=0; j<osize*rsN;j++){
            util[j] = 1.0;
        }
        int isize = in_con[i].size();
        int ind[isize];
        copy(in_con[i].begin(), in_con[i].end(), ind);

        neurons[i] = new Neuron(osize, isize, out_con, ind, util);
        next_exteral[i] = stimuT;
        spiking_time[i] = -1.0;                         //set the initial spiking time smaller than 0, 
                                                        //so that the first spiking can be ascribed to the external stimulus
        next_setting[i] = stimuT+1.0;
        double ccur = current[i];
        double cvol = voltage[i];
        double ntime = next_time(ccur, cvol);
        if (ntime != 0.0){
            next_spiking[i] = event_time[i]+ntime;
        }else{
                next_spiking[i] = stimuT; 
        }
    }

    for (int i=0; i<esize; ++i){
        int ind = extneuron[i];
        next_exteral[ind] = expdis(generator);
    }
}

int main()
{
    intialization();
    vector<int> size_avalanche;  
    vector<double> start_time;                              //starting time of avalenches
    vector<double> end_time;                                //ending time of avalenches
    vector<double> s_spike;                                 //time of spikes
    vector<int> s_ind;                                      //neuron index for spikes
    int num_ava=0;                                          //number of avalanches have been intialized

    auto start = chrono::high_resolution_clock::now();
    double *pmine, *pmins, *pminr;
    int eind, sind, rind;
    double emin=0.0, smin=0.0, rmin=0.0;
    
    double step = 1000;
    while(emin<stimuT)
    {

        pmine = min_element(begin(next_exteral),end(next_exteral));
        eind = distance(begin(next_exteral), pmine);
        emin = next_exteral[eind];
        pmins = min_element(begin(next_spiking),end(next_spiking));
        sind = distance(begin(next_spiking), pmins);
        smin = next_spiking[sind];
        pminr = min_element(begin(next_setting),end(next_setting));
        rind = distance(begin(next_setting), pminr);
        rmin = next_setting[rind];

        if ((rmin<smin)&&(rmin<emin)){     //reset
            double ccur = 0.0;
            double cvol = Vr;
            current[rind] = ccur;
            voltage[rind] = cvol;
            event_time[rind] = rmin;
            next_setting[rind] = stimuT;
        } else if (emin < smin){           //enact external stimuli and update the time for the next
            double pcur = current[eind];
            double pvol = voltage[eind];
            double ptime = event_time[eind];
            double delt = emin-ptime;
            assert(delt>=0.0);
            next_exteral[eind] += expdis(generator);
            double cvol = (pvol+fact*pcur)*exp(iRC*delt)-fact*pcur*exp(taus*delt);
            double ccur = pcur*exp(taus*delt)+we;
            double ntime = next_time(ccur, cvol);
            
            if (ntime != 0.0){
                next_spiking[eind] = emin+ntime;
            }else{
                next_spiking[eind] = stimuT; 
            }
            event_time[eind] = emin;
            current[eind] = ccur;
            voltage[eind] = cvol;
            external_time[eind] = emin;
            next_exteral[eind] += expdis(generator); 
            assert(cvol<=20.0);

        } else{                                  //enact a spike and related update
            num_firing[sind] += 1;
            Neuron* nptr = neurons[sind];        //pointer to current neuron
            int inum = nptr->num_in;
            
            if(sind>=s_neuron){
                int ctype = 0;
                double ctime = -2.0;
                double cspiking = 0.0;
                int *in_con = nptr->in_connection;
                int causal = -1;
                for (int ns=0; ns<inum; ns++){   //find the latest causal spiking
                    cspiking=spiking_time[*(in_con+ns)];
                    if(ctime<cspiking){
                        ctime=cspiking;
                        causal = *(in_con+ns);
                        ctype=avalanche_num[causal-s_neuron];
                    }
                }
                if (sind<t_neuron){
                    ctype=num_ava;
                    size_avalanche.push_back(0);
                    start_time.push_back(smin);
                    end_time.push_back(smin);
                    num_ava += 1;
                }else{
                    num_branching[causal] += 1;
                }
                avalanche_num[sind-s_neuron] = ctype;
                size_avalanche[ctype] += 1;
                end_time[ctype] = smin;
                s_ind.push_back(sind);
                s_spike.push_back(smin);
            }

            event_time[sind] = smin;                  
            current[sind] = 0.0;                                   // current at the time
            voltage[sind] = Vr;                                    // voltage resetting
            next_setting[sind] = smin+taurp;                       // the resetting time
            next_spiking[sind] = stimuT;
            
            int onum = nptr->num_out;
            double *utptr = nptr->synapse_util;
            int *out_con = nptr->out_connection;
            double pstime = spiking_time[sind];
            double sdelt = smin-pstime;
            assert(sdelt>0.0);
            spiking_time[sind] = smin;
            for (int ns=0; ns<onum; ns++){                         //update the neruons this neuron connneted to        
                int oind = *(out_con+ns);
                assert(oind != sind);
                double utill = 0.0;
                for (int nsy=0; nsy<rsN; nsy++){
                    double putill = *utptr;
                    double cutill = 1.0-(1.0-putill)*exp(taur*sdelt);
                    if (uniform(generator)<pr*cutill){
                        utill+=1.0;
                        *utptr = 0.0;
                    }else{*utptr=cutill;}
                    ++utptr;
                }
                double pcur = current[oind];
                double pvol = voltage[oind];
                double ptime = event_time[oind];
                double delt = smin-ptime;
                double cvol = (pvol+fact*pcur)*exp(iRC*delt)-fact*pcur*exp(taus*delt);
                double ccur = pcur*exp(taus*delt)+utill*win;
                double ntime = next_time(ccur, cvol);
                if (ntime != 0.0){
                    next_spiking[oind] = smin+ntime;
                } else {
                    next_spiking[oind] = stimuT;
                }
                event_time[oind] = smin;
                inner_time[oind] = smin;
                current[oind] = ccur;
                voltage[oind] = cvol;
            }
        }
        
        if (emin>step){
            cout<<"----------"<<emin<<endl;
            step *= 10.0;
        }
    }

    auto end = chrono::high_resolution_clock::now(); 
    double time_taken =  chrono::duration_cast<chrono::nanoseconds>(end - start).count(); 
    cout << time_taken << endl; 
    int asize = size_avalanche.size();
    int ssize = s_spike.size();

    int *causal_ava;
    causal_ava = new int[asize];
    double *ava_start;
    ava_start = new double[asize];
    double *ava_end;
    ava_end = new double[asize];
    double *spike_time;
    spike_time = new double[ssize];
    int *spike_ind;
    spike_ind = new int[ssize];

    ofstream num_fire("./num_spike.dat", ios::binary);
    ofstream num_branch("./num_branch.dat", ios::binary);
    ofstream num_avan("./num_ava.dat", ios::binary);
    ofstream startt("./ava_start.dat", ios::binary);
    ofstream endt("./ava_end.dat", ios::binary);
    ofstream times("./stimes.dat", ios::binary);
    ofstream fsind("./sind.dat", ios::binary);

    copy(size_avalanche.begin(), size_avalanche.end(), causal_ava);
    copy(start_time.begin(), start_time.end(), ava_start);
    copy(end_time.begin(), end_time.end(), ava_end);
    copy(s_spike.begin(), s_spike.end(), spike_time);
    copy(s_ind.begin(), s_ind.end(), spike_ind);

    num_fire.write((const char*)&num_firing, sizeof(int)*num_neuron);
    num_branch.write((const char*)&num_branching, sizeof(int)*num_neuron);
    num_avan.write((const char*)causal_ava, sizeof(int)*asize);
    startt.write((const char*)ava_start, sizeof(double)*asize);
    endt.write((const char*)ava_end, sizeof(double)*asize);
    times.write((const char*)spike_time, sizeof(double)*ssize);
    fsind.write((const char*)spike_ind, sizeof(int)*ssize);
    cout<<"Done"<< endl;
    return 0;
}