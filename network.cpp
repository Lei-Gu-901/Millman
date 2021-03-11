#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <fstream>
#include <iterator>
#include <iostream>
#include <assert.h>
#include <math.h>
#include "parameter.h"

using namespace std;
double num = 7.5;                 //mean num of connections for a neuron
int perlayer = num_perlayer;
double cutoff = (double)num/perlayer;
int num_com = perlayer*chight;
double ccutoff = (double)num/num_com;
double icutoff = 0.9*(double)num/num_com;

unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);
std::uniform_real_distribution<double> distribution(0.0,1.0);

int main()
{   
    //generate neurons index by 1) number of neurons in a layer, 2) number of the layer
    int positions[num_neuron][2];
    int counter = 0;
    for (int i0=0; i0<hight; ++i0){
        for (int i1=0; i1<num_perlayer; ++i1){
            positions[counter][0] = i0;
            positions[counter][1] = i1;
            ++counter;
        }
    }
    assert(counter == num_neuron);

    vector <int> ext_neuron;         //the neurons receiving external simulations
    vector <int> connect;            
    int con_inf[num_neuron];

    int border = hight-chight-1;
    size_t isize = sizeof(int);
    for (int nn=0;nn<num_neuron; nn++){
        int origin[2] = {positions[nn][0], positions[nn][1]};
        int nlayer = origin[0];
        if (nlayer == 0){ext_neuron.push_back(nn);}
        // vector<int> tempv;
        counter = 0;
        if (nlayer<border){
            for (int j=0; j<num_neuron;++j){
                int layer = positions[j][0];
                if(layer==nlayer+1){
                    double ran = distribution(generator);
                    if (ran < cutoff) {connect.push_back(j);
                                ++counter;}   
                }
            }
            con_inf[nn]=counter;
        }

        if (nlayer == border){
            for (int j=0; j<num_neuron;++j){
                int layer = positions[j][0];
                if(layer>border){
                    double ran = distribution(generator);
                    if (ran < ccutoff) {connect.push_back(j);
                                ++counter;}   
                }
            }
            con_inf[nn]=counter;
        }
        
        if (nlayer > border){
            for (int j=0; j<num_neuron;++j){
                int layer = positions[j][0];
                if(layer> border && j!=nn){
                    double ran = distribution(generator);
                    if (ran < icutoff) {connect.push_back(j);
                                ++counter;}   
                }
            }
            con_inf[nn]=counter;
        }
    }
    // the connection.dat contains all connection information
    // the con_info.dat tells which segment belongs to which neuron
    // extneuron.dat are those receiving external stimuli
    ofstream con("connection.dat", ios::binary);   
    ofstream inf("con_info.dat", ios::binary);
    ofstream ext("extneuron.dat", ios::binary);
    
    int csize = connect.size();
    int *conarray;
    conarray = new int[csize];
    int esize = ext_neuron.size();
    int *extneuron;
    extneuron = new int[esize];
    copy(connect.begin(), connect.end(), conarray);
    copy(ext_neuron.begin(), ext_neuron.end(), extneuron);
    inf.write((const char*)&con_inf, isize*num_neuron);
    con.write((const char*)conarray, isize*csize);
    ext.write((const char*)extneuron, isize*esize);
    cout << csize << endl;
    return 0;
}

