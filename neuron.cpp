#include "neuron.h"
#include "parameter.h"

Neuron::Neuron(int x, int y, int oind[], int iind[], double util[]): num_out(x), num_in(y), out_connection(new int[x]), in_connection(new int[y]), synapse_util(new double[x*rsN]){
            for (int i=0; i<num_out; i++) {
                out_connection[i] = oind[i];
            }
            for (int i=0; i<num_in; i++) {
                in_connection[i] = iind[i];
            }
            for (int i=0; i<num_out*rsN; i++) {
                synapse_util[i] = util[i];
            }
        }
Neuron::~Neuron() {delete[] synapse_util;
                   delete[] in_connection;
                   delete[] out_connection;
                   }