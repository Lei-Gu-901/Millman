class Neuron{
    public:
    int num_out;                   //number of neurons this neuron conneting to
    int num_in;                    //number of neruons that inject current to this neuron
    int *out_connection;           //array of indices of neurons the neuron connects to
    int *in_connection;            //array of indices of neurons connects to this neuron
    double *synapse_util;          //utility of the synapses
    Neuron(int x, int y, int oind[], int iind[], double util[]);
    ~Neuron();
};