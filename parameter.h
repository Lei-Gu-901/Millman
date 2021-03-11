#ifndef PARAMETER_H
#define PARAMETER_H
    static const int stimuT = 1000000;                       //simulation time of 100s
    static const double we = 95*0.1/3.0;                     //pA  has /C
    static const double win = 50*0.1/3.0;                    //pA  has /C
    static const double threshold = 20.0;                    //mV  
    static const double Vr = 0.0;                            //mV
    static const double iRC = -0.05;                         // 1/RC ms
    static const int num_perlayer = 81;
    static const int hight = 16;                             //number of all layers
    static const int chight = 10;                            //number of compuation layers
    static const int num_neuron = num_perlayer*hight;        // number of neurons in this model
    static const int t_neuron = num_perlayer*(hight-chight); // locate the last transmission neuron
    static const int max_ind = num_neuron-1;
    static const int rsN = 6;                                 // number of relase site
    static const double edt = 40.0/1.0e3;                     // time step for external spiking
    static const double taus = -1.0/5.0;                      // current decay time 5ms
    static const double taur = -1.0/100.0;                    // synapse recover time 100ms
    static const int taurp = 1.0;                             // time steps for refraction
    static const double pr = 0.3;                             // rate for viscle release
#endif