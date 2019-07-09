//
//  Dipolar_MC.cpp
//  
//  Class file for Dipolar Monte Carlo Simulations.
//
//  The purpose of this class file is to create an object for
//  performing Monte Carlo simulations that are based on the
//  dipolar interaction energy as the main energy term
//
//  Created by Charudatta Phatak on 7/8/19.
//
//

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <utility>
#include <random>
#include <thread>
#include <mutex>
//#include "input.h"

using namespace std;

//--------------------------------------------------------------------------------------
// Class definition.
class Dipolar_MC {
    public:
    
    // parameters describing the lattice
    int n_isl = 6; //No. of islands in the lattice.
    vector<double> center_x; //X-center of islands
    vector<double> center_y; //Y-center of islands
    vector<double> angles; //Angle of magnetization or the island.
    vector<double> nn_inds; //List of nearest neighbor indices
    double max_nn_num = 0.0; //Max. number of nearest neighbors
    double max_nn_dist = 0.0; //Max. number of distance for nearest neighbors
    
    //parameters describing MC sim
    double mc_iters = 1000; //Total no. of MC iters
    double eq_iters = 0; //number of iters for equilibriation.
    double temp = 0; //Temp. for MC sim
    
    //other variables for storing thermodynamic parameters
    double lat_energy = 0; //Lattice Energy
    double lat_netmag = 0; //Lattice Net Magnetization.
    double sp_heat = 0; //Sp. heat
    double suscep = 0; //Magnetic Susceptibility.
    int n_lowaccept = 0; //Low energy accepted flips.
    int n_highaccept = 0; //High energy accepted flips.
    int n_noaccept = 0; //No accepted flips.
    
    //Member functions.
    double Calc_Latt_Energy(); //Calculates the total lattice energy
    double Calc_E_Dipole(int site); //Calculates the change in energy due to spin flip of location <site>.
    void MCMC_move(double Beta); // MC_move fnction that actually performs the random spin flip, calculates change in energy, updates the lattice information. It also calculates the thermodynamic parameters of sp_heat, and susceptibility.
    
    Dipolar_MC(int num_isl, vector<double> cen_x, vector<double> cen_y, vector<double> angs, vector<double> num_inds, double max_nn_num, double max_nn_dist);
    Dipolar_MC();
};

//--------------------------------------------------------------------------------------
//Constructor for the Dipolar_MC class.
//
// num_isl: number of islands.
// cens: centers of islands.
// angs: angles of the magnetization.
// num_inds: number of the nearest neighbor indices.
// max_nn_num: max. number of nearest neigbors to consider
// max_nn_dist: max. number of nearest neighbor distance.
//--------------------------------------------------------------------------------------

Dipolar_MC :: Dipolar_MC(int num_isl, vector<double> cen_x, vector<double> cen_y, vector<double> angs, vector<double> num_inds, double max_nn_num, double max_nn_dist){
    
    this->n_isl = num_isl;
    this->center_x = cen_x;
    this->center_y = cen_y;
    this->angles = angs;
    this->nn_inds = num_inds;
    this->max_nn_num = max_nn_num;
    this->max_nn_dist = max_nn_dist;

}

Dipolar_MC :: Dipolar_MC(){
    
    this->n_isl = 0;
    this->center_x.assign (0.0, 1.0);
    this->center_y.assign (0.0,1.0);
    this->angles.assign (0.0,1.0);
    this->nn_inds.assign (0.0,1.0);
    this->max_nn_num = 0.0;
    this->max_nn_dist = 0.0;
    
    cout<<"Blank Initiated. \n";
    
}

//**************************************************************************************
// Main function
//

int main(){
    Dipolar_MC run1 = Dipolar_MC();
}
