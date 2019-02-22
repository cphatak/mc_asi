//
//  afasi_mcsims.cpp
//
//
# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <string>
# include <algorithm>
# include <vector>
# include <chrono>
# include <random>

#define PI 3.14159265
#define MU0 4.0*PI*1e-7
#define kB 1.38e-23

using namespace std;
using namespace std::chrono;

//Function initializations to be used in the program.
int main ( int argc, char *argv[] );
int initialize_lattice ( int m, int n, double a, double s, double cenx[], double ceny[], double angles[] );
double calc_energy ( int n_isl, double *distmap, double cenx[], double ceny[], double magx[], double magy[], double max_nn_dist, double max_nn_num, double energy_array[] );
void timestamp ( );


//****************************************************************

int main ( int argc, char *argv[] )

//****************************************************************
//
// Purpose:
//
// MAIN is the main program for the AFASI lattice simulations
//
// Usage:
//
// afasi_mcsims m n iterations st_temp en_temp red_temp
//
// * M, N, the number of rows and columns in the lattice.
//   Each motif consists of 6 islands.
// * ITERATIONS, the number of iterations to perform per temp step.
// * ST_TEMP, the start temperature parameter for the MC simulation.
// * EN_TEMP, the end temperature parameter for the MC simulation.
// * RED_TEMP, the reduction factor for temperature parameter.
// * A is the lattice parameter
// * S is the island separation
// * SEED, a seed for the random number generator.
//
// Licensing:
//
// This code is distributed under the GNU GPL license.
//
// Modified:
//
// 17 February 2019
//
// Author:
//
// Charudatta Phatak (cd@anl.gov)
{
    // Get starting timepoint
    auto start = high_resolution_clock::now();
    
    int lattice;
    int iterations;
    int m;
    int n;
    int seed;
    double st_temp;
    double en_temp;
    double red_temp;
    double a;
    double s;
    
    
    // Calculated parameters
    int n_temp;
    int n_isl;
    
    timestamp ( );
    cout << "\n";
    cout << "AF_ASI MC sim\n";
    cout << "c++ version\n";
    //
    //  Get Input
    //
    if (1 < argc )
    {
        m = atoi ( argv[1] );
    }
    else
    {
        m = 3;
    }
    if (2 < argc)
    {
        n = atoi ( argv[2] );
    }
    else
    {
        n = 3;
    }
    if (3 < argc)
    {
        iterations = atoi ( argv[3] );
    }
    else
    {
        iterations = 1000;
    }
    if (4 < argc)
    {
        st_temp = atof ( argv[4] );
    }
    else
    {
        st_temp = 1000;
    }
    if (5 < argc)
    {
        en_temp = atof ( argv[5] );
    }
    else
    {
        en_temp = 1;
    }
    if (6 < argc)
    {
        red_temp = atof ( argv[6] );
    }
    else
    {
        red_temp = 0.95;
    }
    if (7 < argc)
    {
        a = atof ( argv[7] );
    }
    else
    {
        a = 350.0;
    }
    if (8 < argc)
    {
        s = atof ( argv[8] );
    }
    else
    {
        s = 120.0;
    }
    if (9 < argc)
    {
        seed = atoi ( argv[9] );
    }
    else
    {
        seed = 1123581119;
    }
    //
    // Compute number of temp. steps.
    //
    n_temp = static_cast<int>(log(en_temp/st_temp)/log(red_temp));
    //
    // Compute the total number of islands.
    //
    n_isl = m*n*6;
    double cenx[n_isl];
    double ceny[n_isl];
    double angles[n_isl];
    //
    // Set the max. distance and max. nearest neighbors
    //
    double max_nn_dist = 500.0;
    double max_nn_num = 9;
    
    //seed random number generator
    srand((unsigned)time(0));
    
    //
    // Print useful information for the MC simulations.
    //
    cout << "\n";
    cout << "Parameters that can be changed by passing arguments in the order \n";
    cout << "------------------------------------------- \n";
    cout << " The number of rows is M = " << m << "\n";
    cout << " The number of columns is N = " << n << "\n";
    cout << " The number of iterations is  = " << iterations << "\n";
    cout << " The start temperature parameter is = " << st_temp << "\n";
    cout << " The end temperature parameter is = " << en_temp << "\n";
    cout << " The reduction factor for temperature parameter is = " << red_temp << "\n";
    cout << " The lattice parameter is A = " << a << "\n";
    cout << " The separation distance is S = " << s << "\n";
    cout << " The random generator seed = " << seed << "\n";
    cout << "------------------------------------------- \n";
    cout << "FIXED/Derived PARAMETERS:\n";
    cout << " Maximum number of N. Neighbors = " << max_nn_num << "\n";
    cout << " Maximum neighbor distance is = " << max_nn_dist << "\n";
    cout << " The total number of islands is = " << n_isl << "\n";
    cout << " Total number of temperature steps is = " << n_temp << "\n";
    cout << "\n";
    //
    // Initialize the lattice
    //
    cout << "Initializing lattice..\n";
    lattice = initialize_lattice(m,n,a,s,cenx,ceny,angles);
    
    //verbose infor for initializing lattice
    //int ii;
    //for (ii = 0; ii < n_isl; ii++)
    //{
    //    cout << cenx[ii] << ", " << ceny[ii] << ", " << angles[ii] << "\n";
    //}
    
    //
    // Compute the distance matrix
    //
    cout << "Computing distance matrix..\n";
    double *distmap = new double[n_isl*n_isl];
    int i,j;
    for (i = 0; i < n_isl; i++)
    {
        for (j = 0; j < n_isl; j++)
        {
            distmap[i + j*n_isl] = sqrt(pow((cenx[i]-cenx[j]),2) + pow((ceny[i]-ceny[j]),2));
           // cout << distmap[i][j] << "\t";
        }
        //cout << "\n";
    }
    
    //
    // Define the magnetization arrays.
    //
    cout << "Computing Magnetization.. \n";
    double magx[n_isl];
    double magy[n_isl];
    int ran_mult = 1;
    
    for (i = 0; i < n_isl; i++)
    {
        double rr = rand() / double (RAND_MAX);
        if (rr < 0.5)
        {
            ran_mult = (-1);
        }
        else
        {
            ran_mult = 1;
        }
        magx[i] = cos(angles[i] * PI / 180) * ran_mult;
        magy[i] = sin(angles[i] * PI / 180) * ran_mult;
        //cout << magx[i] << "\t" << magy[i] << "\n";
    }
    
    //
    // Calculate the energy
    //
    double energy_array[n_isl];
    double curr_energy = 0;
    curr_energy = calc_energy(n_isl, (double *)distmap, cenx, ceny, magx, magy, max_nn_dist, max_nn_num, energy_array);
    cout << "Total energy: " << curr_energy << "\n";
    
    // Here we start the loop over temperature range.
    int temp_i = 0;
    int mc_i = 0;
    double temp = st_temp;
    
    cout << "Temperature: " << temp << "\n";
    
    //variables for storing information.
    double latt_energy[n_temp];
    double latt_spheat[n_temp];
    double latt_netmag[n_temp];
    double latt_susc[n_temp];
    double temp_arr[n_temp];
    
    //random number for partition function
    double low_bound = 0;
    double up_bound = 1;
    std::uniform_real_distribution<double> unif(low_bound,up_bound);
    std::default_random_engine re;
    
    //Multiplier for the energy value/partition function.
    double mult_fac = MU0 * 1e-9 * 10;
    
    for (temp_i = 0; temp_i < n_temp; temp_i ++)
    {
        //define variables for storing MC run info.
        double avg_en = 0;
        double avg_en2 = 0;
        double avg_mag = 0;
        double avg_mag2 = 0;
        long n_accept = 0;
        double new_energy = 0;
        
        //beta value for partition function.
        double beta_val = mult_fac/(kB * temp);
        
        //now loop over MC iters.
        for (mc_i = 0; mc_i < iterations; mc_i++)
        {
            //Loop over total number of islands.
            for (int isl = 0; isl < n_isl; isl++)
            {
                //Pick a random island and change its magnetization.
                int site = rand() % n_isl;
                
                //change the magnetization.
                magx[site] *= (-1);
                magy[site] *= (-1);
                
                //calculate the energy.
                new_energy = calc_energy(n_isl, (double *)distmap, cenx, ceny, magx, magy, max_nn_dist, max_nn_num, energy_array);
                
                //calculate the difference
                double dE = new_energy - curr_energy;
                
                //check if we accept it or not.
                if (dE < 0)
                {
                    n_accept += 1;
                    curr_energy = new_energy;
                }
                
                if (dE > 0)
                {
                    double a_ran_num = unif(re);
                    double part_fun = exp(-dE * beta_val);
                    if (a_ran_num < part_fun)
                    {
                        n_accept += 1;
                        curr_energy = new_energy;
                    }
                    else
                    {
                        magx[site] *= (-1);
                        magy[site] *= (-1);
                    }
                }
                
                //compute thermodynamic variables.
                avg_en += curr_energy;
                avg_en2 += (curr_energy * curr_energy);
                double mag_val = 0;
                for (int i = 0; i < n_isl; i++)
                {
                    mag_val += sqrt(magx[i]*magx[i] + magy[i]*magy[i]);
                }
                avg_mag += mag_val;
                avg_mag2 += (mag_val * mag_val);
                
            }
        }
        //Compute the data
        double cn = 1.0/(iterations * n_isl);
        latt_energy[temp_i] = avg_en * cn;
        latt_netmag[temp_i] = avg_mag * cn;
        latt_spheat[temp_i] = (avg_en2*cn - avg_en*avg_en*cn*cn)/(temp*temp);
        latt_susc[temp_i] = (avg_mag2*cn - avg_mag*avg_mag*cn*cn)/(temp);
        temp_arr[temp_i] = temp;
        //Output values.
        cout << temp << "\t" << latt_energy[temp_i] << "\t" << latt_netmag[temp_i] << "\t" << latt_spheat[temp_i] << "\t" << latt_susc[temp_i] << "\n";
        
        //Update temperature
        temp *= red_temp;
        
    }
    
    
    //end
    // Get ending timepoint
    auto stop = high_resolution_clock::now();
    
    // Get duration. Substart timepoints to
    // get durarion. To cast it to proper unit
    // use duration cast method
    auto duration = duration_cast<seconds>(stop - start);
    
    cout << "Time taken by function: " << duration.count() << " seconds" << "\n";
    delete[] distmap;
    return 0;
}
//
// -----------------From here on are all the functions for the main program ----------------
//
//***************************************************************************

int initialize_lattice ( int m, int n, double a, double s, double cenx[], double ceny[], double angles[] )

//***************************************************************************
//
// Purpose:
//
// INITIALIZE_LATTICE intializes the lattice for AFASI by computing the centers,
// and angles for each island, and returns a list of island numbers.
//
// Usage:
//
// initialize_lattice(m, n, a, s, *cenx, *ceny, *angles)
//
// *M is the number of rows
// *N is the number of columns
// *A is the lattice parameter
// *S is the island separation
// *CENX is an array for holding the center_x of islands
// *CENY is an array for holding the center_y of islands
// *ANGLES is an array for holding the angles of each island
//
// Licensing:
//
// This code is distributed under the GNU GPL license.
//
// Modified:
//
// 11 February 2019
//
// Author:
//
// Charudatta Phatak (cd@anl.gov)
{
    int i = 0;
    int j = 0;
    int count = 0;
    double n_isl = m*n*6;
    
    for (j = 0; j < n; j++)
    {
        for (i = 0; i < m; i++)
        {
            // horizontal islands
            angles[count] = 0;
            angles[count+1] = 0;
            cenx[count] = i*2*a - a/2+j*a;
            ceny[count] = j*sqrt(3)*a - a*sqrt(3)/4+s/2;
            cenx[count+1] = i*2*a - a/2+j*a;
            ceny[count+1] = j*sqrt(3)*a - a*sqrt(3)/4-s/2;
            //first set of rotated islands
            angles[count+2] = 120;
            angles[count+3] = 120;
            cenx[count+2] = i*2*a+a/2 + j*a + sqrt(3)/4*s;
            ceny[count+2] = j*sqrt(3)*a - a*sqrt(3)/4 + s/4;
            cenx[count+3] = i*2*a+a/2 + j*a - sqrt(3)/4*s;
            ceny[count+3] = j*sqrt(3)*a - a*sqrt(3)/4 - s/4;
            //second set of rotated islands
            angles[count+4] = 60;
            angles[count+5] = 60;
            cenx[count+4] = i*2*a + j*a - sqrt(3)*s/4;
            ceny[count+4] = j*sqrt(3)*a + a*sqrt(3)/4 + s/4;
            cenx[count+5] = i*2*a + j*a + sqrt(3)*s/4;
            ceny[count+5] = j*sqrt(3)*a + a*sqrt(3)/4 - s/4;
            count = count + 6;
        }
    }
    return 1;
}

//***************************************************************************

double calc_energy ( int n_isl, double *distmap, double cenx[], double ceny[], double magx[], double magy[], double max_nn_dist, double max_nn_num, double energy_array[] )

//***************************************************************************
//
// Purpose:
//
//  CALC_ENERGY calculates and returns an array consisting of dipolar energies for
//  for each island.
//
// Usage:
//
//  calc_energy(distmap, n_isl, cenx, ceny, magx, magy, max_nn_dist, max_nn_num, energy_array)
//
//  *DISTMAP - distance map calculated as the distance of each island from another
//  *N_ISL - Number of islands.
//  *CENX - Array consisting the x value of centers of each island
//  *CENY - Array consisting the y value of centers of each island
//  *MAGX - Array consisting the x component of magnetization of each island
//  *MAGY - Array consisting of the y component of magnetization of each island
//  *MAX_NN_NUM - Maximum number of nearest neighbors to be considered
//  *MAX_NN_DIST - Maximum number of nearest neighbor distance.
//  *ENERGY_ARRAY - Array consisting of energy of each island.
//
// Licensing:
//
// This code is distributed under the GNU GPL license.
//
// Modified:
//
// 17 February 2019
//
// Author:
//
// Charudatta Phatak (cd@anl.gov)
{
    //variables for counters
    int i;
    int j;
    int cnt;
    int arr_index;
    double dist_val;
    //variables for energy expression
    double si_sj;
    double r_ij;
    double si_rij;
    double sj_rji;
    //total energy variable
    double temp_energy = 0;
    double total_energy = 0;
    
    //defining vector array for sorting and indexing nearest neighbors.
    vector<vector<double> > nn_array;
    nn_array.resize(n_isl);
    for (i = 0; i < n_isl; i++)
    {
        nn_array[i].resize(2);
    }
    
    //loop over each island i, calculate nearest neighbors, sort them, and compute energy.
    for (int i = 0; i < n_isl; i++)
    {
        // get the array of neighbors
        for (int j = 0; j < n_isl; j++)
        {
            //nn_array[j][0] = distmap[i][j];
            nn_array[j][0] = distmap[i + j*n_isl];
            nn_array[j][1] = static_cast<double>(j);
        }
        
        //next we sort all the elements of the array
        std::sort(nn_array.begin(),nn_array.end());

        //loop over the max. nearest neighbors.
        for (cnt = 1; cnt < n_isl; cnt++)
        {
            if (cnt <= max_nn_num)
            {
                
                //assign the array index of the nearest neighbors
                arr_index = static_cast<int>(nn_array[cnt][1]);
                
                dist_val = nn_array[cnt][0];
                //cout << arr_index << "\t" << dist_val << "\n";
                
                if (dist_val <= max_nn_dist)
                {
                    //compute various terms of the energy expression
                    si_sj = magx[i]*magx[arr_index] + magy[i]*magy[arr_index];
                    r_ij = dist_val;
                    si_rij = (cenx[i]-cenx[arr_index])*magx[i] + (ceny[i]-ceny[arr_index])*magy[i];
                    sj_rji = (cenx[arr_index]-cenx[i])*magx[arr_index] + (ceny[arr_index]-ceny[i])*magy[arr_index];
                    temp_energy = (si_sj)/pow(r_ij,3) - (3.0 * si_rij * sj_rji)/pow(r_ij,5);
                    total_energy += temp_energy;
                    //cout << si_sj << ", " << r_ij << ", " << si_rij << ", " << sj_rji << ", " << temp_energy << ", " << total_energy << "\n";
                    energy_array[i] += temp_energy;
                }
  
            }
        }
        
    }
    //cout << "\n finishing \n";
    return total_energy/2.0;
}

//***************************************************************************

void timestamp ( )

//***************************************************************************
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40
    
    static char time_buffer[TIME_SIZE];
    const struct std::tm *tm_ptr;
    size_t len;
    std::time_t now;
    
    now = std::time ( NULL );
    tm_ptr = std::localtime ( &now );
    
    len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );
    
    std::cout << time_buffer << "\n";
    
    return;
# undef TIME_SIZE
}
