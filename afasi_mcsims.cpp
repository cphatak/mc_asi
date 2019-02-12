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

using namespace std;

int main ( int argc, char *argv[] );
int initialize_lattice ( int m, int n, double a, double s, double cenx[], double ceny[], double angles[] );
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
// This code is distributed under the GNU LGPL license.
//
// Modified:
//
// 11 February 2019
//
// Author:
//
// Charudatta Phatak (cd@anl.gov)
{
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
    double n_temp;
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
        m = 5;
    }
    if (2 < argc)
    {
        n = atoi ( argv[2] );
    }
    else
    {
        n = 5;
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
        a = 350;
    }
    if (8 < argc)
    {
        s = atof ( argv[8] );
    }
    else
    {
        s = 120;
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
    n_temp = log(en_temp/st_temp)/log(red_temp);
//
// Compute the total number of islands.
//
    n_isl = m*n*6;
    double cenx[n_isl];
    double ceny[n_isl];
    double angles[n_isl];
    
    cout << "\n";
    cout << " The number of rows is M = " << m << "\n";
    cout << " The number of columns is N = " << n << "\n";
    cout << " The total number of islands is = " << n_isl << "\n";
    cout << " The number of iterations is  = " << iterations << "\n";
    cout << " The start temperature parameter is = " << st_temp << "\n";
    cout << " The end temperature parameter is = " << en_temp << "\n";
    cout << " The reduction factor for temperature parameter is = " << red_temp << "\n";
    cout << " Total number of temperature steps is = " << n_temp << "\n";
    cout << " The random generator seed = " << seed << "\n";
    cout << "\n";
//
// Initialize the lattice
//
    cout << "Initializing lattice..\n";
    lattice = initialize_lattice(m,n,a,s,cenx,ceny,angles);
    int ii;
    for (ii = 0; ii < n_isl; ii++)
    {
        cout << cenx[ii] << ", " << ceny[ii] << ", " << angles[ii] << "\n";
    }
//
// Compute the distance matrix
//
    double distmap[n_isl][n_isl];
    int i,j;
    for (i = 0; i < n_isl; i++)
    {
        for (j = 0; j < n_isl; j++)
        {
            distmap[i][j] = sqrt(pow((cenx[i]-cenx[j]),2) + pow((ceny[i]-ceny[j]),2));
            cout << distmap[i][j] << "\n";
        }
    }
    
    cout << "\n";
    timestamp ( );
    
    return 0;
}
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
// This code is distributed under the GNU LGPL license.
//
// Modified:
//
// 11 February 2019
//
// Author:
//
// Charudatta Phatak (cd@anl.gov)
{
    int *lattice;
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
        }
    }
    return 1;
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
