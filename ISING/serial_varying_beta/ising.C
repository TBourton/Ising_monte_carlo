/**
* ISING MODEL CODE
* Tom Bourton, Swansea University
* 701329@swansea.ac.uk
**/

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>

#define numbins 25

using namespace std;

//prototypes
void read_input(void);
void print_log(void);
double rndnum(void);
void init_spin(void);
void rini(int iran);
void initweight(void);
int neighbour(int a, int b);
int staple(int a, int b);
int Heatbath(int a, int b);
void do_update(void);
void do_measurements(int a);
void write_measures();
void jack_error(double *obs, double *avg, double *err);
float elapsed_time;

// definition of variables
int **spin;
unsigned int thermalisation;
unsigned int size;
unsigned int measurements;
int randomseed=1812;
int istart;
double beta;
double weight[5];
double volume;
double *magn, *ene;
double avg_ene, err_ene;
double avg_magn, err_magn;
ofstream outfile("measures.dat"); 

// Start of the program
int main()
{
  clock_t start, finish;
  int l, m, n;
  double avg_ene, err_ene;
   
  start = clock();

  // Reads input and prints log
  read_input();
  print_log();

  //Sanity check
  if ((measurements/numbins*numbins) != measurements) {
    cout << "Error: the number of measurements must be divisible by numbins" << endl;
    cout << "We have numbins = " << numbins <<  " and " << measurements  << " requested measurements" << endl;
    exit(4);
  }

  // Lattice initialisation
  rini(randomseed);
  init_spin();
  initweight();

  // Thermalisation
  for ( l=0 ; l<thermalisation ; l++)
    do_update() ;

  // Initialisation of observables
  ene = new double [measurements]; 
  magn = new double [measurements]; 
  volume = (double)size*(double)size;

  // Measures
  for ( l=0 ; l<measurements ; l++) {
    do_update() ;
    do_measurements(l);
  }

  //Analysis
  jack_error(ene,&avg_ene,&err_ene);
  jack_error(magn,&avg_magn,&err_magn);
  cout << "Avg ene : " << avg_ene << " +/- " << err_ene << endl;
  cout << "Avg magn: " << avg_magn << " +/- " << err_magn << endl;

  // Save the measures on a file
  write_measures();
  
  finish = clock();
   
  elapsed_time = ((float)finish - (float)start) / (float)CLOCKS_PER_SEC;
  
  cout << "Computed in time " << elapsed_time << " seconds" << endl;

  return 1;
}

void read_input()
{
  cin >> size;
  cin >> beta;
  cin >> thermalisation;
  cin >> measurements;
  cin >> istart;
  if ( (istart != 1) && (istart != 0) ) {
    cout << "Error: istart can only be either 0 (cold) or 1 (hot)" << endl;
    exit(3);
  }  
}

void print_log()
{
  cout << "Ising model on a 2D square lattice" << endl ;
  cout << "Size of the lattice: " << size << endl ;
  cout << "Beta: " << beta << endl;
  cout << "Thermalisation sweeps: " << thermalisation << endl;
  cout << "Number of measurements: " << measurements << endl;
  if (istart == 1) {
    cout << "Start: hot" << endl;
  }
  else {
    cout << "Start: cold" << endl;
  }
}


void init_spin()
{
  int i,j;

  spin = new int* [size] ; 
  for ( i=0 ; i< size ; i++) {
    spin[i] = new int [size];
    for ( j=0 ; j < size ; j++ ) {
      switch (istart) {
      case 0: 
	spin[i][j] = 1 ;
	break;
      case 1:	
	spin[i][j] = rndnum() < 0.5 ? 1 : -1;
	break;
      default:
	cout << "Unknown initialisation parameter" << endl;
	exit(3);
	break;
      }
    }
  }
}

void initweight () {
  int i, j;
  
  for (i=0; i<5; i++) {
    j = 2*(i - 2) ;
    weight[i] = 1.0/(1.0 + exp(2*j*beta)) ;
  }
}

int neighbour(int i, int k) {
  int j;
  
  j = i - (i/k)*k ;
  if (j < 0) j = k + j; 

  return j;

}


int staple (int i, int j) {
  int k;
  
  int ifor = neighbour(i+1,size);
  int jfor = neighbour(j+1,size);
  int iback = neighbour(i-1,size);
  int jback = neighbour(j-1,size);
 
  k = spin[ifor][j] + spin[iback][j] + spin[i][jfor] + spin[i][jback];

  return k;
      
}

int Heatbath(int i, int j) {
  int k, m;
  int trialspin;
  double rran;

  
  k = staple(i,j);
  m = k/2 + 2; 

  rran = rndnum();

  if (rran < weight[m]) {
    trialspin = -1 ; 
  }
  else {
    trialspin = 1 ;
  }

  return trialspin;

}

void do_update() 
{
  int i,j;

  for ( i=0 ; i< size ; i++) {
    for ( j=0 ; j < size ; j++ ) {
      spin[i][j] =  Heatbath(i,j) ;
    }
  }
}

void do_measurements(int k)
{
  int i,j;

  magn[k] = 0.0;
  ene[k] = 0.0;
  
  for ( i=0 ; i< size ; i++) {
    for ( j=0 ; j < size ; j++ ) {
      int ifor = neighbour(i+1,size);
      int jfor = neighbour(j+1,size);
      magn[k] += spin[i][j] ;
      ene[k] += spin[i][j]*(spin[ifor][j] + spin[i][jfor]);
    }
  }

  magn[k] = fabs(magn[k])/volume ; 
  ene[k] /= volume ;
}

void write_measures()
{
  int l;

  for ( l=0 ; l < measurements ; l++) {
    outfile << ene[l] << "\t" << ene[l]*ene[l] << "\t" << \
      magn[l] << "\t" << magn[l]*magn[l] << "\n" ;
  }
}

void jack_error(double *obs, double *avg, double *err)
{
  double inv_volum = 1.0/volume;
  int l,k,slice;
  double bin[numbins], jackbins[numbins];
  double sumbins;

  *avg = 0.0;
  
  //Compute the simple average
  for ( l=0 ; l < measurements ; l++) {
    *avg += obs[l];
  }  

  *avg/=measurements;
  
  //Bin the data
  slice = measurements/numbins;
  sumbins = 0.0;
  for (l = 0; l< numbins; l++) {
    bin[l] = 0.0;
    for (k=0; k < slice; k++) {
      bin[l] += obs[l*slice + k] ; 
    }
    bin[l]/=slice;
    sumbins += bin[l];
  }
  //Form the jack-knife bins
  for (l = 0; l< numbins; l++) {
    jackbins[l] = (sumbins - bin[l])/(numbins - 1);
  }  
  
  //Compute the jack-knife error
  *err = 0.0;
  for (l = 0; l< numbins; l++) {
    *err += (*avg - jackbins[l])*(*avg - jackbins[l]);
  }
  *err *= ((double)(numbins - 1))/(double)numbins;
  *err = sqrt(*err);

}
