/**
* ISING MODEL MPI CODE
* Tom Bourton, Swansea University
* 701329@swansea.ac.uk
**/

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define numbins 25

using namespace std;

//prototypes
void read_input(void);
void print_log(void);
double rndnum(void);
void init_spin_MPI(void);
void rini(int iran);
void initweight(void);
int neighbour(int a, int b);
int staple(int a, int b);
int Heatbath(int a, int b);
void do_update_bulk_MPI(void);
void do_update_row1_MPI(void);
void do_updatetot_MPI(int *id);
void do_measurements_MPI(int a);
void write_measures_MPI();
void jack_error_MPI(double *observable, double *average, double *error);

// definition of variables
int **spin;
unsigned int thermalisation;
unsigned int size, sizenp,sizenp1, sizenp2;
unsigned int measurements;
int randomseed=1812, seedid;
int istart;
double beta;
double weight[5];
double volume, volumenp;
double *magn, *ene, *magntot, *enetot;
double average_energy, average_magn;
double error_energy, error_magn; 
int numprocs;
float elapsed_time;
ofstream outfile("measures.dat");

// Start of the program
int main(int argc, char *argv[])
{
  unsigned int l;
  int myid;
  clock_t start, finish;

  start = clock();

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  // Reads input and prints log
  if(myid==0){
    read_input();
    print_log();
  }
  
  MPI_Bcast(&size, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  MPI_Bcast(&beta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&thermalisation, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  MPI_Bcast(&measurements, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  
  sizenp=size/numprocs;
  sizenp1=sizenp+1;
  sizenp2=sizenp+2;

  //test
  if((size/numprocs)*numprocs-size !=0){
    printf("Error! size has to be multiple of numprocs\n");
    exit(11);
  }

  // Lattice initialisation
  seedid=randomseed*(1 + myid); //new seed for every myid
  rini(seedid);

  if(myid==0){
    initweight();
  }
  MPI_Bcast(weight, 5, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  init_spin_MPI();

  // Thermalisation
  for ( l=0 ; l<thermalisation ; l++){
    do_updatetot_MPI(&myid);
  }


  // Initialisation of observables
  ene = new double [measurements]; 
  enetot = new double [measurements]; 
  magn = new double [measurements];
  magntot = new double [measurements];
 
  volume = (double)size*(double)size;
  volumenp=(double)size*(double)sizenp;

  
  // Measures
  for ( l=0 ; l<measurements ; l++) {
    do_updatetot_MPI(&myid);
    do_measurements_MPI(l);
  }
  
  MPI_Reduce(magn, magntot, measurements, MPI_DOUBLE, MPI_SUM, 0, 
	     MPI_COMM_WORLD);

  MPI_Reduce(ene, enetot, measurements, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  // Save the measures on a file
  if(myid==0){
    write_measures_MPI();
  }

  //Analysis; we can perform analysis in serial as there is no benefit to running in parallel
 
 if (myid ==0)
 {
  jack_error_MPI(magntot, &average_magn, &error_magn);
  jack_error_MPI(enetot, &average_energy, &error_energy);

  cout << "The Average Energy is " << average_energy << " +/- " << error_energy << endl;
  cout << "The Average Magnetisation is " << average_magn << " +/- " << error_magn << endl;
 }

  MPI_Finalize();
  
 if (myid ==0)
 {
  finish = clock();
  elapsed_time = ((float)finish - (float)start)/(float)CLOCKS_PER_SEC;
  cout << "Computed in time " << elapsed_time << " seconds" << endl;
 }
  return 0;
}

void jack_error_MPI(double *observable, double *average, double *error)
{
   int l,i; //counter	
   int binwidth; //the width of the bins
   double bin[numbins]; 
   double jackbin[numbins]; // jacknife bins
   double sum_bins; //sum of all bins
   
  
   //Compute simple average
   *average = 0.0;

   for (l=0; l<measurements; l++)
   {
    *average = *average + observable[l];
   }
   *average= *average/measurements;

  //Bin the data serially
  binwidth = measurements/numbins; //define # of data points per bin
  sum_bins = 0.0;
  for (l=0; l<numbins; l++)
  {
    bin[l] = 0.0;
    for (i=0; i<binwidth; i++)
    {
       bin[l] = bin[l] + observable[l*binwidth + i];
    }

  bin[l] = bin[l]/binwidth; //Take partial average
  sum_bins =  sum_bins + bin[l]; 
  }
  
  //Jackknife bins
  for (i=0; i<numbins; i++)
  {  
    jackbin[i] = (sum_bins - bin[l])/(numbins - 1);
  }
  
  //Calculate error
  *error = 0.0;
  for (i=0; i<numbins;i++)
  {
    *error = *error + ((*average - jackbin[i])*(*average - jackbin[i]));
  } 
  
  *error = *error * (((double)(numbins - 1)) / ((double)numbins));
  *error = sqrt(*error);

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


void init_spin_MPI()
{
  unsigned int i,j;

  spin = new int* [sizenp2] ; 
  //  for ( i=1 ; i<= sizenp ; i++) {
  for ( i=0 ; i< sizenp2 ; i++) {
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

  i=0; //buffer i=0
  for ( j=0 ; j < size ; j++ ) {
    spin[i][j] = 0 ;
  }
  
  i=sizenp1; //buffer i=sizenp+1
  for ( j=0 ; j < size ; j++ ) {
    spin[i][j] = 0 ;
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


int staple_MPI (int i, int j) {
  int k;
  
  //  int ifor = neighbour(i+1,size);
  int ifor = i+1;
  int jfor = neighbour(j+1,size);
  //  int iback = neighbour(i-1,size);
  int iback = i-1;
  int jback = neighbour(j-1,size);
 
  if(i==0 || i==(signed)sizenp1){
    printf("Error in staple_MPI \n");
    exit(11);
  }

  k = spin[ifor][j] + spin[iback][j] + spin[i][jfor] + spin[i][jback];

  return k;
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
  unsigned int i,j;

  for ( i=0 ; i< size ; i++) {
    for ( j=0 ; j < size ; j++ ) {
      spin[i][j] =  Heatbath(i,j) ;
    }
  }
}

void do_update_bulk_MPI() 
{
  unsigned int i,j;

  for ( i=2 ; i<= sizenp ; i++) {
    for ( j=0 ; j < size ; j++ ) {
      spin[i][j] =  Heatbath(i,j) ;
    }
  }
}


void do_update_row1_MPI() 
{
  unsigned int i,j;

  i=1;
  for ( j=0 ; j < size ; j++ ) {
    spin[i][j] =  Heatbath(i,j) ;
  }
}


void do_measurements(int k)
{
  unsigned int i,j;

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

void do_measurements_MPI(int k)
{
  unsigned int i,j;

  magn[k] = 0.0;
  ene[k] = 0.0;
  
  for ( i=1 ; i<= sizenp ; i++) {
    for ( j=0 ; j < size ; j++ ) {
      //      int ifor = neighbour(i+1,size);
      int ifor = i+1;
      int jfor = neighbour(j+1,size);
      magn[k] += spin[i][j] ;
      ene[k] += spin[i][j]*(spin[ifor][j] + spin[i][jfor]);
    }
  }

  magn[k] = fabs(magn[k])/volumenp ;
  magn[k] /=numprocs;
  ene[k] /= volumenp ;
  ene[k] /= numprocs;
}

void write_measures()
{
  unsigned int l;

  for ( l=0 ; l < measurements ; l++) {
    outfile << ene[l] << "\t" << ene[l]*ene[l] << "\t" << \
      magn[l] << "\t" << magn[l]*magn[l] << "\n" ;
  }
}

void write_measures_MPI()
{
  unsigned int l;

  for ( l=0 ; l < measurements ; l++) {
    outfile << enetot[l] << "\t" << enetot[l]*enetot[l] << "\t" << \
      magntot[l] << "\t" << magntot[l]*magntot[l] << "\n" ;
  }
}


void do_updatetot_MPI(int *id)
{
  unsigned int j;
  int up_id, dn_id, myparity;
  MPI_Status status;
  
  up_id=(*id+1)%numprocs;
  dn_id=(*id+numprocs-1)%numprocs;
  
  myparity=(*id)%2;
  
  if(numprocs==1){
    myparity=3;
  }
  
  
  if(myparity==0){
    MPI_Send(spin[1], size, MPI_INT,dn_id,10, MPI_COMM_WORLD);
    MPI_Recv(spin[sizenp1], size, MPI_INT,up_id,10, MPI_COMM_WORLD, &status);
  }
  else if (myparity==1){
    MPI_Recv(spin[sizenp1], size, MPI_INT,up_id,10, MPI_COMM_WORLD, &status);
    MPI_Send(spin[1], size, MPI_INT,dn_id,10, MPI_COMM_WORLD);
  }
  else if (myparity==3){
    for ( j=0 ; j < size ; j++ ){
      spin[sizenp1][j]=spin[1][j];
      //      spin[0][j]=spin[sizenp][j];
    }
  }
  
  do_update_bulk_MPI();
  
  if(myparity==0){
    MPI_Send(spin[sizenp], size, MPI_INT,up_id,11, MPI_COMM_WORLD);
    MPI_Recv(spin[0], size, MPI_INT,dn_id,11, MPI_COMM_WORLD, &status);
  }
  else if (myparity==1){
    MPI_Recv(spin[0], size, MPI_INT,dn_id,11, MPI_COMM_WORLD, &status);
    MPI_Send(spin[sizenp], size, MPI_INT,up_id,11, MPI_COMM_WORLD);
  }
  else if (myparity==3){
    for ( j=0 ; j < size ; j++ ){
      //      spin[sizenp1][j]=spin[1][j];
      spin[0][j]=spin[sizenp][j];
    }
  }
  
  do_update_row1_MPI();
  
}
//*************************************************************************    
