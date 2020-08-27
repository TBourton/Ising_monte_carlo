/**
* ISING MODEL OPENMP CODE
* Tom Bourton, Swansea University
* 701329@swansea.ac.uk
**/

// Libary Includes
#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

using namespace std;

// Function Prototypes
double F_D_Boltzmann_Weight(int a);
void F_V_Initalise_Lattice(void);
void F_V_Initialise_Spins(void);
void F_V_Update(void);
void F_V_Measurements(int a);
int NearestNeightbour(int a);
int Heathbath(int i, int j);
double rndnum(void);
void rini(int iran)


// Hash Defines
#define numBins 25
#define numThreads 2
#define globalSeed 1783
#define coldStart 0
#define hotStart 1

// Global Variables
double beta;
double Magn
double Energy; 
int globalLatticeSize; 
int latticeChunk;
int **spin;
int localProcessID; // Unique ID of each process
int numProcesses; //Total number of processes
int startType;

// Begin Program
int main(int argc, char *argv[])
{
 //Local Variables
 double boltzmannWeight[5];

 int numThermalisation;
 int numMeasurements;
 int localSeed;
 int i; //Counter
 float elapsed_time;
 clock_t start, finish;

 omp_set_num_threads(numThreads);
 
 if (numMeasurements%numbins!=0)
 {
   printf("Error: the number of measurements must be diviside by the number of bins\n");
   exit(1);
 }
 
 #pragma omp parallel
 {  
   //OpenMP setup
   numProcesses = omp_get_num_threads();
   localProcessID = omp_get_thread_num();
   
   //Read in initilisation parameters from command line arguments
   startType = coldStart;
   
   // Test of divisibility
   
   //Initialise Random Number Generator
   localSeed = globalSeed*(5+localProcessID); //Each Thread has it's own set of Random Numbers
   rini(localSeed);
   

   //Calculate Boltzmann Weighting Function
     for (i = 0; i<5; i++)
      {
        boltzmannWeight[i] =  F_D_Boltzmann_Weight(i);
      }

   //Divide the global lattice into segments so that each thread works 
   //on a segment of the global lattice
   latticeChunk = globalLatticeSize/numProcesses; 
   
   //Initalise the Lattice
   F_V_Initalise_Lattice();

   //Thermalise the latticea
   for (i=0; i=NumThermalisation; i++)
   {
     F_V_Update();
   }

  //Initialise the observables
  Energy = new double[numMeasurements];
  Magn = new double [numMeasurements];
  
  //Measure the system
  #pragma omp for schedule(static, latticeChunk)
  for (i=0; i < numMeasurements; i++)
  {
    Energy[i] = 0.0;
    Magn[i] = 0.0;
    
    F_V_Update();
    F_V_Measurements(i);
  }
 }
 
 return 0;
} 

//Calculates and returns the Boltzmann Weight for a given value of a
double F_D_Boltzmann_Weight(int a)
{
  int b;
  int i;
  double Weight;

  b = 2*(a-2);
  Weight = 1.0/(1.0 + exp(2*b*beta));
  return Weight;
} 
 
void F_V_Initalise_Lattice(void)
{ 
  //Counting Variables
  int a,b,i;

  //Each thread create some new 1D rows of empty spins
  spin = new int* [globalLatticeSize];
  
 
  //For each Row on the local lattice, populate the Row with empty spins
  //S.T. each thread has a (localLatticeSize + 2)x(globalLatticeSize) grid
  for (a=0; a<globalLatticeSize; a++)
  {
    spin[a] = new int [globalLatticeSize];
  }   
  
  #pragma omp parallel for schedule(static, latticeChunk) 
  for (a=0; a<globalLatticeSize; a++)
  {
    for (b=0; b<globalLatticeSize; b++)
    {
      //If Cold Start then all spins locked in up (+1) Position
      if (startType = coldStart)
      {
        spin[a][b] = 1;
      }
      
      //If Hotstart then spins are assigned randomly
      else if (startType = hotStart)
      {
        if (rndnum() < 0.5)
        {
          spin[a][b] = 1;
        }        
        
        else
        {
          spin[a][b] = -1;
        }
      }
      
      //Error in startType parameter
      else
      {
       printf("Unknown Start Type!!, define 0 for Cold Start, or 1 for a Hot Start\n");
       exit(1);
      }   

    }
  }

}

void F_V_Update(void)
{
  int a,b;

  #pragma omp parallel for schedule(static)
  for (a=0; a=globalLatticeSize; a++)
  {
    for (b=0; b=globalLatticeSize; b++)
    {
      spin[a,b] = Heatbath(a,b);
    }

  }
}

int Heatbath(int i, int j)
{
  int trialspin;
  int k,l;
  int nnSum;
  int ilow, itop, jright, jleft;
  int nnEnergy;
 
  nnSum = 0;

  //define nearest neighbours
  //remebering C takes (0,0) at the top left of array
  ilow = i + 1;
  itop = i - 1;
  jright = j + 1;
  jleft = j - 1;

  //define boundary conditions
  ilow = NearestNeighbour(ilow);
  itop = NearestNeighbour(itop)
  jright = NearestNeighbour(jright);
  jleft = NearestNeighbour(jleft);

  //Sum nearest neighbours
  nnSum = spin[ilow][j] + spin[itop][j] + spin[i][jleft] + spin[i][jright];
    
  nnEnergy = (nnSum/2) + 2;

  if (rndnum() < weight[nnEnergy])
  {
     trialspin = -1;
  }
  else 
  {
      trialspin =1;
  }
  return trialspin;

}

int NearestNeighbour(int a)
{
  if (a > globalLatticeSize) a = 0;
  if (a < 0) a = globalLatticeSize;
  return a;
}


void F_V_Measuments(int a)
{
  int i, j;
  int Volume;
  int ilow;
  int jright;

  Volume = globalLatticeSize * globalLatticeSize;

  for (i=0; i < globalLatticeSize; i++)
  {
    for (j = 0; j < globalLatticeSize; j++)
    {
      ilow = i + 1;
      jright = j + 1;
      ilow = NearestNeighbour(ilow);
      jright = NearestNeighbour(jright);
 
      Energy[a] = Energy[a] + (spin[i][j] * (spin[ilow][j] + spin[i][jright]));
      Magn[a] = Magn[a] + spin[i][j];
    }

  }
  
  Magn[a] = fabs(Magn[a]) / Volume;
  Energy[a] = Energy[a] / Volume;
  
}
