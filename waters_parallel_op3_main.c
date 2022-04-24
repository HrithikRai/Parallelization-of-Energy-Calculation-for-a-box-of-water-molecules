// waters.c OPTION 3 - ONLY THE ROOT THREAD READS THE DATA AND THEN BROADCASTS IT TO ALL THE OTHER THREADS.
// All the energies calculated with the parallel code were matched against the one calculated with the serial code.



#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <time.h> // for CPU time
#include <sys/time.h> //for gettimeofday
#include <mpi.h> // Message Passing Interface for Parallelizing the Job
#include <string.h>
#define LENGTH 80

// global variables
const int maxnum=10000;
double r[maxnum][3][3],rcutsq=1.44,L;


// r(number of molecule, atom 0=O,1=H,2=H, coordinate 0=x,1=y,2=z)

double sqr(double a){return a*a;}

double energy12(int i1,int i2){
// ============================
  int m,n,xyz;
  double shift[3],dr[3],mn[3],r6,distsq,dist,ene=0;
  const double sig=0.3166,eps=0.65,eps0=8.85e-12,e=1.602e-19,Na=6.022e23,q[3]={-0.8476,0.4238,0.4238};
  double elst,sig6;
  elst=e*e/(4*3.141593*eps0*1e-9)*Na/1e3,sig6=pow(sig,6);

  // periodic boundary conditions
  for(xyz=0;xyz<=2;xyz++){
   dr[xyz]=r[i1][0][xyz]-r[i2][0][xyz];shift[xyz]=-L*floor(dr[xyz]/L+.5); //round dr[xyz]/L to nearest integer
   dr[xyz]=dr[xyz]+shift[xyz];
  }
  distsq=sqr(dr[0])+sqr(dr[1])+sqr(dr[2]);
  if(distsq<rcutsq){ // calculate energy if within cutoff
    r6=sig6/pow(distsq,3);
    ene=4*eps*r6*(r6-1.); // LJ energy
    for(m=0;m<=2;m++){
      for(n=0;n<=2;n++){
        for(xyz=0;xyz<=2;xyz++) mn[xyz]=r[i1][m][xyz]-r[i2][n][xyz]+shift[xyz];
        dist=sqrt(sqr(mn[0])+sqr(mn[1])+sqr(mn[2]));
        ene=ene+elst*q[m]*q[n]/dist;
    } }
  }
  return ene;
}

long *decider(const int num_procs,int j, long n){
    
    long *arr = malloc(sizeof(long)*num_procs);
    
    for(int i =0; i<num_procs; i++){
        
        long count = 0;
        long sum = 0;
        
        do{
            sum =  sum + j;
            j--;
            count++;
        }while(sum<n && j >0);
        
        arr[i] = count;
        
    }
    return arr;
}

long *decider2(int num_procs,int j, long n){
    
    
    long *arr2 = malloc(sizeof(long)*num_procs);
    long s = 0;

    for(int i =0; i<num_procs; i++){
        
        long count = 0;
        long sum = 0;
        
        do{
            sum =  sum + j;
            s++;
            j--;
            count++;
        }while(sum<n && j >0);
        
        
        arr2[i] = s;
        
    }
    return arr2;
}

int
main(int argc,char **argv){
  
  int num_procs,rank;
  // GOAL 1 - ADD MPI COMMANDS
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  struct timeval begin, end;
  gettimeofday(&begin, 0);


  int i,j,natoms,nmol; 
  bool case_zero;
  
  double energy=0,dtime;
  FILE *fp;
  char line[LENGTH],nothing[LENGTH],name[20];
  clock_t cputime0,cputime1,cputime2,cputime3,cputime4; /* clock_t defined in <time.h> and <sys/types.h> as int */

  // GOAL 2 - PARALLELIZE THE READING
  cputime0 = clock();                   // assign initial CPU time (IN CPU CLOCKS)
  if (rank==0){
      
        
        printf("Program to calculate energy of water\n");
        printf("Input NAME of configuration file\n");
        
        
        scanf("%s",name);                     // reading of filename from keyboard by thread 0
        fp=fopen(name, "r");                  // thread 0 opens the file and begins reading from HDD
        fgets(line, LENGTH,fp);               // skip first line
        fgets(line, LENGTH,fp); sscanf(line,"%i",&natoms);

        nmol=natoms/3; printf("Number of molecules %i\n",nmol);

        for (i=0;i<nmol;i++){
            for(j=0;j<=2;j++){
            fgets(line, LENGTH,fp);
            sscanf(line, "%s %s %s %lf %lf %lf",nothing,nothing,nothing, &r[i][j][0],&r[i][j][1],&r[i][j][2]);
        } }
        printf("first line %lf %lf %lf\n",r[0][0][0],r[0][0][1],r[0][0][2]);
        fscanf(fp, "%lf",&L);                      // read box size
        
        gettimeofday(&end, NULL);
        printf("Box size %lf\n",L);
        fclose(fp);
        cputime1= clock()-cputime0;              // calculate  cpu clock time as difference of times after-before
        printf("Time required to read data is %lf\n",(double)cputime1/CLOCKS_PER_SEC);
        gettimeofday(&end, 0);
        long thread_sec = end.tv_sec - begin.tv_sec;
        long thread_microsec = end.tv_usec - begin.tv_usec;
        double read_elapsed = thread_sec + thread_microsec*1e-6;
        
  }

  MPI_Bcast(&nmol,1,MPI_INT,0,MPI_COMM_WORLD);        // Broadcast the number of molecules to other threads
  MPI_Bcast(&L,1,MPI_DOUBLE,0,MPI_COMM_WORLD);        // Broadcast the Box size to other threads
  MPI_Bcast(&r,maxnum*9,MPI_DOUBLE,0,MPI_COMM_WORLD); // Broadcast the 3d array to other threads


  cputime1= clock()-cputime0;                       
  
  printf("first line read by thread %d is  %lf %lf %lf\n",rank,r[0][0][0],r[0][0][1],r[0][0][2]);
  printf("L for thread %d is %lf\n",rank,L);
  printf("Nmol for thread %d is %d\n",rank,nmol);


  long n = (long)((nmol)*(nmol-1))/(2*num_procs);
  long *arr = decider(num_procs,nmol-1,n);
  long *arr2 = decider2(num_procs,nmol-1,n);


  cputime2 = clock();                           // assign  CPU time (IN CPU CLOCKS) before Partial Energy Calculations

// GOAL 3 - PARALLELIZE THE ENERGY CALCULATION
  long npair = 0;
  for(i=arr2[rank]-arr[rank];i<arr2[rank];i++){ // calculate energy as sum over the pairs allocated to respective thread.
    for(j=i+1;j<nmol;j++){
      // GOAL 4 - keeping tracks of pairs formed for further reductions
      npair++;
      energy=energy+energy12(i,j); 
    } 
  }


  cputime3 = clock()-cputime2;            


  cputime4 = clock() - cputime0;

 
  // TIME REDUCTIONS
  int c1 = cputime1;
  int c3 = cputime3;
  int c4 = cputime4;
  int mintime1 = 0;
  MPI_Reduce(&c1,&mintime1,1,MPI_INT,MPI_MIN,0,MPI_COMM_WORLD);
  int mintime3 = 0;
  MPI_Reduce(&c3,&mintime3,1,MPI_INT,MPI_MIN,0,MPI_COMM_WORLD);
  int mintime4 = 0;
  MPI_Reduce(&c4,&mintime4,1,MPI_INT,MPI_MIN,0,MPI_COMM_WORLD);

  int maxtime1 = 0;
  MPI_Reduce(&c1,&maxtime1,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
  int maxtime3 = 0;
  MPI_Reduce(&c3,&maxtime3,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
  int maxtime4 = 0;
  MPI_Reduce(&c4,&maxtime4,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);


  double calc_walltime = 0;
  // MPI_Reduce(&calc_elapsed,&calc_walltime,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

  // GOAL 5 -PAIRS REDUCTION, Load imbalnce calculated below with the final results from the root thread
  long minnpair = 0;
  MPI_Reduce(&npair,&minnpair,1,MPI_INT,MPI_MIN,0,MPI_COMM_WORLD);
  long maxnpair = 0;
  MPI_Reduce(&npair,&maxnpair,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);

  // ENERGY REDUCTION
  double total_energy = 0;
  printf("Partial Energy Calculated by thread %d is %lf\n",rank,energy);
  MPI_Reduce(&energy, &total_energy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  

  if(rank==0){

    printf("Total energy %lf \n ",total_energy);
    printf("Energy per molecule %lf \n",total_energy/nmol);

    printf("Minimum number of pairs = %d , Maximum Number of pairs = %d\n",minnpair,maxnpair);
    printf("The Load Imbalance is %lf \n",2.*(maxnpair-minnpair)/(0.+maxnpair+minnpair));
    printf("Minimum Time taken for reading and broadcasting -> %lf \n ",(double)mintime1/CLOCKS_PER_SEC);
    printf("Minimum Time taken for Partial Energy Calculation -> %lf \n ",(double)mintime3/CLOCKS_PER_SEC);
    printf("Minimum Time taken from reading to end of energy Calculation -> %lf \n ",(double)mintime4/CLOCKS_PER_SEC);
    printf("Maximum Time taken for reading and broadcasting -> %lf \n ",(double)maxtime1/CLOCKS_PER_SEC);
    printf("Maximum Time taken for Partial Energy Calculation -> %lf \n ",(double)maxtime3/CLOCKS_PER_SEC);
    printf("Maximum Time taken from reading to end of energy Calculation -> %lf \n ",(double)maxtime4/CLOCKS_PER_SEC);
    
  }
  MPI_Barrier(MPI_COMM_WORLD);
  gettimeofday(&end, 0);
  long thread_sec = end.tv_sec - begin.tv_sec;
  long thread_microsec = end.tv_usec - begin.tv_usec;
  double thread_elapsed = thread_sec + thread_microsec*1e-6;
  
  printf("I'm thread %d and time measured: %.3f seconds.\n",rank, thread_elapsed);
  

  MPI_Finalize();

}




