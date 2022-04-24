// waters.c OPTION 1 - SERIAL READING WHERE EVERY THREAD READS ITS OWN DATA



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
const int maxnum=100000;
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

long *decider(const int num_procs,long j, long n){
    
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

long *decider2(int num_procs,long j, long n){
    
    
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
  
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  int i,j,natoms,nmol; 
  bool case_zero;
  
  double energy=0,dtime;
  FILE *fp;
  char line[LENGTH],nothing[LENGTH],name[20];
  clock_t cputime0,cputime1,cputime2,cputime3,cputime4; /* clock_t defined in <time.h> and <sys/types.h> as int */
  
  // Calculating Reading time
  struct timeval begin, end;
  gettimeofday(&begin, 0);
  
  
  if (rank==0){
      
        printf("Program to calculate energy of water\n");
        printf("Input NAME of configuration file\n");
        scanf("%s",name); // reading of filename from keyboard by thread 0
        
  }
  cputime0 = clock();                             // assign  CPU time (IN CPU CLOCKS) for reading and broadcasting.
  //gettimeofday(&start, NULL);                     // returns structure with time in s and us (read_microsec)
  MPI_Bcast(&name,6,MPI_CHAR,0,MPI_COMM_WORLD);   // Broadcast the name of the file to other threads
  MPI_Barrier(MPI_COMM_WORLD);                    // Hold up until thread zero has broadcasted the name of the file to other threads...
  fp=fopen(name, "r");                            // opening of file and beginning of reading from HDD
  fgets(line, LENGTH,fp);                         // skip first line
  fgets(line, LENGTH,fp); sscanf(line,"%i",&natoms);

  nmol=natoms/3; printf("Number of molecules %i\n",nmol);

  for (i=0;i<nmol;i++){
    for(j=0;j<=2;j++){
      fgets(line, LENGTH,fp);
      sscanf(line, "%s %s %s %lf %lf %lf",nothing,nothing,nothing, &r[i][j][0],&r[i][j][1],&r[i][j][2]);
  } }
  printf("first line %lf %lf %lf\n",r[0][0][0],r[0][0][1],r[0][0][2]);
  fscanf(fp, "%lf",&L); // read box size
  
  // Stop measuring reading time and calculate the elapsed time
  gettimeofday(&end, 0);
  long read_seconds = end.tv_sec - begin.tv_sec;
  long read_microsec = end.tv_usec - begin.tv_usec;
  double read_elapsed = read_seconds + read_microsec*1e-6;

  printf("Reading time measured: %.3f seconds.\n", read_elapsed);

  printf("Box size %lf\n",L);

  long n = (nmol)*(nmol-1)/(2*num_procs);
  long *arr = decider(num_procs,nmol-1,n);
  long *arr2 = decider2(num_procs,nmol-1,n);

  printf("n = %d\n",n);
  for(int i =0 ; i < num_procs ; i++){
    printf("%d %d\n",arr[i],arr2[i]);
  }

  cputime2 = clock();                    // assign  CPU time (IN CPU CLOCKS) before Partial Energy Calculations
  //gettimeofday(&start, NULL);            // returns structure with time in s and us (read_microsec)
  long npair = 0;
  for(i=arr2[rank]-arr[rank];i<arr2[rank];i++){ // calculate energy as sum over the pairs allocated to respective thread.
    for(j=i+1;j<nmol;j++){
      npair++;
      energy=energy+energy12(i,j); 
    } 
  }
  cputime3 = clock()-cputime2;            // calculate  cpu clock time as difference of times after-before
  //gettimeofday(&end, NULL);


  cputime4 = clock() - cputime0;

 
  // TIME REDUCTIONS
  int c1 = cputime1;int c3 = cputime3;int c4 = cputime4;
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

  int alltime = 0;
  MPI_Reduce(&c4,&alltime,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);

  
  // PAIRS REDUCTION
  int minnpair = 0;
  MPI_Reduce(&npair,&minnpair,1,MPI_INT,MPI_MIN,0,MPI_COMM_WORLD);
  int maxnpair = 0;
  MPI_Reduce(&npair,&maxnpair,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);

  // ENERGY REDUCTION
  double total_energy = 0;
  MPI_Reduce(&energy, &total_energy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  

  if(rank==0){

    printf("Total energy %lf \n ",total_energy);
    printf("Energy per molecule %lf \n",total_energy/nmol);

    printf("Minimum number of pairs = %d , Maximum Number of pairs = %d\n",minnpair,maxnpair);
    printf("The Load Imbalance is %lf \n",(double)(maxnpair-minnpair)/((maxnpair+minnpair)/2));

    printf("Total CPU time taken by all threads -> %lf \n",(double)alltime/CLOCKS_PER_SEC);
    printf("Minimum Time taken for reading and broadcasting -> %lf \n ",(double)mintime1/CLOCKS_PER_SEC);
    printf("Minimum Time taken for Partial Energy Calculation -> %lf \n ",(double)mintime3/CLOCKS_PER_SEC);
    printf("Minimum Time taken from reading to end of energy Calculation -> %lf \n ",(double)mintime4/CLOCKS_PER_SEC);
    printf("Maximum Time taken for reading and broadcasting -> %lf \n ",(double)maxtime1/CLOCKS_PER_SEC);
    printf("Maximum Time taken for Partial Energy Calculation -> %lf \n ",(double)maxtime3/CLOCKS_PER_SEC);
    printf("Maximum Time taken from reading to end of energy Calculation -> %lf \n ",(double)maxtime4/CLOCKS_PER_SEC);
    
  }
  
  //dtime = ((end.tv_sec  - start.tv_sec)+(end.tv_usec - start.tv_usec)/1e6);
  //printf("Elapsed wall time for thread %d is: %f\n",rank, dtime);

  fclose(fp);

  MPI_Finalize();


}




