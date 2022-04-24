# Parallelization-of-Energy-Calculation-for-a-box-of-water-molecules
In this project I have parallelized the massive energy calculation using technologies like MPI and OpenMP. Detailed description in the pdf file.

# MPI (Message Passing Interface based approach)
Here the project is parallelized on two levels.
1. The reading of the massive data.
2. The Energy calculation loop.

The three files namely - waters_parallel_op1.c , waters_parallel_op2.c , waters_parallel_op3_main.c
corresponds to the three methods of reading the data. Detailed description in the pdf file named - MPI_parallel.pdf

# OpenMP based approach
Here I have parallelized the energy calculation part using Intel's OpenMP.
Detailed description in the pdf file named - openmp waters.pdf, code - waters_openmp.c

# Observations :
1. OpenMP comes with a ease of programming as compared to MPI Parallelization.
2. But the operations happens inside a black box with not much control over the 
parallelization. Also upon every run of the program we get different values for time and load 
imbalance.
3. Compared to the MPI Implementation, The wall time was low for the 'auto' scheduler setting 
in OpenMP but since the MPI parallel code was hardcoded, we had more control over the 
algorithm to minimize the Load Imbalance.
4. MPI provides us with much better control over the parallel environment that's why no 
setting of OpenMP was able to beat the load imbalance from the MPI parallelization.
5. The best setting with time in focus -> guided scheduler with chunk size 100 (This may vary 
for every execution because the number of npairs were different for different runs.)
6. The best setting with Load Imbalance in focus -> dynamic with chunk size 100.(This may vary 
too for every run.)
7. Unlike MPI where we were to obtain similar results for a particular setting, in OpenMP the 
results varied. Like, for the dynamic scheduler where the runtime figures out which thread is 
free and can take up the next chunk, we dont have control over this part. 
8. The overhead contribution is different for different schedulers and chunk sizes.
9. Too large chunk size may not justify the load balancing (thread with more work may get 
stuck while others wait.) and too small chunk size will contribute to larger overhead values for 
scheduling like dynamic
