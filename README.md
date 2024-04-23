The Laplace equation solved and its corresponding boundary conditions are from the 5-4 Application section (page 167) of "Computational Fluid Dynamics (Vol. 1)" 4th Edition by Klaus A. Hoffmann (Author), Steve T. Chiang. 

The equation is discretized on a staggered mesh (see section 6-2 : "Introduction to Computational Fluid Dynamics, An: The Finite Volume Method" 2nd Edition by H. Versteeg and W. Malalasekera) and different Krylov Space methods have been used to solve the resultant equation and their convergent rate is compared. 

The resultant matrix from the discretization of the equation is presented in the Compressed Sparse Row (CSR) format (see 3.4 Storage Schemes: "Iterative Methods for Sparse Linear Systems" 2nd Edition by Yousef Saad) for the efficiency in the computation of matrix-vector multiplication.
