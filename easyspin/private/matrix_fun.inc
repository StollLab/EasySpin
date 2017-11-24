/*
// =======================
// Code adapted from mmx package on File Exchange
// http://www.mathworks.com/matlabcentral/fileexchange/37515-mmx-multithreaded-matrix-operations-on-n-d-matrices
// =======================

// ==================================================
// straightforward implementations of matrix multiply
// ==================================================
*/
void multAB(double *A, const int rA, const int cA, 
            double *B, const int cB,
            double *C) {
   int i, j, k;
   double *a, *c, tmp;  
   for( i=0; i<cB; i++ ){
      for( k=0; k<cA; k++ ){
         tmp = B[i*cA+k];
         a   = A + k*rA;
         c   = C + i*rA;
         for( j=0; j<rA; j++ ){
            c[j] += tmp * a[j];
         }
      }
   }   
}

/*
// =============================
// multiply:   C = op(A) * op(B)
// =============================
*/
void mulMatMat(double *A, const int rA, const int cA, 
               double *B, const int rB, const int cB,
               double *C) {
      multAB(A, rA, cA, B, cB, C);
}
