#include "InvMat.h"
#include <unistd.h> 
#include <iostream>
#include <string.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
using namespace std;


InvMat::InvMat(int n)
{
   mDim=n;
   mCovMat=new double*[n];
   for(int i=0; i<n; i++)
       mCovMat[i]=new double[n];
 
   mInvMat=new double*[n];
   for(int i=0; i<n; i++)
       mInvMat[i]=new double[n];
}
   
InvMat::~InvMat()
{

}


void
InvMat::eigenDecompose(double** inputMat, double** eigenVec, double* eigenVal)
{

  gsl_matrix * m = gsl_matrix_alloc (mDim, mDim);
  
  gsl_vector *eval = gsl_vector_alloc (mDim);
  gsl_matrix *evec = gsl_matrix_alloc (mDim, mDim);
  
  gsl_eigen_symmv_workspace * w = 
    gsl_eigen_symmv_alloc (mDim);

  for (int i = 0; i < mDim; i++)
    for (int j = 0; j < mDim; j++)
      gsl_matrix_set (m, i, j, inputMat[i][j]);
  
  gsl_eigen_symmv (m, eval, evec, w); 
  
  
  for(int i=0; i<mDim; i++)
  {
     eigenVal[i]=gsl_vector_get (eval, i);
     for(int j=0; j<mDim; j++)
     {
        eigenVec[i][j]=gsl_matrix_get (evec, i, j);
     }
  }
  
  gsl_vector_free (eval);
  gsl_matrix_free (evec);
  gsl_eigen_symmv_free (w);
  gsl_matrix_free (m);
}

void
InvMat::deconvolution(double** inputMat)
{
    double** e_vec   = new double* [mDim];
    double*  e_val   = new double [mDim];
    for(int i=0; i<mDim; i++) 
    {
        e_vec[i] =new double[mDim];
    }
    eigenDecompose(mCovMat, e_vec, e_val);
    
    for(int i=0; i<mDim; i++)
    {
       for(int j=0; j<mDim; j++)
       {
           mInvMat[i][j]=0;
           for(int k=0; k<mDim; k++)
           {
             mInvMat[i][j]+=(e_val[k]/(1+e_val[k]))*e_vec[i][k]*e_vec[j][k];
           }
       }
    }
}
void
InvMat::read(char* filename)
{
  ifstream input(filename);

  if(!input)
  {
    cerr<<"Can not open file "<<filename<<endl;
    exit(0);
  }

  char buf[100000];
  char* token;
  for(int i=0; i<mDim; i++)
  {
     input.getline(buf, 100000);
     token=strtok(buf, ", \t\n");
     for(int j=0; j<mDim; j++)
     {
        mCovMat[i][j]=atof(token);   
        token=strtok(NULL, ", \t\n");
     }
  }
}


void
InvMat::write(char* filename)
{
  ofstream output(filename);

  if(!output)
  {
    cerr<<"Can not open file "<<filename<<endl;
    exit(0);
  }

  char buf[100000];
  char* token;
  for(int i=0; i<mDim; i++)
  {
     for(int j=0; j<mDim; j++)
     {
        output<<mInvMat[i][j]<<" ";
     }
     output<<endl;
  }
}


void
InvMat::inverse(double eta)
{
    mEta=eta;
    inverseByRegulizing(mCovMat);
}


void 
InvMat::cholesky(double **a, double **l, int n)
{
  /* Cholesky-Banachiewicz algorithm to make A triangular (A = L*Lt)*/ 
  double sum, *dd;

  dd = new double[n];

  for(int i=0;i<n;i++) {

     for(int j=i;j<n;j++) {

       /* get sum on and off-diagonal */
        sum = a[i][j]; 
        for (int k=0; k<i; k++) {
	         sum -= l[i][k] * l[j][k];
	    }

        if(i == j) {
	       if(sum <= 0.0) {
	          printf ("Error: matrix is not positive definite\n");
	           exit(0); 
	        }
	        dd[i] = sqrt(sum);
        }
        else {
	        l[j][i] = sum/dd[i];
        } 
     }
  } 
  /* finally save the diagonal */
  for(int i=0;i<n;i++)
       l[i][i] = dd[i];
  delete [] dd;
}



double**
InvMat::inverseByRegulizing(double** mat)
{

    double** beta   = new double* [mDim];
    double** dia    = new double* [mDim];
    double*  eye    = new double  [mDim];
    double*   ro    = new double  [mDim]; 
    double*   bo    = new double  [mDim]; 
    double** inverse= new double* [mDim]; 
    for(int i=0; i<mDim; i++) 
    {
        beta[i] =new double[mDim];
         dia[i] =new double[mDim];
     inverse[i] =new double[mDim];
    }

    for (int i=0; i<mDim; i++) 
    {
        double sum=0;
		for (int j=i; j<mDim; j++) 
        {
            
            
            if(i==j)
                sum=mEta;
            else sum=0;

			for (int k=0; k<mDim; k++) 
            {
				sum += mat[i][k]*mat[k][j];
			}
			beta[i][j] = sum;
            beta[j][i] = sum;
		}		
	}		 			

	/** diagonalize matrix (cholesky) */
    cholesky(beta, dia, mDim);
    printf ("#Cholesky decomposition done. Solving linear equations...\n");
    
    for (int i=0; i<mDim; i++) 
    { 	    	
       for (int s=0; s<mDim; s++)
       {
	      eye[s]=0;
	   }
	   eye[i]=1.0;  /* column of identity matrix*/
	    	
		   /* L*bo = e     ---- solve for b */
       int r;
       double* dvec;
       double xx;
       for (r=0; r<mDim; r++) 
       {
           xx = eye[r];
           dvec = dia[r];
           for (int p=0; p<r; p++) 
           {
               xx -= dvec[p] * bo[p];
           }
           bo[r] = xx/dvec[r];
       }     
		   	  
		   /* Lt*ro = bo     ---- solve for w */
	   for (r=mDim-1; r>=0; r--)
       {
           xx = bo[r];
           dvec = dia[r];
           for (int p=mDim-1; p>r; p--) 
           {
			 xx -= dia[p][r] * ro[p];  /*swap indices, as here L is transposed*/ 
		   }
		   ro[r] = xx/dvec[r];
		   mInvMat[i][r] = ro[r]; 
       }	  
    }

    //for(int i=0; i<mDim; i++)
    //{
    //   for(int j=0; j<mDim; j++)
    //       cout<<inverse[i][j]<<" ";
    //   cout<<endl;
   // }
    for(int i=0; i<mDim; i++) 
    {
      delete   beta[i];
      delete    dia[i];
    }

    delete [] beta;
    delete []  dia;
    delete []  eye;
    delete []   ro;
    delete []   bo;

    return mInvMat;
}


