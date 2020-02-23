#include "Inversion.h"
#include <unistd.h> 
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#include <iostream>
#include <string.h>
#include <math.h>
#include "MICParaMgr.h"
#include "CovMat.h"
using namespace std;

extern "C" { void glasso_(int *, double *, double *, int *, int *, int *, int *, double *, int *, double *, double *, int *, double *, int *); }

Inversion::Inversion(JointMat* aData)
{
   mData=aData;
   mMargin=mData->margin();
   mJointProbMatrix=mData->joint();
   
   CovMat*          cov=new CovMat(aData);
   mCovMat=cov->getCovMat();
}



Inversion::~Inversion()
{

}

/*
void
Inversion::readCovMatrix(char* filename, int n)
{
  mCovMat=new double[n];
  ifstream input(filename);

  if(!input)
  {
    cerr<<"Can not open file "<<filename<<endl;
    exit(0);
  }

  char buf[100000];
  char* token;
  for(int i=0; i<n; i++)
  {
     input.getline(buf, 100000);
     token=strtok(buf, " \t\n");
     for(int j=0; j<n; j++)
     {
        mCovMat[i][j]=atof(token);   
        token=strtok(NULL, " \t\n");
     }
  }
}


void
Inversion::writeCovMatrix(char* filename, int n)
{
  ofstream output(filename);

  if(!output)
  {
    cerr<<"Can not open file "<<filename<<endl;
    exit(0);
  }

  char buf[100000];
  char* token;
  for(int i=0; i<n; i++)
  {
     for(int j=0; j<n; j++)
     {
        output<<mInvMat[i][j]<<" ";
     }
     output<<endl;
  }
}
*/

void
Inversion::inverse()
{
   if(MICParaMgr::_INVERSE_OPT_==_REMOVE_ONE_)
   {
     invertByRemovingLastColumns();
   }
   else if(MICParaMgr::_INVERSE_OPT_==_LASSO_)
   {
     invertByLasso();   
   }
   else if(MICParaMgr::_INVERSE_OPT_==_REGULARIZATION_)
   {
     invertByLSRegularize();
   }
}
void
Inversion::invertByRemovingLastColumns()
{//each column should have at least two columns

   //int size     = mData->getMaxHiddenSize();
   int tot_size = mData->numClus();
   int num_cols = mData->numCols();

   int non_red_col=tot_size-num_cols;

   double** reducedCovMat = new double* [non_red_col];

   for(int i=0; i<non_red_col; i++)
   {
      reducedCovMat[i]=new double[non_red_col];
   }//for

   for(int i=0; i<num_cols; i++)
   {
       int dim1=mData->dim(i);
       int off1=mData->offset(i);

       for(int j=0; j<num_cols; j++)
       {
         int dim2=mData->dim(j);
         int off2=mData->offset(j);

         for(int k1=0; k1<dim1-1; k1++)
         {
            for(int k2=0; k2<dim2-1; k2++)
            {
               reducedCovMat[off1-i+k1][off2-j+k2]=mCovMat[off1+k1][off2+k2];
            }
         }
       }
    }
    
    double** mat=inverseMat(reducedCovMat, non_red_col);


   for(int i=0; i<non_red_col; i++)
   {
      delete [] reducedCovMat[i];
   }//for
   delete [] reducedCovMat;


    mInvMat=new double *[tot_size];

    for(int i=0; i < tot_size; i++)
    {
       mInvMat[i]=new double[tot_size];
       for(int j=0; j< tot_size; j++)
       {
           mInvMat[i][j]=0;
       }
    }
    
    for(int i=0; i<num_cols; i++)
    {
       int dim1=mData->dim(i);
       int off1=mData->offset(i);

       for(int j=0; j<num_cols; j++)
       {
         int dim2=mData->dim(j);
         int off2=mData->offset(j);

         for(int k1=0; k1<dim1-1; k1++)
         {
            for(int k2=0; k2<dim2-1; k2++)
            {
               if(mInvMat[off1+k1][off2+k2]!=mInvMat[off1+k1][off2+k2])
                 cout<<"error in covarmatrix "<<endl;//=mJointProbMatrix[off1+k1][off2+k2]-mMargin[i][k1]*mMargin[j][k2]; 
                mInvMat[off1+k1][off2+k2]=mat[off1-i+k1][off2-j+k2];
            }
         }
         
       }
    }
}


void
Inversion::invertByLSRegularize()
{
   int tot_size = mData->numClus();
   int num_cols = mData->numCols();
   int non_red_col=tot_size;//-num_cols;

    
    double** mat=inverseByRegulizing(mCovMat, non_red_col, 0.1);//inverseMat(mCovMat, non_red_col);
    mInvMat=new double *[tot_size];

    for(int i=0; i < tot_size; i++)
    {
       mInvMat[i]=new double[tot_size];
       for(int j=0; j< tot_size; j++)
       {
           mInvMat[i][j]=0;
       }
    }
    
    for(int i=0; i<num_cols; i++)
    {
       int dim1=mData->dim(i);
       int off1=mData->offset(i);

       for(int j=0; j<num_cols; j++)
       {
         int dim2=mData->dim(j);
         int off2=mData->offset(j);

         for(int k1=0; k1<dim1; k1++)
         {
            for(int k2=0; k2<dim2; k2++)
            {
               mInvMat[off1+k1][off2+k2]=mat[off1+k1][off2+k2];
               
               if(mInvMat[off1+k1][off2+k2]!=mInvMat[off1+k1][off2+k2])
                 cout<<"error in covarmatrix "<<endl;//=mJointProbMatrix[off1+k1][off2+k2]-mMargin[i][k1]*mMargin[j][k2]; 
            }
         }
         
       }
    }
}



void
Inversion::invertByLasso()
{

    int ndim = mData->numClus();
    int seqlen = mData->numCols();

    
    double* rho =new double[ndim*ndim];//= f_matrix_calloc(ndim, sizeof(double));
    double* ww  =new double[ndim*ndim];// f_matrix_calloc(ndim, sizeof(double));
    double* wwi =new double[ndim*ndim];// f_matrix_calloc(ndim, sizeof(double));

    double lastfnzero=0.0;
    double rhodefault=-1;
    
    int nids, s, nseqs, ncon, opt, filtflg=0, approxflg=0, initflg=0, debugflg=0, diagpenflg=1;
    int apcflg=1, maxit=10000, npair, nnzero, niter, jerr, rawscflg = 1, pseudoc = 1, minseqsep = 5;
    unsigned int *wtcount;//ccount[MAXSEQLEN];
    double thresh=1e-4, del, **pcmat, *pcsum, pcmean, pc;
    double sum, score, (**pab)[150][150], **pa, wtsum, lambda,  fnzero, r2, targfnzero = 0.0, scsum, scsumsq, mean, sd, zscore, ppv;    
    double *weight, idthresh = -1.0, maxgapf = 0.9;

    
    double trialrho;
    // Guess at a reasonable starting rho value if undefined 
    if (rhodefault < 0.0)
	trialrho = maximum(0.001, 1.0 /1000);
    else
	trialrho = rhodefault;

    double rfact = 0.0;

    cerr<<"Entering invert by lasso"<<endl;
    bool shrinkflg=true;


    
    double* cmat=new double[ndim*ndim];
    double* tempmat=new double[ndim*ndim];

    for(int i=0; i<seqlen; i++)
    {
      int dim1=mData->dim(i);
      int off1=mData->offset(i);
        
      for(int j=0; j<seqlen; j++)
      {
         int dim2=mData->dim(j);
         int off2=mData->offset(j);
          
         for(int a=0; a<dim1; a++)
         {
            for(int b=0; b<dim2; b++)
            {
              cmat[(off1+a)*ndim+(off2+b)]=mJointProbMatrix[off1+a][off2+b]-mMargin[i][a]*mMargin[j][b];
            }
         }
      }
    }//for


    if (shrinkflg)
    {
        double smean=0;
        for (int i=0; i<ndim; i++)
	      smean += cmat[i*ndim+i];
	    smean /= (double)ndim;
        lambda = 0.1;

	  for (;;)
	  {
	     memcpy(tempmat, cmat, ndim*ndim*sizeof(double));
	    
	    /* Test if positive definite using Cholesky decomposition */
	    if (!test_cholesky(tempmat, ndim))
		    break;
	    
	    for (int i=0; i<seqlen; i++)
        {
          int dim1=mData->dim(i);
          int off1=mData->offset(i);
		  for (int j=0; j<seqlen; j++)
          {
            int dim2=mData->dim(j);
            int off2=mData->offset(j);
		    for (int a=0; a<dim1; a++)
            {
			   for (int b=0; b<dim2; b++)
               {
                  //if(i==j&& a==b) cerr<<cmat[(i*(ALPHABET_SIZE-1)+a)*ndim + j*(ALPHABET_SIZE-1)+b]<<" ";
			      if (i!= j || a!=b)
				     cmat[(off1+a)*ndim + off2+b] *= (1.0 - lambda);
			      else //if (a == b)
				     cmat[(off1+a)*ndim + off2+b] = smean * lambda + (1.0 - lambda) * cmat[(off1+a)*ndim + off2+b];
                  
                  //if(i==j && a==b) cerr<<cmat[(i*(ALPHABET_SIZE-1)+a)*ndim + j*(ALPHABET_SIZE-1)+b]<<" "<<endl;
               }
            }
          }
        }
        cerr<<"ha ha ha "<<smean<<endl;
	  }
    }
    



    for (;;)
    {
      cerr<<"entering the infinite loops"<<endl;
	  for (int i=0; i<ndim; i++)
	     for (int j=0; j<ndim; j++)
		 rho[i*ndim + j] = trialrho;
      for (int i=0; i<seqlen; i++)
      {
        int dim = mData->dim(i);
        int off = mData->offset(i);

        for (int a=0; a<dim; a++)
	 	    for (int b=0; b<dim; b++)
            {
               if(a!=b) rho[(off+a)*ndim + off +b] = 1e9;
            }

      }
	//for (int i=0; i<seqlen; i++)
	//    for (int j=0; j<seqlen; j++)
	//	for (int a=0; a<(ALPHABET_SIZE-1); a++)
	//	    for (int b=0; b<(ALPHABET_SIZE-1); b++)
	//		if ((a != b && i == j) || pa[i][ALPHABET_SIZE] > maxgapf || pa[j][ALPHABET_SIZE] > maxgapf)
	//		    rho[(i*ALPHABET_SIZE+a)*ndim + j*ALPHABET_SIZE+b] = 1e9;
	
    
	// All matrices are symmetric so no need to transpose before/after calling Fortran code 

    cerr<<"before call the graph lasso"<<endl;
	glasso_(&ndim, cmat, rho, &approxflg, &initflg, &debugflg, &diagpenflg, &thresh, &maxit, ww, wwi, &niter, &del, &jerr);
     cerr<<"calling the graph lasso"<<endl;
	if (targfnzero <= 0.0)
	    break;
    int i, j;
	for (npair=nnzero=i=0; i<ndim; i++)
	    for (j=i+1; j<ndim; j++,npair++)
		if (wwi[i*ndim+j] != 0.0)
		    nnzero++;

	fnzero = (double) nnzero / npair;

//      printf("rho=%f fnzero = %f\n", trialrho, fnzero);

	// Stop iterating if we have achieved the target sparsity level 
	if (fabs(fnzero - targfnzero)/targfnzero < 0.01)
	    break;
    cerr<<"the value at this iteration is "<<fabs(fnzero - targfnzero)/targfnzero<<endl;	
	if (fnzero == 0.0)
	{
	    // As we have guessed far too high, halve rho and try again 
	    trialrho *= 0.5;
	    continue;
	}
	
	if (lastfnzero > 0.0 && fnzero != lastfnzero)
	{
//	    printf("fnzero=%f lastfnzero=%f trialrho=%f oldtrialrho=%f\n", fnzero, lastfnzero, trialrho, trialrho/rfact);
	    
	    rfact = pow(rfact, log(targfnzero / fnzero) / log(fnzero / lastfnzero));

//	    printf("New rfact = %f\n", rfact);
	}

	lastfnzero = fnzero;

	// Make a small trial step in the appropriate direction 

	if (rfact == 0.0)
	    rfact = (fnzero < targfnzero) ? 0.9 : 1.1;
	
	trialrho *= rfact;
    }
    
    mInvMat=new double *[ndim];
    for(int i=0; i<ndim; i++)
    {
       mInvMat[i]=new double[ndim];
    }

    for(int i=0; i<seqlen; i++)
    {
       int dim1=mData->dim(i);
       int off1=mData->offset(i);
       for(int j=0; j<seqlen; j++)
       {
          int dim2=mData->dim(j);
          int off2=mData->offset(j);
           
          for(int a=0; a<dim1; a++)
          {
             for(int b=0; b<dim2; b++)
             {
               int i_prime=(off1+a);
               int j_prime=(off2+b);
               mInvMat[off1+a][off2+b]=wwi[(off1+a)*ndim+(off2+b)];
             }
          }
       }
    }
}

bool
Inversion::test_cholesky(double *a, const int n) 
{
    int i, j, k;
    double sum;
    static double *diag;

    if (diag == NULL)
	diag = new double[n];//(n, sizeof(double));

    for (i=0; i<n; i++)
    {
	  for (j=i; j<n; j++)
	  {
	    sum = a[i*n+j];

	    for (k=i-1; k >= 0; k--)
		   sum -= a[i*n+k]*a[j*n+k];

	    if (i == j)
	    {
		    cerr<<sum<<"in test cholestry"<<endl;
            if (sum <= 0.0)
		    return true;

		diag[i] = sqrt(sum);
	    }
	    else
		a[j*n+i] = sum / diag[i];
	  }
    }
    
    return false;
}

double**
Inversion::inverseMat(double** mat, int size)
{

	gsl_matrix * m = gsl_matrix_alloc (size, size);
	gsl_matrix * inverse = gsl_matrix_alloc (size, size);
	gsl_permutation * perm = gsl_permutation_alloc (size);
    int s;
	// Fill the matrix m
	for(int i=0; i<size; i++)
    {
       for(int j=0; j<size; j++)
       {
          gsl_matrix_set(m, i, j, mat[i][j]);
       }
    }

	gsl_linalg_LU_decomp (m, perm, &s);
	gsl_linalg_LU_invert (m, perm, inverse);
   
    double** rev=new double *[size];
    for(int i=0; i<size; i++)
    {
       rev[i]=new double[size];;
    }

  	for(int i=0; i<size; i++)
    {
       for(int j=0; j<size; j++)
       {
          rev[i][j]=gsl_matrix_get(inverse, i, j);
       }
    }
    gsl_matrix_free(m);
    gsl_matrix_free(inverse);
    gsl_permutation_free(perm);
    
    return rev;
}

void 
Inversion::cholesky(double **a, double **l, int n)
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
Inversion::inverseByRegulizing(double** mat, int size, double eta)
{

    double** beta   = new double* [size];
    double** dia    = new double* [size];
    double*  eye    = new double  [size];
    double*   ro    = new double  [size]; 
    double*   bo    = new double  [size]; 
    double** inverse= new double* [size]; 
    for(int i=0; i<size; i++) 
    {
        beta[i] =new double[size];
         dia[i] =new double[size];
     inverse[i] =new double[size];
    }

    for (int i=0; i<size; i++) 
    {
        double sum=0;
		for (int j=i; j<size; j++) 
        {
            if(i==j)
                sum=eta;
            else sum=0;

			for (int k=0; k<size; k++) 
            {
				sum += mat[i][k]*mat[k][j];
			}
			beta[i][j] = sum;
            beta[j][i] = sum;
		}		
	}		 			

	/** diagonalize matrix (cholesky) */
    cholesky(beta, dia, size);
    //printf ("#Cholesky decomposition done. Solving linear equations...\n");
    
    for (int i=0; i<size; i++) 
    { 	    	
       for (int s=0; s<size; s++)
       {
	      eye[s]=0;
	   }
	   eye[i]=1.0;  /* column of identity matrix*/
	    	
		   /* L*bo = e     ---- solve for b */
       int r;
       double* dvec;
       double xx;
       for (r=0; r<size; r++) 
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
	   for (r=size-1; r>=0; r--)
       {
           xx = bo[r];
           dvec = dia[r];
           for (int p=size-1; p>r; p--) 
           {
			 xx -= dia[p][r] * ro[p];  /*swap indices, as here L is transposed*/ 
		   }
		   ro[r] = xx/dvec[r];
		   inverse[i][r] = ro[r]; 
       }	  
    }

    //for(int i=0; i<size; i++)
    //{
    //   for(int j=0; j<size; j++)
    //       cout<<inverse[i][j]<<" ";
    //   cout<<endl;
   // }
    for(int i=0; i<size; i++) 
    {
      delete   beta[i];
      delete    dia[i];
    }

    delete [] beta;
    delete []  dia;
    delete []  eye;
    delete []   ro;
    delete []   bo;

    return inverse;
}




/*
double**
Inversion::inverseByEigen(int size)
{
  double data[] = { 1.0  , 1/2.0, 1/3.0, 1/4.0,
                    1/2.0, 1/3.0, 1/4.0, 1/5.0,
                    1/3.0, 1/4.0, 1/5.0, 1/6.0,
                    1/4.0, 1/5.0, 1/6.0, 1/7.0 };


  gsl_matrix * m = gsl_matrix_alloc (size, size);
  gsl_matrix * inverse = gsl_matrix_alloc (size, size);
  gsl_permutation * perm = gsl_permutation_alloc (size);

  gsl_matrix_view m 
    = gsl_matrix_view_array (data, 4, 4);

  gsl_vector *eval = gsl_vector_alloc (size);
  gsl_matrix *evec = gsl_matrix_alloc (size, size);

  gsl_eigen_symmv_workspace * w = 
    gsl_eigen_symmv_alloc (size);
  
  gsl_eigen_symmv (&m.matrix, eval, evec, w);

  gsl_eigen_symmv_free (w);

  gsl_eigen_symmv_sort (eval, evec, 
                        GSL_EIGEN_SORT_ABS_ASC);
  
  {
    int i;

    for (int i = 0; i < size ; i++)
    {
        double eval_i 
           = gsl_vector_get (eval, i);
        gsl_vector_view evec_i 
           = gsl_matrix_column (evec, i);

        printf ("eigenvalue = %g\n", eval_i);
        printf ("eigenvector = \n");
        gsl_vector_fprintf (stdout, 
                            &evec_i.vector, "%g");
      }
  }

  gsl_vector_free (eval);
  gsl_matrix_free (evec);


}*/

double**
Inversion::getInvCorr()
{
   
   if(MICParaMgr::_INV_CORR_OPT_==_INV_COV_MAT_FIRST_)
   {
     invertByLSRegularize();
     int tot_size = mData->numClus();
     mInvCorrMat=new double *[tot_size];
     for(int i=0; i<tot_size; i++) mInvCorrMat[i]=new double[tot_size];
     
     for(int i=0; i<tot_size; i++)
     { 
        for(int j=0; j<tot_size; j++)
        {
           mInvCorrMat[i][j]=mInvMat[i][j]/sqrt(mInvMat[i][i]*mInvMat[j][j]);
        }
     }//for

     //should delete the mInvMat here!!!!!!
   }
   else if(MICParaMgr::_INV_CORR_OPT_==_INV_CORR_MAT_FIRST_)
   {
      double** cov_saved=mCovMat; //save the mCotMat first
      int tot_size = mData->numClus();

      mCovMat    =new double *[tot_size];
      for(int i=0; i<tot_size; i++) 
      {
        mCovMat[i]    =new double[tot_size];
      }

      for(int i=0; i<tot_size; i++)
      { 
        for(int j=0; j<tot_size; i++)
        {
           mCovMat[i][j]=cov_saved[i][j]/sqrt(cov_saved[i][i]*cov_saved[j][j]);
        }
        invertByLSRegularize();
        mInvCorrMat=mInvMat;
      }//for

      for(int i=0;  i<tot_size; i++) delete [] mCovMat[i];
      delete [] mCovMat;
      mCovMat=cov_saved;
   }
   return mInvCorrMat;
}
