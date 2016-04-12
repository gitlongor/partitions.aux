#include "config.h"

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>

extern "C" {
	/*
		Zoghbi and stojmenovic. Fast Algorithms for Generating Integer Partitions.
	*/
	SEXP ZS1(SEXP nR, SEXP startxR, SEXP kR, SEXP outR, SEXP nsolsR, SEXP min1stR)
	{
		int n; 	n=*INTEGER(nR);
		int m;  m=LENGTH(startxR);
		int k;  k=*INTEGER(kR);
		int h=1;
		int *x0; x0=INTEGER(startxR);
		int *out; out=INTEGER(outR);
		unsigned int nsols; nsols=(unsigned int)*INTEGER(nsolsR);
		unsigned int min1st; min1st=(unsigned int)*INTEGER(min1stR);
		unsigned int nsol;
		unsigned int * x;
		x = new unsigned int[n+1];

		for(int i=0; i<m; ++i) {
			*(out++) = x[i+1] = x0[i];
			if(x0[i]>1) h=i+1;
		}
		for(int i=m; i<k; ++i) {
			x[i+1]=1;
			*(out++)=0;
		}

		int r,t;
		nsol=1;
		unsigned long long niter=0;
//		while(x[1]!=1){
		while(x[1]>=min1st){
			++niter;
			if(x[h]==2){
				++m; x[h--]=1;
			}else{
				r=x[h]-1; t=m-h+1; x[h]=r;
				while(t>=r){
					x[++h]=r; t-=r;
				}
				if(t==0) m=h;
				else{
					m=h+1;
					if(t>1) x[++h]=t;
				}
			}

			if(m > k) continue;
			for(int i=1; i<=m; ++i) *(out++)=x[i] ;
			for(int i=m+1; i<=k; ++i)  *(out++)=0;

			if(++nsol==nsols) break;
		}
		delete[] x;

		Rprintf("niter=%d\n", niter);
		if(nsol < nsols){
			SEXP ans = PROTECT(allocMatrix(INTSXP, k, nsol));
			memcpy(INTEGER(ans), INTEGER(outR), sizeof(int) * k * nsol);
			UNPROTECT(1);
			return(ans);
		}else {
			SEXP ans=PROTECT(allocVector(LGLSXP, 1));
			(LOGICAL(ans))[0] = 1;
			UNPROTECT(1);
			return(ans);
		}
	}
}
