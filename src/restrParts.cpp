#include "config.h"

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>

#define MAX0(x) ((x)>0 ? (x) : 0)
#define MIN(x,y)  ((x)>(y) ? (y) : (x))
#define CEILQ(x,y) ((x)/(y) + ((x) % (y) != 0))

extern "C" {
	SEXP restrParts(SEXP xR, SEXP ctR, SEXP minctR, SEXP maxctR, SEXP ctallowR, SEXP valuesR, SEXP nextvalR, SEXP diffvalsR, SEXP outR, SEXP nsolsR)
	{
		int * x; x=INTEGER(xR);
		int * ct; ct=INTEGER(ctR);
		int * minct; minct=INTEGER(minctR);
		int * maxct; maxct=INTEGER(maxctR);
		int * ctallow; ctallow=INTEGER(ctallowR);
		int * values; values=INTEGER(valuesR);
		int * nextval; nextval=INTEGER(nextvalR);
		int * diffvals; diffvals=INTEGER(diffvalsR);
		int * out; out=INTEGER(outR);
		unsigned int nsols; nsols=(unsigned int)INTEGER(nsolsR)[0];
		unsigned int nvals; nvals=(unsigned int)(LENGTH(valuesR)-1);
		
		unsigned int nsol=0;
		unsigned int niter=0;
		unsigned int lev = 1;
		bool nextLev=true;
		do{
			if (nextLev){
				while(values[lev] > x[lev-1] && lev < nvals){
					// fastforward
					ct[lev]=maxct[lev]=minct[lev]=0;
					ctallow[lev]=ctallow[lev-1];
					x[lev]=x[lev-1];
					++lev;
				}
				maxct[lev] = MIN(ctallow[lev-1], (int)(x[lev-1]/values[lev]));
				minct[lev] = MAX0(CEILQ( x[lev-1L]-ctallow[lev-1L]*nextval[lev], diffvals[lev] )); 
				nextLev = false;
				
				ct[lev] = maxct[lev] + 1;
			}
			if(ct[lev] <= minct[lev]) {
				--lev; 
				goto nextiter;
			} else --ct[lev];
			
			++niter;
			x[lev] = x[lev-1] - ct[lev] * values[lev];
			ctallow[lev] = ctallow[lev-1] - ct[lev];
/*			
			Rprintf("\nx =\n");
			for(unsigned int j=1; j<=nvals; ++j) Rprintf("%d ", x[j]);

			Rprintf("\ncounts =\n");
			for(unsigned int j=1; j<=nvals; ++j) Rprintf("%d ", ct[j]);

			Rprintf("\nmin =\n");
			for(unsigned int j=1; j<=nvals; ++j) Rprintf("%d ", minct[j]);
			Rprintf("\nmax =\n");
			for(unsigned int j=1; j<=nvals; ++j) Rprintf("%d ", maxct[j]);
			Rprintf("\n rem =\n");
			for(unsigned int j=1; j<=nvals; ++j) Rprintf("%d ", ctallow[j]);
			
			Rprintf("\nlev=%d\n", lev);
*/
			if (x[lev] == 0){
				// found a partition
				++nsol;
				for(unsigned int j=1; j<=lev; ++j)
					for(int i=0; i<ct[j]; ++i)
						*(out++)=values[j];
				for(int j=0; j<ctallow[lev]; ++j) *(out++)=0;
				if(nsol == nsols) break;
//				Rprintf("nsol=%d\n", nsol);
			}else if(lev < nvals){
				++lev;
				nextLev = true;
			}

nextiter:			
			while(ct[lev] == minct[lev] && !nextLev && lev>0)
				--lev; // fastrewind

			R_CheckUserInterrupt();
		}while(lev>0);
	
//		Rprintf("niter=%d\n", niter);
		if(nsol < nsols){
			SEXP ans = PROTECT(allocMatrix(INTSXP, ctallow[0], nsol));
			memcpy(INTEGER(ans), INTEGER(outR), sizeof(int) * ctallow[0] * nsol);
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
