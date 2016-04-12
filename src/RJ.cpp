#include "config.h"

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>

extern "C" {

SEXP RJa0(SEXP mR, SEXP l1R, SEXP l2R, SEXP nR, SEXP numR, SEXP outR)
/*
	z=0: 0 for arbitrary choice of parts; 1 for distinct parts;
	m: number to be partitioned;
	l1: lower bound of each part; 
	l2: upper bound of each part;
	n: number of parts;
	num: number of partitions to find;
	out: preallocated space to hold results. 
*/
{
	//const int z=0; // z=*INTEGER(zR); 
	// in this case z=0; y[i]=0; j=0 (except for loop idx); 
	int m; m=*INTEGER(mR); 
	int l1; l1=*INTEGER(l1R);
	int l2; l2=*INTEGER(l2R); 
	int n; n=*INTEGER(nR); 
	int nsols; nsols=*INTEGER(numR);
	int *out; out=INTEGER(outR); 
	int num;
	
	int *x, i, j;
	x=new int[n+1]; 
	//y=new int[n+1]; 
	
	num = 0; m -= n * l1 ; l2 -= l1;  // diff bounds
	if (m>=0 && m <= n*l2 ){
		for(i=1; i<=n; ++i) x[i] = l1; // init lower bnd;
		i=1; // l2 -= z*(n-1);
		do{
//a1:			
			while (m > l2) m -= (x[i++] = l2); // keep alloc upper bnd; 
			x[i] = m; ++num;  // last part; 
			for(int tmp=1; tmp<=n; ++tmp) *(out++) = x[tmp]; 
			if (num == nsols) break; // goto end;
			
			if (i<n && m>1) { // possible to increase # of parts? 
							  // if so, do it; 
							  // but after this, no way to extend it; 
							  // have to backtrack. 
				m = 1; --x[i++]; x[i] = 1; 
				++num; 
				for(int tmp=1; tmp<=n; ++tmp) *(out++) = x[tmp]; 
				if (num == nsols) break; //goto end;
			} 
			for (j=i-1; j>=1; --j) { // backtracking
				l2 = x[j] - 1; ++m;  // go back to location j and consider the next smaller value; 
				if (m <= (n-j) * l2) { // setting all >j locations to current save value as reduced x[j], can we get a solution? 
					x[j] = l2;  
					break; // goto a1;
				}else{  // otherwise, keep going back
					m += l2; x[i] = 0; i = j; // i always points to j's next location.  At the end: i=1.
				}
			}
		}while(i>1);
	}
//end: 
	delete[] x;
//	delete[] y;
	
	if(num < nsols){
		*INTEGER(numR) = num;
		SEXP ans = PROTECT(allocMatrix(INTSXP, n, num));
		memcpy(INTEGER(ans), INTEGER(outR), sizeof(int) * n * num);
		UNPROTECT(1);
		return(ans);
	}else {
		SEXP ans=PROTECT(allocVector(LGLSXP, 1));
		(LOGICAL(ans))[0] = 1;
		UNPROTECT(1);
		return(ans);
	}
}



SEXP RJa1(SEXP mR, SEXP l1R, SEXP l2R, SEXP nR, SEXP numR, SEXP outR)
/*
	z=1: 0 for arbitrary choice of parts; 1 for distinct parts;
	m: number to be partitioned;
	l1: lower bound of each part; 
	l2: upper bound of each part;
	n: number of parts;
	num: number of partitions to find;
	out: preallocated space to hold results. 
*/
{
	const int z=1; // z=*INTEGER(zR); 
	int m; m=*INTEGER(mR); 
	int l1; l1=*INTEGER(l1R);
	int l2; l2=*INTEGER(l2R); 
	int n; n=*INTEGER(nR); 
	int nsols; nsols=*INTEGER(numR);
	int *out; out=INTEGER(outR); 
	int num;
	
	int *x, *y, i, j;
	x=new int[n+1]; 
	y=new int[n+1]; 
	
	num = 0; j = z*n*(n-1); m -= n * l1 - j / 2; l2 -= l1;
	if (m>=0 && m <= n*l2 -j){
		for(i=1; i<=n; ++i)
			x[i] = y[i] = l1 + z*(n-i);
		i=1; l2 -= z*(n-1);
a1:
		if (m > l2){
			m -= l2; 
			x[i] = y[i] + l2;
			++i;
			goto a1;
		}
		x[i] = y[i] + m; ++num; 
		for(int tmp=1; tmp<=n; ++tmp) *(out++) = x[tmp]; 
		if (num == nsols) goto end;
		if (i<n && m>1) {
			m = 1; --x[i++]; x[i] = y[i] + 1; 
			++num; 
			for(int tmp=1; tmp<=n; ++tmp) *(out++) = x[tmp]; 
			if (num == nsols) goto end;
		}else 
		for (j=i-1; j>=1; --j) {
			l2 = x[j] - y[j] -1; ++m;
			if (m <= (n-j) * l2) {
				x[j] = y[j] + l2;
				goto a1;
			}
			m += l2; x[i] = y[i]; i = j; 
		}
	}
end: 
	delete[] x;
	delete[] y;
	
	if(num < nsols){
		*INTEGER(numR) = num;
		SEXP ans = PROTECT(allocMatrix(INTSXP, n, num));
		memcpy(INTEGER(ans), INTEGER(outR), sizeof(int) * n * num);
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