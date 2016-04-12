#include "config.h"

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>

extern "C" {

SEXP RJa(SEXP zR, SEXP mR, SEXP l1R, SEXP l2R, SEXP nR, SEXP numR, SEXP outR)
/*
	z: 0 for arbitrary choice of parts; 1 for distinct parts;
	m: number to be partitioned;
	l1: lower bound of each part; 
	l2: upper bound of each part;
	n: number of parts;
	num: number of partitions to find;
	out: preallocated space to hold results. 
*/
{
	int z; z=*INTEGER(zR); 
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
		}
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