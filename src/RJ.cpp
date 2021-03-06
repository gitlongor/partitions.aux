#include "config.h"

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>


#include "integer_arith.h"
#define MIN(x,y)  ((x)>(y) ? (y) : (x))


extern "C" {

#ifndef NOTHING
// slightly rewritten version w/o goto statements and optimized for the case of z=0: 
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
	// in this case z=0; y[i]=l1; j=0 (except for loop idx);
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
			while (m > l2) {
				m -= l2; 
				x[i++] = l1 + l2; // keep alloc upper bnd;
			}
			x[i] = l1 + m; ++num;  // last part;
			for(int tmp=1; tmp<=n; ++tmp) *(out++) = x[tmp];
			if (num == nsols) break; // goto end;

			if (i<n && m>1) { // possible to increase # of parts?
							  // if so, do it;
							  // but after this, no way to extend it;
							  // have to backtrack.
				m = 1; --x[i++]; x[i] = l1 + 1;
				++num;
				for(int tmp=1; tmp<=n; ++tmp) *(out++) = x[tmp];
				if (num == nsols) break; //goto end;
			}
			for (j=i-1; j>=1; --j) { // backtracking
				l2 = x[j] - l1 - 1; ++m;  // go back to location j and consider the next smaller value;
				if (m <= (n-j) * l2) { // setting all >j locations to current save value as reduced x[j], can we get a solution?
					x[j] = l1 + l2;
					break; // goto a1;
				}else{  // otherwise, keep going back
					m += l2; x[i] = l1; i = j; // i always points to j's next location.  At the end: i=1.
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
#else

// original implementation for comparison	
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
	const int z=0; // z=*INTEGER(zR);
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
#endif


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

SEXP RJa0ss(SEXP mR, SEXP l1R, SEXP l2R, SEXP nR, SEXP ssR, SEXP numR, SEXP outR)
/*
	z=0: 0 for arbitrary choice of parts; 1 for distinct parts;
	m: number to be partitioned;
	l1: lower bound of each part;
	l2: upper bound of each part;
	n: number of parts;
	ss: upper bound on sum of squares;
	num: number of partitions to find;
	out: preallocated space to hold results.
*/
{
	//const int z=0; // z=*INTEGER(zR);
	// in this case z=0; y[i]=l1; j=0 (except for loop idx);
	int m; m=*INTEGER(mR);
	int l1; l1=*INTEGER(l1R);
	int l2; l2=*INTEGER(l2R);
	int n; n=*INTEGER(nR);
	int ss; ss=*INTEGER(ssR);
	int nsols; nsols=*INTEGER(numR);
	int *out; out=INTEGER(outR);
	int num;

	int *x, i, j; 
	const int l1sq = l1 * l1, l12= 2 * l1;
	int toAdd;
	x=new int[n+1];
	//y=new int[n+1];

	num = 0; m -= n * l1 ; l2 -= l1;  // diff bounds
	ss -= n * l1sq; 
	if (m>=0 && m <= n*l2 && ss>=0){
		for(i=1; i<=n; ++i) x[i] = l1; // init lower bnd;
		
		i=1; // l2 -= z*(n-1);
		do{
a1:
			// keep allocating current allowed upper bound
			while (m > (toAdd = MIN(l2, isqrt(l1sq+ss)-l1)) && toAdd > 0 && i<n) 
			{// toAdd holds allowed upper bound
				m -= toAdd; ss -= toAdd * ( toAdd + l12); 
				x[i++] = l1 + toAdd; 
			}
			if (toAdd != 0 && i<=n && m<=toAdd) {
				x[i] = l1 + m; ++num;  // last part;
				for(int tmp=1; tmp<=n; ++tmp) *(out++) = x[tmp];
				if (num == nsols) break; // goto end;
			
				if (i<n && m>1) { // possible to increase # of parts?
								  // if so, do it;
								  // but after this, no way to extend it w/o affecting decreasing order -- have to backtrack then.
					// note that this never decrease ss: 
					// original sum of squares at location i and i+1 is 
					//		(l1+toAdd)^2 + l1^2; 
					// new sum of squares at these 2 locations is 
					//		(l1+toAdd-1)^2 + (l1+1)^2 = (original ss) -2*(toAdd-1).
					ss -= (x[i]-1)*(x[i]-1); // ss=ss+x[i]^2-(x[i]-1)^2
					m = 1; --x[i++]; x[i] = l1 + 1; 
					++num;
					for(int tmp=1; tmp<=n; ++tmp) *(out++) = x[tmp];
					if (num == nsols) break; //goto end;
				}
			}
			for (j=i-1; j>=1; --j) { // backtracking
				// tentatively reduce x[j] 1-by-1 and check if there could be solutions: 
				const int n_j = n-j; 
				const int m0 = m; 
				const int zUpper = x[j]-l1 - idiv_ceil(x[j]-l1+m0, n_j+1) ;
				for(int z=1; z <= zUpper; ++z){
					//further check ss: 
					//if all x[i]..x[n] are assigned to values to minimize their sum of squares, will that exceed allowed ss? 
					++m;
					const int k = m / n_j; const int rem = m % n_j; 
					if(m*(l12+k) + (k+1)*rem <= ss+2*z*x[j] - z*z) {
						// both necessary conditions (m and ss) are met
						ss += 2*z*x[j] - z*z; // ss + x[j]^2 - (x[j]-z)^2
						x[j] -= z; 
						l2 = x[j] - l1; 
						goto a1;
					}
				}
				// decreasing x[j] is not possible, keep going back
					x[i--]=l1; // i always points to j's next location.  At the end: i=1.
					ss += x[i]*x[i] - l1sq; 
					m = m0 + x[i] - l1; 
			}
		}while(i>1);
	}
//end:
	delete[] x;

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


SEXP RJb(SEXP muR, SEXP mR, SEXP vR, SEXP nR, SEXP numR, SEXP outR)
/*
	mu: array of upper bounds for each allowed value v;
	m: number to be partitioned;
	v: array of allowed values for each part;
	n: number of parts;
	num: number of partitions to find;
	out: preallocated space to hold results.
*/
{
	int * mu;
	int m; m=*INTEGER(mR);
	int * v;
	int n; n=*INTEGER(nR);
	int nsols; nsols=*INTEGER(numR);
	int *out; out=INTEGER(outR);
	int num;
	int r; r=LENGTH(vR);

	mu=new int[n+1]; memcpy(mu+1, INTEGER(muR), sizeof(int) * r);
	v =new int[n+1]; memcpy(v+1,  INTEGER(vR),  sizeof(int) * r);

	int *x, *y, *ii, i, j, k, ll, lr;
	x=new int[n+1];
	y=new int[n+1];
	ii=new int[n+1];

	j=1; k=mu[1]; ll=v[1]; num=0;
	for(i=n; i>0; --i){
		x[i] = y[i] = ll; --k; m -= ll;
		if(k == 0){
			if(j == r) goto b3;
			k = mu[++j]; ll = v[j];
		}
	}
	lr = v[r]; ll = v[1];
	if (m < 0 || m > n*(lr-ll)) goto b3;
	if (m ==0){
		num = 1;
		for(int tmp=1; tmp<=n; ++tmp) *(out++) = x[tmp];
		if (num == nsols) goto end;
		goto b3;
	}

	i = 1; m += y[1];
b1:
	for (j = mu[r]; j>0; --j){
		if (m <= lr) goto b2;
		x[i] = lr; ii[i] = r - 1;
		m -= lr - y[++i];
	}
	--r;
b2:
	for (j=r; v[r] > m; --r);
	lr = v[r];
	if (m == lr){
		x[i] = lr; ++num;
		for(int tmp=1; tmp<=n; ++tmp) *(out++) = x[tmp];
		if (num == nsols) goto end;
		--r; lr = v[r];
	}
	k = y[i];
	if (lr > k && m - lr <= (n - i) * (lr - ll)) goto b1;
	else x[i] = k;
	for (--i; i>0; --i){
		r = ii[i]; lr = v[r]; m += x[i] - k; k = y[i];
		if (lr > k && m - lr <= (n - i) * (lr - ll)) goto b1;
		else x[i] = k;
	}
b3:
end:
	delete[] mu;
	delete[] v;
	delete[] x;
	delete[] y;
	delete[] ii;

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
