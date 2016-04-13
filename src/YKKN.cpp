#include "config.h"

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>

#ifdef NOTHING
// this is the initial O(n^2) storage version of the algorithm;
class ykkn {
	public:
		int * * x; // current status;
		int k;   // number of allowed parts;
		int n;	 // number to be partitioned;
		int * m; 	 // number of parts in the current partition;
		//int * out;

		int * output (int row, int * out, bool pad0=false) const
		{
			int mp1 = m[row]+1;
			memcpy(out, x[row], sizeof(int)*mp1);
			if (pad0) memset(out + mp1, 0, sizeof(int) * (k-mp1));
			return out + k;
		}
		int * findAllChildren1(int row, int * out);
		int * allParts1(int * out);

		ykkn (int nvalue)
		{
			ykkn_(nvalue, nvalue);
		}
		ykkn (int nvalue, int kvalue)
		{
			ykkn_(nvalue, kvalue);
		}
		void ykkn_ (int nvalue, int kvalue)
		{
			n = nvalue;
			k = kvalue;
			x = new int*[k];
			for(int i = 0; i < k; ++i){
				x[i] = new int[k];
			}
			x[0][0]=n;
			m = new int[k];
			m[0]=0;
		}
		~ykkn ()
		{
			delete[] m;
			for(int i=0; i < k; ++i) delete[] x[i];
			delete[] x;
		}
	private:
		void AmPlus1(int row)
		{
			memcpy(x[row+1], x[row], sizeof(int) * (m[row] + 1));
			m[row+1] = m[row]+1;
			--x[++row][0];
			x[row][m[row]]=1;
		}
		void Am(int row)
		{
			memcpy(x[row+1], x[row], sizeof(int) * (m[row] + 1));
			m[row+1]=m[row];
			--x[++row][0];
			++x[row][m[row]];
		}
};

int * ykkn::findAllChildren1(int row, int * out)
{
	R_CheckUserInterrupt();
	out = output (row, out);
	if (x[row][0] > x[row][1]){
		AmPlus1(row);
//Rprintf("m[%d]=%d\n", row, m[row]);
//Rprintf("m[%d]=%d\n", row+1, m[row+1]);
//for(int i=0; i<=m[row+1]; ++i)
//	Rprintf("x[%d][%d]=%d\n", row+1, i, x[row+1][i]);
		out = findAllChildren1(row+1, out);
		if (x[row][m[row]-1] > x[row][m[row]] &&
			(m[row]>1 || x[row][m[row]-1] - x[row][m[row]] > 1)) {
				Am(row);
				out = findAllChildren1(row+1, out);
		}
	}
	return out;
}

int * ykkn::allParts1(int * out)
{
	x[0][0]=n;
	m[0]=0;
	x[1][0]=n-1;
	x[1][1]=1;
	m[1]=1;
	out = output(0, out);
	out = findAllChildren1(1, out);
	return out;
}

#else

// this is the O(n) storage version; 
class ykkn {
	public:
		int * x; // current status;
		int k;   // number of allowed parts;
		int n;	 // number to be partitioned;
		int m; 	 // number of parts in the current partition;
		//int * out;

		inline void output (bool pad0=false) 
		{
			int mp1 = m+1;
			memcpy(out_, x, sizeof(int)*mp1);
			if (pad0) memset(out_ + mp1, 0, sizeof(int) * (k-mp1));
			out_ += k;
		}
		void findAllChildren1();
		void findAllChildren1Stack();
		void findAllChildren2(int depth);
		int * allParts1(int * out);
		int * allParts2(int * out);

		ykkn (int nvalue)
		{
			ykkn_(nvalue, nvalue);
		}
		ykkn (int nvalue, int kvalue)
		{
			ykkn_(nvalue, kvalue);
		}
		void ykkn_ (int nvalue, int kvalue)
		{
			n = nvalue;
			k = kvalue;
			x = new int[k];
			
			x[0]=n;
			m=0;
		}
		~ykkn ()
		{
			delete[] x;
		}
	private:
		int * out_;
		#define AMPLUS1 --x[0]; x[++m]=1;
		#define AM --x[0]; ++x[m];
		#define INVAMPLUS1 ++x[0]; x[m--]=0;
		#define INVAM ++x[0]; --x[m];
};

void ykkn::findAllChildren1()
{// algorithm 1 based on recursion
	R_CheckUserInterrupt();
	output();
	if (x[0] > x[1]){
		AMPLUS1
		findAllChildren1();
		INVAMPLUS1
		if (x[m-1] > x[m] &&
			(m>1 || x[m-1] - x[m] > 1)) {
				AM
				findAllChildren1();
				INVAM
		}
	}
}

void ykkn::findAllChildren1Stack()
{// algorithm 1 based on manual stack management
	bool * stack = new bool[k];
	unsigned int depth = 0;
start:	
	R_CheckUserInterrupt();
	output();
	if (x[0] > x[1]){
		AMPLUS1
		stack[depth++] = true ;
		goto start; //findAllChildren1();
backTrue:
		INVAMPLUS1
		if (x[m-1] > x[m] &&
			(m>1 || x[m-1] - x[m] > 1)) {
				AM
				stack[depth++] = false;
				goto start; //findAllChildren1();
backFalse:
				INVAM
		}
	}
	if (depth==0) {
		delete[] stack;
		return;
	}
	if (stack[--depth]) goto backTrue; else goto backFalse;
}

void ykkn::findAllChildren2(int depth)
{
	R_CheckUserInterrupt();
	bool even = (depth % 2 == 0);
	if (even) output();
	if (x[0] > x[1]){
		AMPLUS1
		findAllChildren2(depth+1);
		INVAMPLUS1
		if (x[m-1] > x[m] &&
			(m>1 || x[m-1] - x[m] > 1)) {
				AM
				findAllChildren2(depth + 1);
				INVAM
		}
	}
	if (!even) output();
}

int * ykkn::allParts1(int * out)
{
	out_ = out;
	x[0]=n;
	m=0;
	output (); 
	x[0]=n-1;
	x[1]=1;
	m=1;
	//findAllChildren1();
	findAllChildren1Stack();
	return out_;
}
int * ykkn::allParts2(int * out)
{
	out_ = out;
	x[0]=n;
	m=0;
	output (); 
	x[0]=n-1;
	x[1]=1;
	m=1;
	findAllChildren2(1);
	return out_;
}

	
#endif

extern "C" {

SEXP ykknAllParts(SEXP nR, SEXP outR, SEXP methodR)
{
	int * out = INTEGER(outR);
	int n = *INTEGER(nR);
	ykkn ykknTree(n);
	int * out0; 
	out0 =  (*INTEGER(methodR) == 1) ? ykknTree.allParts1(out) : ykknTree.allParts2(out) ;
	
	SEXP ans=PROTECT(allocVector(LGLSXP, 1));
	LOGICAL(ans)[0] = out0 > out ? 1 : 0;
	UNPROTECT(1);
	return ans;
}

}
