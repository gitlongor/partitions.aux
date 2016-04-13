#include "config.h"

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>

class ykkn {
	public:
		int * * x; // current status;
		int k;   // number of allowed parts;
		int n;	 // number to be partitioned;
		int * m; 	 // number of parts in the current partition;
		//int * out;

		int * output (int row, int * out, bool pad0=false) const
		{
			int * thisRow = x[row];
			int i;
			for(i=0; i<=m[row]; ++i) {
				*(out++)=thisRow[i];
//				Rprintf("%d ", thisRow[i]);
			}
//			Rprintf("\n");
			if (pad0) for(; i<k; ++i) *(out++)=0;
			return (out);
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
	out = output (row, out, false);
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
	out = output(0, out, false);
	out = findAllChildren1(1, out);
	return out;
}

extern "C" {

SEXP ykknAllParts1(SEXP nR, SEXP outR)
{
	int * out = INTEGER(outR);
	int n = *INTEGER(nR);
	ykkn ykknTree(n);
	int * out0 = ykknTree.allParts1(out);
	SEXP ans=PROTECT(allocVector(LGLSXP, 1));
	LOGICAL(ans)[0] = out0 > out ? 1 : 0;
	UNPROTECT(1);
	return ans;
}

}
