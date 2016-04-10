#include "config.h"
#ifdef HAVE_CXX1X

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>

#include <comb/partition-gen.h>

extern "C" {

SEXP fxt_partition_gen(SEXP n, SEXP k, SEXP values, SEXP nout)
{
	SEXP out;
	int i;
	int intk, intn, intnout, nvals;
	ulong ulk;
	int const *  valp;
	valp = INTEGER(values);

	intk = *INTEGER(k); ulk=(ulong)intk;
	intn = *INTEGER(n); //uln=(ulong)intn;
	intnout =*INTEGER(nout);
	nvals = LENGTH(values);

	out = PROTECT(allocMatrix(INTSXP, intk, intnout));
	int * outp;
	outp = INTEGER(out);
	for(i=0; i<intk*intnout; ++i) *(outp++)=R_NaInt;

	ulong *pv;
	pv = new ulong[nvals];
	for(i=0; i<nvals; ++i) pv[i]=(ulong) valp[i];

	partition_gen pp((ulong) intn, (ulong) nvals, pv);

	int pci; ulong checkCounts, cti, rowi; //, tmp;
	for ( pp.next(); outp>INTEGER(out); pp.next() ){
		for(checkCounts=0, pci=0; pci<nvals; ++pci) checkCounts += pp.pc_[pci];
		R_CheckUserInterrupt();
		
		if(checkCounts > ulk) continue;
/*
		Rprintf("\nCounts = %d\n%d * %d", checkCounts, pp.pc_[0], pp.pv_[0]);
		for(int tmp=1; tmp<nvals; ++tmp)
			Rprintf("+ %d * %d", pp.pc_[tmp], pp.pv_[tmp]);
		Rprintf("\n");
*/
		for(pci=nvals-1; pci>=0; --pci)
			for(cti=pp.pc_[pci]; cti>0; --cti){
				*(--outp)=valp[pci];
//				Rprintf("%d\t", valp[pci]);
			}
		for(rowi=checkCounts; rowi < ulk; ++rowi) *(--outp)=0;
//		Rprintf("\nDone.\n");
		
	}

	delete[] pv;
	UNPROTECT(1);
	return out;
}




}

#endif
