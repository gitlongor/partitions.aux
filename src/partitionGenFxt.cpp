#include "config.h"
#ifdef HAVE_CXX1X

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>

#include <comb/partition-gen.h>

extern "C" {

/*
 fxt_partition_gen is based on the partition-en in fxt library.
 It finds all partitions with each part in "values" (excluding 0).
 Partitions with more than k parts are then ignored.
*/

SEXP fxt_partition_gen(SEXP n, SEXP k, SEXP values, SEXP nout)
{
	SEXP ans, tmpans;
	int i;
	int intk, intn, intnout, nvals;
	ulong ulk;
	int const *  valp;
	valp = INTEGER(values);

	intk = *INTEGER(k); ulk=(ulong)intk;
	intn = *INTEGER(n); //uln=(ulong)intn;
	intnout =*INTEGER(nout);
	nvals = LENGTH(values);

	tmpans = PROTECT(allocMatrix(INTSXP, intk, intnout));
	int * outp;
	outp = INTEGER(tmpans);

	ulong *pv;
	pv = new ulong[nvals];
	for(i=0; i<nvals; ++i) pv[i]=(ulong) valp[i];

	partition_gen pp((ulong) intn, (ulong) nvals, pv);

	int pci; ulong checkCounts, cti, rowi,  nsols, unout=intnout;

//	for ( pp.next(); outp>INTEGER(tmpans); pp.next() ){
	for ( nsols=0; nsols < unout;  ){
		pp.next();
		for(checkCounts=0, pci=0; pci<nvals; ++pci) checkCounts += pp.pc_[pci];
		R_CheckUserInterrupt();

		if(checkCounts > ulk ) continue;
		if (checkCounts == 0) break;
		++nsols;
/*
		Rprintf("\nCounts = %d\n%d * %d", checkCounts, pp.pc_[0], pp.pv_[0]);
		for(int tmp=1; tmp<nvals; ++tmp)
			Rprintf("+ %d * %d", pp.pc_[tmp], pp.pv_[tmp]);
		Rprintf("\n");
*/
		for(pci=nvals-1; pci>=0; --pci)
			for(cti=pp.pc_[pci]; cti>0; --cti){
				*(outp++)=valp[pci];
//				Rprintf("%d\t", valp[pci]);
			}
		for(rowi=checkCounts; rowi < ulk; ++rowi) *(outp++)=0;
//		Rprintf("\nDone.\n");

	}
	delete[] pv;

	ans=PROTECT(allocMatrix(INTSXP, intk, (int)nsols));
	memcpy(INTEGER(ans), INTEGER(tmpans), sizeof(int) * intk * nsols);
	
	UNPROTECT(2);
	return ans;
}




}

#endif
