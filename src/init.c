#include "config.h"
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


#ifdef HAVE_CXX1X
	extern SEXP fxt_partition_gen(SEXP, SEXP, SEXP, SEXP);
//	extern SEXP anyDupAtomMatHash(SEXP, SEXP, SEXP);
//	extern SEXP grpDupAtomMatHash(SEXP, SEXP, SEXP);
//	extern int initHash(void);


static R_CallMethodDef callMethods[]  = {
  {"fxt_partition_gen", (DL_FUNC) fxt_partition_gen, 4},
//  {"anyDupAtomMat", (DL_FUNC) &HASHED_NAME(anyDupAtomMat), 3},
//  {"grpDupAtomMat", (DL_FUNC) &HASHED_NAME(grpDupAtomMat), 3},
//  {"dbl_dig", (DL_FUNC) &dbl_dig, 0},
  {NULL, NULL, 0}
};
	

void R_init_partitions_aux(DllInfo *info)
{
   R_registerRoutines(info, NULL, callMethods, NULL, NULL);
   R_useDynamicSymbols(info, FALSE);
}
#endif
