#include "config.h"
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


#ifdef HAVE_CXX1X
	extern SEXP fxt_partition_gen(SEXP, SEXP, SEXP, SEXP);
	extern SEXP restrParts(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
	extern SEXP ZS1(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
	extern SEXP RJa0(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
	extern SEXP RJa1(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
	extern SEXP RJb(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

//	extern int initHash(void);


static R_CallMethodDef callMethods[]  = {
  {"fxt_partition_gen", (DL_FUNC) fxt_partition_gen, 4},
  {"restrParts", (DL_FUNC) restrParts, 10},
  {"ZS1", (DL_FUNC) ZS1, 6},
  {"RJa0", (DL_FUNC) RJa0, 6},
  {"RJa1", (DL_FUNC) RJa1, 6},
  {"RJb", (DL_FUNC) RJb, 6},
//  {"dbl_dig", (DL_FUNC) &dbl_dig, 0},
  {NULL, NULL, 0}
};
	

void R_init_partitions_aux(DllInfo *info)
{
   R_registerRoutines(info, NULL, callMethods, NULL, NULL);
   R_useDynamicSymbols(info, FALSE);
}
#endif
