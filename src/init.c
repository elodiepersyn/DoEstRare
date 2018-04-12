#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void aSum_stat(void *, void *, void *, void *, void *, void *, void *, void *);
extern void CAST_stat(void *, void *, void *, void *, void *, void *, void *);
extern void DoEstRare_stat(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void Sum_stat(void *, void *, void *, void *, void *, void *, void *);
extern void wSum_stat(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP rMFNCHypergeo_c(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
  {"aSum_stat",      (DL_FUNC) &aSum_stat,       8},
  {"CAST_stat",      (DL_FUNC) &CAST_stat,       7},
  {"DoEstRare_stat", (DL_FUNC) &DoEstRare_stat, 11},
  {"Sum_stat",       (DL_FUNC) &Sum_stat,        7},
  {"wSum_stat",      (DL_FUNC) &wSum_stat,      13},
  {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
  {"rMFNCHypergeo_c", (DL_FUNC) &rMFNCHypergeo_c, 5},
  {NULL, NULL, 0}
};

void R_init_DoEstRare(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}