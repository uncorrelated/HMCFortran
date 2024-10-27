#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#define MANY_ARGS 0

/*
    **r_fn_ptr: Rの関数
    **envir_ptr: Rの環境（名前空間）
*/
double call_r_objf(int *np, double **p, void **r_fn_ptr, void **envir_ptr)
{
    SEXP        ans, s, t, r;
    unsigned	i;
    double      rv;

    /* Abstract Syntax Treeを構築．p, hp, X, y を渡す */

    /* 命令用の領域を確保; tとsには同じポインタが入る */
    t = s = PROTECT(allocList(2));

    /* LANGSXPを指定 */
    SET_TYPEOF(s, LANGSXP);

    /* 最初の領域には関数のポインターを入れる/SETCAR(t, install("r_callback"));とすれば関数名の指定も可能 */
    SETCAR(t, *r_fn_ptr);

    /* tを次の領域に移動させる; sは最初の領域のまま */
    t = CDR(t);

    /* 次の領域にはベクターpを入れる */
	PROTECT(r = allocVector(REALSXP, *np)); 
	for(i = 0; i<*np; i++)
	    REAL(r)[i] = (*p)[i];

    /* R用の変数を引数にセット */
    SETCAR(t, r);
    /* pと名前をつけておく */
    SET_TAG(t, install("p"));

    /* 環境*envir_ptrでevalを行う。代わりに.R_GlobalEnvを指定しても呼び出せるが、呼ばれた関数がローカル変数を読めなくなる。*/
    PROTECT(ans = eval(s, *envir_ptr));
    // PROTECT(ans = eval(R_fcall, *envir_ptr));

    /* 戻り値の型を確認 */
    if(!isReal(ans)){
        UNPROTECT(3);
        error("The length of return of 'objf' must be of numeric.");
    }

    /* Rの仕様で戻り値はベクトル */
    if(0==length(ans)){
        UNPROTECT(3);
        error("The length of return of 'objf' should be one.");
    }
    rv = REAL(ans)[0];

    /* PROTECTを外して、ガーベッジコレクターが領域を開放できるようにする */
    UNPROTECT(3);

    return rv;
}

/* グラディエント関数を呼ぶ */
void call_r_objfg(int *np, double **p, void **r_fn_ptr, void **envir_ptr, double **g)
{
    SEXP        ans, s, t, r;
    unsigned	i;

    t = s = PROTECT(allocList(2));
    SET_TYPEOF(s, LANGSXP);
    SETCAR(t, *r_fn_ptr);
    t = CDR(t);
	PROTECT(r = allocVector(REALSXP, *np)); 
	for(i = 0; i<*np; i++)
	    REAL(r)[i] = (*p)[i];
    SETCAR(t, r);
    SET_TAG(t, install("p"));
    PROTECT(ans = eval(s, *envir_ptr));
    if(!isReal(ans)){
        UNPROTECT(3);
        error("The length of return of 'objfg' must be of numeric.");
    }
    if(length(ans) < *np){
        UNPROTECT(3);
        error("The length of return of 'objfg' is less than the number of parameters.");
    }
    for(i = 0; i < length(ans); i++){
        (*g)[i] = REAL(ans)[i];
    }

    UNPROTECT(3);
}

double *fromRtoC_dbl_v(SEXP a){
    int n, i;
    double *b;

    n = length(a);
    b = (double *)R_alloc(n, sizeof(double));
    // Rprintf("n:%d %p\n", n, b);
    for(i = 0; i < n; i++)
        b[i] = REAL(a)[i];
    return b;
}

double fromRtoC_dbl(SEXP a){
    if(0 == length(a))
        error("fromRtoC_dbl: 0 length");
    return REAL(a)[0];
}

static int fromRtoC_int(SEXP a){
    if(0 == length(a))
        error("fromRtoC_int: 0 length");
    return INTEGER(a)[0];
}

int *fromRtoC_int_v(SEXP a){
    int n, i;
    int *b;

    n = length(a);
    b = (int *)R_alloc(n, sizeof(int));
    // Rprintf("n:%d %p\n", n, b);
    for(i = 0; i < n; i++)
        b[i] = INTEGER(a)[i];
    return b;
}

/* Fortranのプロシージャを関数として宣言 */
void hmc_fit_rusrf(int *N, double *sample, int *np, double *init, 
    void *r_objf, void *r_objf_envir, void *r_objfg, void *r_objfg_envir, double *h, 
    double *epsilons, int *BI, int *L, int *randlength, double *Sigma, int *constrain, int *seed, double *accept_r);

/* Rから呼び出されるCの関数/MCMC用 */
SEXP invoker(SEXP r_p, SEXP r_N, 
    SEXP r_objf, SEXP r_objf_envir, SEXP r_objfg, SEXP r_objfg_envir, SEXP r_h,
    SEXP r_epsilons, SEXP r_adjustEpsilonsN, SEXP r_L, SEXP r_randlength, SEXP r_Sigma, SEXP r_constrain, SEXP r_seed){
double  *p, h, *epsilons, *Sigma, accept_r, *sample;
int     np, seed, i, adjustEpsilonsN, randlength, L, constrain, N;
SEXP    r, m, ar;

    /* 引数が関数か確認 */
    if(!isFunction(r_objf))
        error("'r_objf' must be a function");
    if(!isEnvironment(r_objf_envir))
        error("'r_objf_envir' should be an environment");
    if(isNull(r_objfg)){
        r_objfg = r_objfg_envir = NULL;
    } else {
        if(!isFunction(r_objfg))
            error("'r_objfg' should be a function");
        if(!isEnvironment(r_objfg_envir))
            error("'r_objfg_envir' should be an environment");
    }
 
    np = length(r_p);
    h = fromRtoC_dbl(r_h);
    p = fromRtoC_dbl_v(r_p);
    N = fromRtoC_int(r_N);
    epsilons = fromRtoC_dbl_v(r_epsilons);
    adjustEpsilonsN = fromRtoC_int(r_adjustEpsilonsN);
    L = fromRtoC_int(r_L);
    randlength = fromRtoC_int(r_randlength);
    Sigma = fromRtoC_dbl_v(r_Sigma);
    constrain = fromRtoC_int(r_constrain);
    seed = fromRtoC_int(r_seed);
    sample = (double *)R_alloc(N*np, sizeof(double));
    // sample = calloc(N*np, sizeof(double));

    /* Rの関数と環境のポインターを引数に、Fortranのプロシージャを呼ぶ */
    hmc_fit_rusrf(&N, sample, &np, p, &r_objf, &r_objf_envir, &r_objfg, &r_objfg_envir, &h,
        epsilons, &adjustEpsilonsN, &L, &randlength, Sigma, &constrain, &seed, &accept_r);

    PROTECT(r = allocVector(VECSXP, 2));

    PROTECT(m = allocMatrix(REALSXP, N, np));
    for(i = 0; i < N*np; i++)
        REAL(m)[i] = sample[i];
    SET_VECTOR_ELT(r, 0, m);

    PROTECT(ar = allocVector(REALSXP, 1));
    REAL(ar)[0] = accept_r;
    SET_VECTOR_ELT(r, 1, ar);

    UNPROTECT(3);

    // free(sample);
    return r;
}

/* 以下の関数はRのoptimで同等になるため、本質的に不要 */

void log_ml_rusrf(int *nr, int *nc, int *np, double *p, 
    void *r_objf, void *r_objf_envir, void *r_objfg, void *r_objfg_envir, double *h, 
    double *lg_mlf, int *info);

/* Rから呼び出されるCの関数/対数周辺尤度計算用 */
SEXP invoke_lg_mlf(SEXP r_nr, SEXP r_nc, SEXP r_p, SEXP r_objf, SEXP r_objf_envir, SEXP r_objfg, SEXP r_objfg_envir, SEXP r_h){
    double  *p, h, lg_mlf;
    int     nr, nc, np, info; 
    SEXP    r, r_lg_mlf, r_info;

    /* 引数が関数か確認 */
    if(!isFunction(r_objf))
        error("'r_objf' must be a function");
    if(!isEnvironment(r_objf_envir))
        error("'r_objf_envir' should be an environment");
    if(isNull(r_objfg)){
        r_objfg = r_objfg_envir = NULL;
    } else {
        if(!isFunction(r_objfg))
            error("'r_objfg' should be a function");
        if(!isEnvironment(r_objfg_envir))
            error("'r_objfg_envir' should be an environment");
    }
 
    nr = fromRtoC_int(r_nr);
    nc = fromRtoC_int(r_nc);
    np = length(r_p);
    h = fromRtoC_dbl(r_h);
    p = fromRtoC_dbl_v(r_p);

    /* Rの関数と環境のポインターを引数に、Fortranのプロシージャを呼ぶ */
    log_ml_rusrf(&nr, &nc, &np, p, &r_objf, &r_objf_envir, &r_objfg, &r_objfg_envir, &h, &lg_mlf, &info);

    PROTECT(r = allocVector(VECSXP, 2));

    PROTECT(r_lg_mlf = allocVector(REALSXP, 1));
    REAL(r_lg_mlf)[0] = lg_mlf;
    SET_VECTOR_ELT(r, 0, r_lg_mlf);

    PROTECT(r_info = allocVector(REALSXP, 1));
    REAL(r_info)[0] = info;
    SET_VECTOR_ELT(r, 1, r_info);

    UNPROTECT(3);

    return r;
}

void optim_rusrf(int *np, double *p, void *r_objf, void *r_objf_envir, void *r_objfg, void *r_objfg_envir, int *nbd, double *u, double *l, double *h, double *f, double *hessian, int *info);

/* Rから呼び出されるCの関数/L-BFGS-B用 */
// SEXP r_objfg, SEXP r_objfg_envir
SEXP invoke_optim(SEXP r_p, SEXP r_objf, SEXP r_objf_envir, SEXP r_objfg, SEXP r_objfg_envir, SEXP r_h, SEXP r_nbd, SEXP r_u, SEXP r_l){
    double  *p, f, h, *hessian, *u, *l;
    int     np, i, j, *nbd, info;
    SEXP    r, r_op, r_f, r_hessian, r_info;

    /* 引数が関数か確認 */
    if(!isFunction(r_objf))
        error("'r_objf' must be a function");
    if(!isEnvironment(r_objf_envir))
        error("'r_objf_envir' should be an environment");
    if(isNull(r_objfg)){
        r_objfg = r_objfg_envir = NULL;
    } else {
        if(!isFunction(r_objfg))
            error("'r_objfg' should be a function");
        if(!isEnvironment(r_objfg_envir))
            error("'r_objfg_envir' should be an environment");
    }

    np = length(r_p);
    if(np != length(r_nbd) || np != length(r_u) || np != length(r_l))
        error("'p', 'nbd', 'u' and 'l' should have the same number of elements.");
    h = fromRtoC_dbl(r_h);
    p = fromRtoC_dbl_v(r_p);
    nbd = fromRtoC_int_v(r_nbd);
    u = fromRtoC_dbl_v(r_u);
    l = fromRtoC_dbl_v(r_l);
    hessian = (double *)R_alloc(np * np, sizeof(double));

    /* Rの関数と環境のポインターを引数に、Fortranのプロシージャを呼ぶ */
    optim_rusrf(&np, p, &r_objf, &r_objf_envir, &r_objfg, &r_objfg_envir, nbd, u, l, &h, &f, hessian, &info);
  
    PROTECT(r = allocVector(VECSXP, 4));

    PROTECT(r_op = allocVector(REALSXP, np));
    for(i = 0; i < np; i++)
        REAL(r_op)[i] = p[i];
    SET_VECTOR_ELT(r, 0, r_op);

    PROTECT(r_f = allocVector(REALSXP, 1));
    REAL(r_f)[0] = f;
    SET_VECTOR_ELT(r, 1, r_f);

    PROTECT(r_hessian = allocMatrix(REALSXP, np, np));
    for(j = 0; j < np; j++)
        for(i = 0; i < np; i++)
            REAL(r_hessian)[i + j*np] = hessian[i + j*np];
    SET_VECTOR_ELT(r, 2, r_hessian);

    PROTECT(r_info = allocVector(REALSXP, 1));
    REAL(r_info)[0] = info;
    SET_VECTOR_ELT(r, 3, r_info);

    UNPROTECT(5);

    return r;
}
