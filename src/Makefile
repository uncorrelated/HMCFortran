DEBUG=-fbacktrace -fbounds-check -ffpe-trap=invalid,zero,overflow,underflow,inexact,denormal
FC=gfortran $(DEBUG)
CC=gcc -O3 -mtune=native -march=native -mfpmath=both
R_INC=/usr/share/R/include
SO=so
SHSRC=/var/tmp
dSFMT_if_path=$(SHSRC)/dSFMT_F03_interface/dSFMT_F03_interface/
dSFMT_path=$(SHSRC)/dSFMT_F03_interface/dSFMT-src-2.2.3/
dSFMT_c_obj=dSFMT.o dSFMT_str_state.o
dSFMT_if_c_obj=dSFMT_utils.o
dSFMT_c_src=$(addprefix $(dSFMT_path), $(patsubst %.o, %.c, $(dSFMT_c_obj))) $(addprefix $(dSFMT_if_path), $(patsubst %.o, %.c, $(dSFMT_if_c_obj)))
dSFMT_if_obj=dSFMT_interface.o
dSFMT_if_src=$(addprefix $(dSFMT_if_path), $(patsubst %.o, %.F90, $(dSFMT_if_obj)))
dSFMT=$(dSFMT_c_obj) $(dSFMT_if_obj) $(dSFMT_if_c_obj)
DDSFMT_MEXP=-DDSFMT_MEXP=19937
lbfgsb_path=$(SHSRC)/lbfgsb/src/
lbfgsb_files=lbfgsb_blas_module.F90 lbfgsb_kinds_module.F90 lbfgsb_linpack_module.f90 lbfgsb.f90
lbfgsb_src=$(addprefix $(lbfgsb_path), $(lbfgsb_files))
lbfgsb_obj=$(patsubst %.f90, %.o, $(patsubst %.F90, %.o, $(lbfgsb_files)))
COMMON=hmc_fit.o mvnorm.o imputation.o $(dSFMT) $(lbfgsb_obj)
o_names=objf_bayesian_lm objf_bayesian_logit objf_bayesian_ologit objf_bayesian_mlogit objf_poisson_exp objf_r_usr_f test_ologit_objf test_mlogit_objf test_poisson_objf r_if_bayesian_lm r_if_bayesian_logit r_if_bayesian_ologit r_if_bayesian_mlogit r_if_bayesian_poisson_exp r_if_r_usr_f objf_bayesian_imp_lm r_if_bayesian_imp_lm objf_bayesian_imp_logit r_if_bayesian_imp_logit objf_bayesian_imp_ologit r_if_bayesian_imp_ologit objf_imp_poisson_exp r_if_imp_poisson_exp
so_names=hmcfortran test_mlogit_objf test_ologit_objf test_poisson_objf

ALL: $(addsuffix .$(SO), $(so_names))

hmcfortran.$(SO): r_if_bayesian_lm.o objf_bayesian_lm.o r_if_bayesian_imp_lm.o objf_bayesian_imp_lm.o \
	r_if_bayesian_logit.o objf_bayesian_logit.o r_if_bayesian_ologit.o r_if_bayesian_mlogit.o objf_bayesian_ologit.o objf_bayesian_mlogit.o \
	r_if_bayesian_imp_logit.o objf_bayesian_imp_logit.o \
	r_if_bayesian_imp_ologit.o objf_bayesian_imp_ologit.o \
	r_if_bayesian_poisson_exp.o objf_poisson_exp.o r_if_imp_poisson_exp.o objf_imp_poisson_exp.o \
	r_if_r_usr_f.o objf_r_usr_f.o call_r_func.o $(COMMON)
	R CMD SHLIB $^ -o $@

call_r_func.o: ./call_r_func.c
	$(CC) -c $^ -o $@ -I$(R_INC)

test_mlogit_objf.$(SO): test_mlogit_objf.o objf_bayesian_mlogit.o $(COMMON)
	R CMD SHLIB $^ -o $@

test_ologit_objf.$(SO): test_ologit_objf.o objf_bayesian_ologit.o $(COMMON)
	R CMD SHLIB $^ -o $@

test_poisson_objf.$(SO): test_poisson_objf.o objf_poisson_exp.o $(COMMON)
	R CMD SHLIB $^ -o $@

imputation.o: ./imputation.f90 
	$(FC) -c $< -o $@ -fPIC

define AT_O
	$(FC) -c $(addsuffix .f90, $(addprefix ./, $(1))) -o $(addsuffix .o, $(1)) -fPIC

endef

$(addsuffix .o, $(o_names)): $(addsuffix .f90, $(addprefix ./, $(o_names))) $(COMMON)
	$(foreach fn, $(o_names), $(call AT_O, $(fn)))

hmc_fit.o: ./hmc_fit.f90 $(dSFMT_if_obj) $(lbfgsb_obj) mvnorm.o
	$(FC) -c $< -o $@

mvnorm.o: ./mvnorm.f90
	$(FC) -c $< -o $@

$(dSFMT_if_obj): $(dSFMT_if_src)
	$(FC) -c $^ $(DDSFMT_MEXP)

define AT_C
	$(CC) -c $(1) -O3 -fPIC $(DDSFMT_MEXP) -I$(dSFMT_path)

endef

$(dSFMT_c_obj) $(dSFMT_if_c_obj): $(dSFMT_c_src)
	$(foreach fn, $(dSFMT_c_src), $(call AT_C, $(fn)))

define AT_F
	$(FC) -c $(1)

endef

$(lbfgsb_obj): $(lbfgsb_src)
	$(foreach fn, $(lbfgsb_src), $(call AT_F, $(fn)))

clean:
	rm -f *.so *.o ./mod/*.mod penmp # ./src/*.o ./src/*.so
