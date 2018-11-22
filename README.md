# the traces of Hecke operators

**trace_gamma0.sage**
use trace formulas to calculate the traces of Hecke operators
- trace_gamma0(n, k, N, chi=False, verbose=False)
- trace_gamma0_new(n, k, N, chi=False, verbose=False)


**trace_gamma0_alg.sage**
find zero levels $N$ of the traces of Hecke operator $T_n$ given weight $k$ and congruence subgroups $\Gamma(N, chi)$
- trace_gamma0_alg(n, k, chi = trivial_character(1), newform_case = False, verbose = true)

**trace_gamma0_alg2.sage**
to be continue
- list_generating_fun(n, N, verbose=True, zero_filter=False)
- auto_check_fun_coeff_mod(mod_list, rational_fun, end_num, begin_num=False, speed=5, verbose=True)
- find_period_of_coeff_mod(mod_list, rational_fun, verbose=True)
- find_zeros_in_coeff_mod(mod_list, rational_fun, verbose=True)

