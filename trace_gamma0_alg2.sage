
load("trace_gamma0_alg.sage")


def pre_cal_c(N, n, chi):
    pre_cal_dic = {}
    es_set = gen_es_set(n)
    for s in es_set:
        zero_piv = 1
        for f in divisors(tt(s, n)):
            pre_cal_dic[(s,f)] = c(s, f, N, n, chi)
            if pre_cal_dic[(s,f)] != 0:
                zero_piv = 0
        if zero_piv == 1:
            pre_cal_dic[s] = 0
        else:
            pre_cal_dic[s] = 1
    us_set = gen_us_set(n)
    for s in us_set:
        for f in divisors(tt(s, n)):
            pre_cal_dic[(s,f)] = c(s, f, N, n, chi)
    return pre_cal_dic


def gen_generating_fun(n, N, eps_dic, precal_c=False):
    # N, n coprime, n not a square
    # get the generating function given n, N and epsilon_s's
    chi = trivial_character(N)
    if precal_c == False:
        precal_c = pre_cal_c(N, n, chi)
    R = sigma(n)
    s_set = gen_es_set(n)
    d_set = list(filter(lambda x: x < sqrt(n), divisors(n)))
    for s in s_set:
        tmp = 0
        for f in divisors(tt(s, n)):
            tmp += hij_b(s, f, n) * precal_c[(s,f)]
        R -= eps_dic[abs(s)] * tmp * (n * x + 1) / (n^2 * x^2 + (2*n-s^2) * x + 1) / 2
    for d in d_set:
        tmp = 0
        s = n / d + d
        for f in divisors(tt(s, n)):
            tmp += hij_b(s, f, n) * precal_c[(s,f)]
        R += tmp * (2 * d^2) / (n - d^2) / (d^2 * x - 1) 
    R = R.simplify_full()
    return R


def list_generating_fun_search(n, N, eps_dic, s_num, current_s, output, verbose, zero_filter, precal_c):
    # the DFS part of list_generating_fun
    if current_s == s_num:
        if zero_filter:
            allzero_piv = 1
            for it in range(s_num):
                if eps_dic[it] > 0:
                    allzero_piv = 0
                    break
            if allzero_piv == 1:
                return 0
        gen_fun = gen_generating_fun(n, N, eps_dic, precal_c)
        if verbose:
            print("N0 = {}, eps = ".format(N)),
            for k in range(s_num):
                print(k, eps_dic[k]),
            print("")
            print("generating function:"),
            print(gen_fun)
            print("")
        output[tuple([N] + eps_dic.values())] = gen_fun
    else:
        if zero_filter and (precal_c[current_s] == 0):
            eps_dic[current_s] = -100
            list_generating_fun_search(n, N, eps_dic.copy(), s_num, current_s + 1, output, verbose, zero_filter, precal_c)
        else:
            eps_dic[current_s] = 0
            list_generating_fun_search(n, N, eps_dic.copy(), s_num, current_s + 1, output, verbose, zero_filter, precal_c)
            eps_dic[current_s] = 1
            list_generating_fun_search(n, N, eps_dic.copy(), s_num, current_s + 1, output, verbose, zero_filter, precal_c)
    return 0


def list_generating_fun(n, N, verbose=True, zero_filter=False):
    # get all the generating function given n and N
    output = {}
    s_num = len(gen_es_set(n, neg=False))
    for N0 in divisors(N):
        precal_c = pre_cal_c(N0, n, trivial_character(N0))
        eps_dic = {}
        list_generating_fun_search(n, N0, eps_dic, s_num, 0, output, verbose, zero_filter, precal_c)
    return output


def check_fun_coeff_mod(mod_list, rational_fun, term_num=50, verbose=True):
    # print the first "term_num" terms modulo each of "mod_list" of the expansion coefficients of "rational_fun"
    coeff_list = rational_fun.series(x, term_num).coefficients(sparse=False)
    out_dic = {}
    for it in mod_list:
        output = [( QQ(coeff) % it) for coeff in coeff_list]
        if verbose:
            print("coeffs mod {} are:".format(it))
            print(output)
        out_dic[(rational_fun, it)] = output
    return out_dic


def check_funs_coeff_mod(mod_list, rational_funs=False, n=False, N=False, term_num=50, verbose=True):
    # apply "check_fun_coeff_mod" to a given list of rational functions or the list given by (n, N)
    if (n != False) and (N != False):
        rational_funs = list_generating_fun(n, N, zero_filter=True)
    out_dic = {}
    for rat_fun in rational_funs:
        out_dic.update(check_fun_coeff_mod(mod_list, rat_fun, term_num, verbose))
    return out_dic












