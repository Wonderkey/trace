
load("trace_gamma0_alg.sage")


def gen_generating_fun(n, N, eps_dic):
    # N, n coprime, n not a square
    # get the generating function given n, N and epsilon_s's
    chi = trivial_character(N)
    R = sigma(n)
    s_set = gen_es_set(n)
    d_set = list(filter(lambda x: x < sqrt(n), divisors(n)))
    for s in s_set:
        tmp = 0
        for f in divisors(tt(s, n)):
            tmp += hij_b(s, f, n) * c(s, f, N, n, chi)
        R -= eps_dic[abs(s)] * tmp * (n * x + 1) / (n^2 * x^2 + (2*n-s^2) * x + 1) / 2
    for d in d_set:
        tmp = 0
        s = n / d + d
        for f in divisors(tt(s, n)):
            tmp += hij_b(s, f, n) * c(s, f, N, n, chi)
        R += tmp * (2 * d^2) / (n - d^2) / (d^2 * x - 1) 
    R = R.simplify_full()
    return R


def list_generating_fun_search(n, N, eps_dic, s_num, current_s, output):
    # the DFS part of list_generating_fun
    if current_s == s_num:
        print("N0 = {}, eps = ".format(N)),
        for k in range(s_num):
            print(k, eps_dic[k]),
        print("")
        print("generating function:"),
        gen_fun = gen_generating_fun(n, N, eps_dic)
        print(gen_fun)
        print("")
        output[tuple([N] + eps_dic.values())] = gen_fun
    else:
        eps_dic[current_s] = 0
        list_generating_fun_search(n, N, eps_dic.copy(), s_num, current_s + 1, output)
        eps_dic[current_s] = 1
        list_generating_fun_search(n, N, eps_dic.copy(), s_num, current_s + 1, output)
    return 0


def list_generating_fun(n, N):
    # get all the generating function given n and N
    output = {}
    s_num = len(gen_es_set(n, neg=False))
    for N0 in divisors(N):
        list_generating_fun_search(n, N0, {}, s_num, 0, output)
    return output









