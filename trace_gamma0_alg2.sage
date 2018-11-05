
load("trace_gamma0_alg.sage")


def gen_generating_fun(n, N, eps_dic):# N, n coprime, n not a square
    chi = trivial_character(N)
    R = sigma(n)
    s_set = list(filter(lambda x: x >= 0, gen_ss_set(n)))
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







