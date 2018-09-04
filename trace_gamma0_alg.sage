
load("trace_gamma0.sage")


def M_dic(n):
    s_set = gen_s_set(n)
    M_dict = {}
    for s in s_set:
        if s >= 0:
            det = s ^ 2 - 4 * n
            for p in prime_divisors(det):
                if gcd(p, n) == 1:
                    if M_dict.has_key(p):
                        M_dict[p] = max(Integer(det).ord(p) + 1, M_dict[p])
                    else:
                        M_dict[p] = Integer(det).ord(p) + 1
    return M_dict


def pre_cal(n, k, M_dic):
    s_set = gen_s_set(n)
    aa = {} #a_dic
    for s in s_set:
        aa[s] = hij_a(s,k,n)
    bb = {} #b_dic
    c_dic = {} #c_dic
    for s in s_set:
        for f in divisors(tt(s,n)):
            bb[(s,f)] = hij_b(s, f, n)
            for p in M_dic.keys():
                p = Integer(p)
                for it in range(1, M_dic[p] + 1):
                    c_dic[(s,f,p,it)] = cc(s,f,p^it,n,p,trivial_character(p^it))
    return (aa, bb, c_dic)


def find_zero_levels(n, ss_set, s_0, num, M, a_dic, b_dic, c_dic):
    print(num)
    if num == 0: # calculation part
        res = []
        s_set = gen_s_set(n)
        for N0 in divisors(M):
            tmp_tr = 0
            for s in s_set:
                if s not in s_0:
                    tmp = 0
                    for f in divisors(tt(s,n)):
                        ccc = 1
                        for p in prime_divisors(N0):
                            ccc *= c_dic[(s, f, p, N0.ord(p))] 
                        tmp += b_dic[(s,f)] * ccc
                    tmp_tr += a_dic[s] * tmp
            if type(tmp_tr) is type(1+I):
                tmp_tr = tmp_tr.simplify_full()
            print(s_0,N0)
            print(tmp_tr)
            if tmp_tr == 0:
                res.append((N0,s_0))
        return res
    else: #search part
        s = ss_set[num - 1]
        res = find_zero_levels(n, ss_set, copy(s_0), num - 1, M, a_dic, b_dic, c_dic)
        s_0.append(s)
        res.extend(find_zero_levels(n, ss_set, copy(s_0), num - 1, M, a_dic, b_dic, c_dic))
        return res


def trace_gamma0_alg(n,k): # k>2
    M_dict = M_dic(n)
    M = 1
    for p in M_dict.keys():
        M *= p ^ (M_dict[p])
    # Step 1
    (a_dic, b_dic, c_dic) = pre_cal(n, k, M_dict) #pre-calculate
    s_0 = []
    ss_set = gen_ss_set(n)
    res = find_zero_levels(n, ss_set, s_0, len(ss_set), M, a_dic, b_dic, c_dic) #DFS
    # Step 2
    # end
    print("zero-level N_0 and epsilon_s = 0 terms are:")
    for it in res:
        print(it)
    return












