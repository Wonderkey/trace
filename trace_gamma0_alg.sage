
load("trace_gamma0_new.sage")


def M_dic(n, chi, newform_case = False):
    s_set = gen_s_set(n)
    M_dict = {}
    mo = chi.modulus()
    for p in prime_divisors(mo):
        if gcd(p, n) == 1:
            M_dict[p] = 2 * mo.ord(p) + 1
    for s in s_set:
        if s >= 0:
            det = s ^ 2 - 4 * n
            for p in prime_divisors(det):
                if gcd(p, n) == 1:
                    if M_dict.has_key(p):
                        M_dict[p] = max(Integer(det).ord(p) + 1, M_dict[p])
                    else:
                        M_dict[p] = Integer(det).ord(p) + 1
    if newform_case:
        for p in M_dict.keys():
            M_dict[p] += 2
    return M_dict


def c_new(s, f, N, n):
    tmp = 0
    for d in divisors(Integer(N)):
        tmp += arith_beta(1, d) * c(s,f,N/d,n,trivial_character(N))
    return tmp


def pre_cal(n, k):
    s_set = gen_s_set(n)
    aa = {} #a_dic
    for s in s_set:
        aa[s] = hij_a(s,k,n)
    bb = {} #b_dic
    for s in s_set:
        for f in divisors(tt(s,n)):
            bb[(s,f)] = hij_b(s, f, n)
    return (aa, bb)


def find_zero_levels(val, n, ss_set, s_0, num, M, a_dic, b_dic, chi):
    #print(num)
    if num == 0: # calculation part
        res = []
        s_set = gen_s_set(n)
        til_N = chi.modulus()
        for tmp_N0 in divisors(M/til_N):
            N0 = tmp_N0 * til_N
            tmp_chi = chi.extend(N0)
            tmp_tr = 0
            for s in s_set:
                if (s not in s_0) and (-s not in s_0):
                    tmp = 0
                    for f in divisors(tt(s,n)):
                        tmp += b_dic[(s,f)] * c(s, f, N0, n, tmp_chi)
                    tmp_tr += a_dic[s] * tmp
            if type(tmp_tr) is type(1+I):
                tmp_tr = tmp_tr.simplify_full()
            #print(s_0,N0)
            #print(tmp_tr)
            if val == -1: #k > 2
                if tmp_tr == 0:
                    res.append((N0,s_0))
            else: #k == 2
                tmp_sig = sigma(n)
                for it in range(0, val - len(prime_divisors(N0)) + 1):
                    if tmp_tr * (2 ^ it) == tmp_sig:
                        res.append((N0, it, s_0))
        return res
    else: #search part
        s = ss_set[num - 1]
        res = find_zero_levels(val, n, ss_set, copy(s_0), num - 1, M, a_dic, b_dic, chi)
        s_0.append(s)
        res.extend(find_zero_levels(val, n, ss_set, copy(s_0), num - 1, M, a_dic, b_dic, chi))
        return res


def find_zero_levels_newform(val, n, ss_set, s_0, s_1, num, M, a_dic, b_dic, chi = trivial_character(1)):   #only_trivial_char_version_available
    #print(num)
    if num == 0: # calculation part
        res = []
        s_set = gen_s_set(n)
        til_N = chi.modulus()
        for tmp_N0 in divisors(M/til_N):
            N0 = tmp_N0 * til_N
            tmp_chi = chi.extend(N0)
            tmp_tr = 0
            if len(s_0) != 0:
                eps = 0
            else:
                eps = 1
            for s in s_set:
                if s in ss_set:
                    if (s not in s_0) and (-s not in s_0):
                        tmp = 0
                        for f in divisors(tt(s,n)):
                            tmp += b_dic[(s,f)] * c_new(s, f, N0, n)
                        if (s not in s_1) or (-s not in s_1):
                            tmp_tr -= a_dic[s] * tmp
                        else:
                            tmp_tr += a_dic[s] * tmp
                else:
                    if eps == 1:
                        tmp = 0
                        for f in divisors(tt(s,n)):
                            tmp += b_dic[(s,f)] * c(s, f, N0, n, tmp_chi)
                        tmp_tr += a_dic[s] * tmp
            if type(tmp_tr) is type(1+I):
                tmp_tr = tmp_tr.simplify_full()
            #print(s_0,N0)
            #print(tmp_tr)
            if val == -1: #k > 2
                if tmp_tr == 0:
                    res.append((N0,eps,s_0,s_1))
          #  else: #k == 2
          #      tmp_sig = sigma(n)
          #      for it in range(0, val - len(prime_divisors(N0)) + 1):
          #          if tmp_tr * (2 ^ it) == tmp_sig:
          #              res.append((N0, it, s_0))
        return res
    else: #search part
        s = ss_set[num - 1]
        res = find_zero_levels_newform(val, n, ss_set, copy(s_0), copy(s_1), num - 1, M, a_dic, b_dic, chi)
        res.extend(find_zero_levels_newform(val, n, ss_set, copy(s_0)+[s], copy(s_1), num - 1, M, a_dic, b_dic, chi))
        res.extend(find_zero_levels_newform(val, n, ss_set, copy(s_0), copy(s_1)+[s], num - 1, M, a_dic, b_dic, chi))
        return res


def trace_gamma0_alg(n, k, chi = trivial_character(1), newform_case = False, verbose = true): 
    M_dict = M_dic(n,chi,newform_case)
    M = 1
    for p in M_dict.keys():
        M *= p ^ (M_dict[p])
    # Step 1
    (a_dic, b_dic) = pre_cal(n, k) #pre-calculate
    s_0 = []
    ss_set = list(filter(lambda x: x>=0, gen_ss_set(n)))
    if k == 2 and chi.is_trivial(): #haven't consider newforms
        tmp = 3
        for d in divisors(n):
            if d ^ 2 > n:
                break
            tmp = max(tmp, (n - d^2).ord(2))
        val = len(prime_divisors(M)) + sigma(n).ord(2) + tmp
        res = find_zero_levels(val, n, ss_set, s_0, len(ss_set), M, a_dic, b_dic, chi) 
    else: #DFS
        if newform_case:
            s_1 = []
            res = find_zero_levels_newform(-1, n, ss_set, s_0, s_1, len(ss_set), M, a_dic, b_dic, chi)
        else:
            res = find_zero_levels(-1, n, ss_set, s_0, len(ss_set), M, a_dic, b_dic, chi) 
    # Step 2
    # end
    if verbose:
        if newform_case:
            print("zero-level N_0, epsilon and epsilon_s = 0,1 terms are:")
        else:
            print("zero-level N_0 and epsilon_s = 0 terms are:")
        for it in res:
            print(it)
    return





