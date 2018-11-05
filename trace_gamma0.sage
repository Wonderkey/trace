# Some functions in number theory

def roots_of_unity_num(d):
    if d == -4:
        return 4
    elif d == -3:
        return 6
    else:
        return 2


def h0(n):
    return (Integer(n).class_number()) * 2 / roots_of_unity_num(n)


def phi_1(N):
    res = N
    for p in prime_divisors(N):
        res = res * (1 + 1 / p)
	return res


def arith_beta(m, n): # n is a positive integer
    if n == 1:
        return 1
    count = 1
    for p in prime_divisors(n):
        d = n.ord(p)
        if m % p == 0:
            if d >= 2:
                return 0
            else:
                count *= (-1)
        else:
            if d >= 3:
                return 0
            elif d == 1:
                count *= (-2)
            else:
                pass
    return count


def squarefull_part(n):
    if n == 1:
        return 1
    count = 1
    for p in prime_divisors(n):
        if Integer(n).ord(p) > 1:
            pass
        else:
            count *= p
    return n / count


#Some functions in the trace formulas

def hij_a(s, k, n):
    #xx = var('xx')
    #yy = xx ^2 - s * xx + n
    #sol = solve(yy == 0, xx)
    #x1 = sol[0].rhs()
    #x2 = sol[1].rhs()
    det = s ^ 2 - 4 * n
    x1 = (s + sqrt(det)) / 2
    x2 = (s - sqrt(det)) / 2
    if det < 0:
        return (x1 ^ (k - 1) - x2 ^ (k - 1)) / (2 * (x1 - x2))
    else:
        return sgn(x1) ^ k * min(abs(x1), abs(x2)) ^ (k-1) / abs(x1 - x2)


def hij_b(s, f, n):
    det = s ^ 2 - 4 * n
    if det < 0:
        return h0(det / f ^ 2)
    else:
        return euler_phi(sqrt(det) / f) / 2


def cc(s, f, N, n, ll, chi):
    A=[]
    B=[]
    count = 0
    v = Integer(N).ord(ll)
    b = Integer(f).ord(ll)
    for x in range(ll ^ (v + b)):
        if (2 * x - s) % (ll ^ b) ==0:
            if (x ^ 2 - s * x + n) % (ll ^ (v + 2 * b)) == 0:
                A.append(x)
    if ((s^2 - 4 * n) / f ^ 2) % ll == 0:
        for x in A:
            if (x ^ 2 - s * x + n) % (ll ^ (v + 2 * b + 1)) == 0:
                B.append(s - x)
    for x in A:
        count += chi(x)
    for x in B:
        count += chi(x)
    return count


def c(s, f, N, n, chi): #require N != 1
    count = 1
    for chip in chi.decomposition():
        p = prime_divisors(chip.modulus())[0]
        count = count * cc(s, f, N, n, p, chip)
    return count


def tt(s,n):
    det = s^2-4*n
    div = divisors(det)
    div.reverse()
    t0=1
    for t in div:
        if det % (t^2) == 0:
            t0 = t
            break
    if ((det / t0^2) % 4 == 1):
        return t0
    else:
        return t0/2


def par(ll, N, cond):
    nu = Integer(N).ord(ll)
    rho = (nu/2).floor()
    ee = Integer(cond).ord(ll)
    if ee >= rho + 1:
        return 2 * (ll^(nu - ee))
    elif nu % 2 == 0:
        return ll^rho + ll^(rho - 1)
    else:
        return 2 * ll^rho


def gen_es_set(n, neg=True): #generate "s" in the elliptic terms
    s_set = [0]
    it = 1
    if neg == True:
        while it^2 - 4 * n < 0:
            s_set.append(it)
            s_set.append(-it)
            it += 1
    else:
        while it^2 - 4 * n < 0:
            s_set.append(it)
            it += 1
    return s_set


def gen_us_set(n, neg=True):
    s_set = []
    if neg == True:
        for it in divisors(n):
            if it^2 >= n:
                break
            s_set.append(n / it + it)
            s_set.append(- n / it - it)
    else:
        for it in divisors(n):
            if it^2 >= n:
                break
            s_set.append(n / it + it)
    return s_set


# Trace formula of Hecke operators

def trace_gamma0(n, k, N, chi=False, verbose=False):
    if chi == False:
        chi = trivial_character(N)
    # part1: abc
    tr = 0
    s_set = gen_es_set(n) + gen_us_set(n)
    for s in s_set:
        tmp = 0
        for f in divisors(tt(s,n)):
            tmp += hij_b(s,f,n) * c(s,f,N,n,chi)
        tr -= hij_a(s,k,n) * tmp
    # part2: n_is_a_square
    if is_square(n):
        n=floor(n)
        #print(parent(sqrt(n)))
        tr = tr + (phi_1(N) / 12 * (k - 1) * n^(k/2 - 1) * chi(sqrt(n)))
        tmp2 = 1
        cond = chi.conductor()
        for p in prime_divisors(N):
            tmp2 = tmp2 * par(p, N, cond)
        tr = tr - (n^(k/2 - 1) * chi(sqrt(n)) * sqrt(n) / 2 * tmp2)
    # part3: k==2 and trivialchar (assume (n,N)==1)
    if k == 2 and chi.is_trivial():
        for d in divisors(n):
            tr += n/d
    # end
    tr = tr.simplify_full()
    if verbose:
        print('the trace of T({}) acting on S_{}({}, chi) is {}'.format(n, k, N, tr))
    return tr


def real_trace_gamma0(n, k, N, chi=False, verbose=False):
    if chi == False:
        chi = trivial_character(N)
    S = ModularSymbols(Gamma0(N),k,+1).cuspidal_subspace()
    tr = S.hecke_matrix(n).trace()
    if verbose:
        print('the REAL trace of T({}) acting on S_{}({}, chi) is {}'.format(n, k, N, tr))
    return tr


def trace_gamma0_new(n, k, N, chi=False, verbose=False):
    if chi == False:
        chi = trivial_character(N)
    tr = 0
    N1 = N / squarefull_part(N)
    cond = chi.conductor()
    chip = chi.restrict(cond)
    for dd in divisors(N / cond):
        d = cond * dd
        for t in divisors(gcd(dd, N1)):
            if n % (t^2) == 0:
                bet = arith_beta(n/(t^2), N/d)
                #print(bet)
                if bet == 0:
                    pass
                else:
                    tr += chip(t) * (t^(k-1)) * bet * trace_gamma0(n/t^2, k, d/t, chi.restrict(d/t), verbose)
                    #print(trace_gamma0(n/(t^2),k,d/t,chip))
    tr = tr.simplify_full()
    if verbose:
        print('the trace of T({}) acting on S_{}^new({}, chi) is {}'.format(n,k,N,tr))
    return tr


def real_trace_gamma0_new(n, k, N, chi=False, verbose=False):
    if chi == False:
        chi = trivial_character(N)
    S = ModularSymbols(Gamma0(N),k,+1).cuspidal_subspace()
    SS = S.new_subspace()
    tr = SS.hecke_matrix(n).trace()
    if verbose:
        print('the REAL trace of T({}) acting on S_{}^new({}, chi) is {}'.format(n,k,N,tr))
    return tr


# Test functions

def test1(n,k,N):
    for kk in range(2,k+1,2):
        for nn in range(1,n+1):
            for NN in range(1,N+1):
                if gcd(nn,NN) == 1:
                    ans1 = trace_gamma0(nn,kk,NN,trivial_character(NN))
                    ans2 = real_trace_gamma0(nn,kk,NN,trivial_character(NN))
                    #print(ans1,ans2)
                    if(ans1 == ans2):
                        print("yes")
                    else:
                        print("no")
                        

def test1_new(n,k,N):
    for kk in range(2,k+1,2):
        for nn in range(3,n+1):
            for NN in range(3,N+1):
                if gcd(nn,NN) == 1:
                    ans1 = trace_gamma0_new(nn,kk,NN,trivial_character(NN))
                    ans2 = real_trace_gamma0_new(nn,kk,NN,trivial_character(NN))
                    if(ans1 != ans2):
                        print(kk, nn, NN)
