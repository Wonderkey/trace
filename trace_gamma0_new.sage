load("trace_gamma0.sage")

def sigma0(n):
    return len(divisors(n))


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
r


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


def trace_gamma0_new(n, k, N, chi, verbose=False):
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
        print('the trace of T({}) acting on newforms weight {} Gamma_0({}) is {}'.format(n,k,N,tr))
    return tr


def real_trace_gamma0_new(n, k, N, chi, verbose=False):
    S = ModularSymbols(Gamma0(N),k,+1).cuspidal_subspace()
    SS = S.new_subspace()
    tr = SS.hecke_matrix(n).trace()
    if verbose:
        print('the REAL trace of T({}) acting on newforms weight {} Gamma_0({}) is {}'.format(n,k,N,tr))
    return tr


def test1_new(n,k,N):
    for kk in range(2,k+1,2):
        for nn in range(3,n+1):
            for NN in range(3,N+1):
                if gcd(nn,NN) == 1:
                    ans1 = trace_gamma0_new(nn,kk,NN,trivial_character(NN))
                    ans2 = real_trace_gamma0_new(nn,kk,NN,trivial_character(NN))
                    if(ans1 != ans2):
                        print(kk, nn, NN)
