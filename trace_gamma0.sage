def order_of_int(l,n): #return m such that l^m|n
    if n == 0:
        raise ValueError('n cannot be 0')
    o=1
    while n % (l^o) ==0:
        o=o+1
    return o-1


def roots_of_unity_num(d):
    if d == -4:
        return 4
    elif d == -3:
        return 6
    else:
        return 2


def h0(n):
    d = n.squarefree_part()
    if mod(d,4) != 1:
        d = d * 4
    m = sqrt(n/d)
    #part 1
    if roots_of_unity_num(n) == 2:
        h = -m / abs(d)
    else:
        h = -m / abs(d) * roots_of_unity_num(d) / 2
    #part 2
    tmp = 0
    for it in range(1,abs(d)):
        tmp = tmp + kronecker(d,it)*it
    h = h * tmp
    #part 3
    for it in prime_divisors(m):
        h = h * (1-kronecker(d,it)/it)
    return h * 2 / roots_of_unity_num(n)


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
    v = order_of_int(ll, N)
    b = order_of_int(ll, f)
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
    tmp = s^2-4*n
    div = divisors(tmp)
    div.reverse()
    t0=1
    for t in div:
        if (t <= sqrt(tmp)) or (tmp % (t^2) == 0):
            t0 = t
            break
    if ((tmp/t0^2) % 4 == 1):
        return t0
    else:
        return t0/2

'''
def V(n,N):
    vec=[]
    s=0
    while (s^2-4*n<0):
        t=Integer(tt(s,n))
        for f in divisors(t):
            vec.append(c(s,f,N,n))
            vec.append(c(-s,f,N,n))
        s=s+1
    for d in n.divisors():
        if d^2>n:
            break
        for f in (n//d-d).divisors():
            vec.append(c(n//d+d,f,N,n))
    return vec
'''

def par(ll, N, cond):
    nu = order_of_int(ll, N)
    rho = (nu/2).floor()
    ee = order_of_int(ll, cond)
    if ee >= rho + 1:
        return 2 * (ll^(nu - ee))
    elif nu % 2 == 0:
        return ll^rho + ll^(rho - 1)
    else:
        return 2 * ll^rho


def phi_1(N):
    res = N
    for p in prime_divisors(N):
        res = res * (1 + 1 / p)
    return res


def trace_gamma0(n,k,N,chi):
    #part1:abc
    tr = 0
    s_set = [] # s^2-4n square or negative
    s_set.append(0)
    #s_set.append(0)
    it = 1
    while it^2 - 4 * n < 0:
        s_set.append(it)
        s_set.append(-it)
        it += 1
    for it in divisors(n):
        if it^2 >= n:
            break
        s_set.append(n / it + it)
        s_set.append(- n / it - it)
    for s in s_set:
        tmp = 0
        for f in divisors(tt(s,n)):
            tmp += hij_b(s,f,n) * c(s,f,N,n,chi)
        tr -= hij_a(s,k,n) * tmp
    #part2:n_is_a_square
    if is_square(n):
        n=floor(n)
        #print(parent(sqrt(n)))
        tr = tr + (phi_1(N) / 12 * (k - 1) * n^(k/2 - 1) * chi(sqrt(n)))
        tmp2 = 1
        cond = chi.conductor()
        for p in prime_divisors(N):
            tmp2 = tmp2 * par(p, N, cond)
        tr = tr - (n^(k/2 - 1) * chi(sqrt(n)) * sqrt(n) / 2 * tmp2)
    #part3:k==2_and_trivialchar
    if k == 2 and chi.is_trivial():
        for d in divisors(n):
            tr += n/d
    #end
    tr = tr.simplify_full()
    #print('the trace of T({}) acting weight {} Gamma_0({}) is {}'.format(n, k, N, tr))
    return tr


def real_trace_gamma0(n,k,N,chi):
    S = CuspForms(Gamma0(N),k)
    tr = S.hecke_matrix(n).trace()
    #print('the REAL trace of T({}) acting weight {} Gamma_0({}) is {}'.format(n, k, N, tr))
    return tr


def test1(n,k,N):
    for kk in range(2,k+1,2):
        for nn in range(1,n+1):
            for NN in range(1,N+1):
                if gcd(nn,NN) == 1:
                    ans1 = trace_gamma0(nn,kk,NN,trivial_character(NN))
                    ans2 = real_trace_gamma0(nn,kk,NN,trivial_character(NN))
                    if(ans1 == ans2):
                        print("yes")
                    else:
                        print("no")


