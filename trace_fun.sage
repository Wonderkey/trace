def orde(l,n):
    o=0
    while n % l^o ==0:
        o=o+1
    return o-1


def cc(s,f,N,n,ll):
    count=0
    A=0
    B=0
    if (s^2-4*n)/f^2 % ll ==0:
        bb=1
    else:
        bb=0
    v=orde(ll,N)
    b=orde(ll,f)
    for x in range(ll^(v+b)):
        if (2*x-s) % ll^b ==0:
            Phi=x^2-s*x+n
            if Phi % ll^(v+2*b)==0:
                A=A+1
                if (b==1 and (Phi % ll^(v+2*b+1)) == 0):
                    B=B+1
    count=A+B
    return count

def c(s,f,N,n):
    count=1
    for x in prime_range(2,sqrt(N)+1):
        if N % x ==0:
            count=count * cc(s,f,N,n,x)
    if is_prime(N):
        count=count * cc(s,f,N,n,N)
    return count

def tt(s,n):
    tmp=s^2-4*n
    div=tmp.divisors()
    div.reverse()
    print div
    t0=1
    for t in div:
        if (t<=sqrt(tmp)) or (tmp % t^2 ==0):
            t0=t
            break
    if ((tmp//t0^2) % 4 == 1):
        return t0
    else:
        return t0//2


def V(n,N):
    vec=[]
    s=0
    while (s^2-4*n<0):
        t=Integer(tt(s,n))
        try:
            for f in divisors(t):
                vec.append(c(s,f,N,n))
                vec.append(c(-s,f,N,n))
            print 1
        except ValueError:
            pass
        s=s+1
    try:
        for d in n.divisors():
            if d^2>n:
                break
            for f in (n//d-d).divisors():
                vec.append(c(n//d+d,f,N,n))
                print 2
    except ValueError:
        pass
    return vec


print c(4,2,11,3)

print c(1,1,11,3)

print V(3,11)
