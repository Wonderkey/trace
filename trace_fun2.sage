def B(u,t,n,N):
    ans=euler_phi(N)/euler_phi(N//u)
    count=0
    for x in range(N):
        if gcd(x,N)==1:
            if ((x^2-t*x+n)%(N*u)==0):
                count=count+1
    ans=ans*count
    return ans

def C(t,u,N,n):
    ans=0
    for d in divisors(u):
        ans=ans+B(u//d,t,n,N)*moebius(d)
    return ans

def Phi(a,d,N):
    ans=0
    for r in divisors(N):
        s=N//r
        if gcd(N,a-d) % gcd(r,s) == 0:
            ans=ans+euler_phi(gcd(r,s))
    return ans
#print Phi(1,3,11)

#print C(0,1,11,3),C(0,2,11,3),C(1,1,11,3),C(2,1,11,3),C(3,1,11,3),C(4,1,11,3),C(4,2,11,3)

#for u in divisors(11):
 #   for t in range(4):
  #      print C(t,u,3,11)


