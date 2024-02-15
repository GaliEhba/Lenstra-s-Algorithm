import time
import random
from math import gcd
def add_ec(P,Q,A,B,N):
    """
    Add 2 points P & Q on the Elliptic curve E over Z_N
    E : y^2 = x^3 + Ax +B
    A,B are in Z_N
    Z_N finte ring Z_N = Z/NZ
    P, Q on E(Z_N)
    """
    O=(0,1,0) # point at infinity
    x1,y1,z1 = P[0],P[1],P[2]
    x2,y2,z2 = Q[0],Q[1],Q[2]
    if z1 > 1 :
        return P
    if z2 > 1 :
        return Q
    if z1 == 0: #point at infinity
        return Q
    elif z2 == 0: #point at infinity
        return P
    elif (z1 == 0 and z2 == 0): #point at infinity
        return O
    elif x1 == x2:
        if y1 == y2: # doubling
            up = (3*x1**2+A)%N    
            down = (2*y1)%N
            try:
                m = (up*pow(down,-1,N))%N     #slope of tangent to E at the point P
                x3 = (m**2 - x1 - x2)%N
                y3 = (m*(x1-x3)-y1)%N
                return(x3,y3,1)
            except:
                #print("Error1 while trying to find inverse of {} modulo {} !".format(down,N))
                return (0,0,down)
        else:
            return O       # Q = -P, P+Q = O
    else:        #adding two distinct points
        up = (y2 - y1)%N
        down = (x2 - x1)%N
        try:
            m = (up*pow(down,-1,N))%N # slope of chord PQ
            x3 = (m**2-x1-x2)%N
            y3 = (m*(x1-x3)-y1)%N
            return(x3,y3,1)
        except:
            #print("Error2 while trying to find inverse of {} modulo {} !".format(down,N))
            return (0,0,down)



def multiply_ec(n,P,A,B,N):
    """
    Double and add algorithm on elliptic curve
    Perform the multiplication [n]P = P+P+...+P  n times on the elliptic curve E over the ring Z_N
    E : y^2 = x^3 + Ax +B
    P : an element of E(Z_N)
    n : a positive integer
    """
    O = (0,1,0) # point at infinity
    Q = O  # Initialize Q
    if P[2] == 0:  #inf
        return O
    if P[2]>1:     #error on inversion
        return P
    while n>0:
        if P[2] == 0:  #inf
            return O
        if P[2]>1:     #error on inversion
            return P
        if n%2 == 1:
            Q = add_ec(P,Q,A,B,N)
        n = n//2
        P = add_ec(P,P,A,B,N)
    return Q



def erastosthens_sieve(n):
    list1 = [True for i in range(n+1)]
    list2 = []
    for p in range(2,n+1):
        if list1[p]:
            list2.append(p)
            for i in range(p,n+1,p): # erase all multiples of p
                list1[i] = False
    return list2 # list of all primes less than n



def lenstras_factor(N,bound):
    """
    function try to factorize an integer N using elliptic curve of the form 
    y^2 = x^3 + Ax + B over the field Z_N
    bound: maximun value of iteration before giving up usualy max e_iq_i , #E(Z_N)=prod (q_i)^e_i
    """
    g = N # randomize the choice of curve and an initial point P
    
    while g == N:  # make sure we will choose good parameters
        P = (random.randint(0,N-1), random.randint(0,N-1),1) # choose a random point P
        A = random.randint(0,N-1) # choose random coefficient A
        B = (P[1]**2 - P[0]**3 - A*P[0])%N # compute B to get P on E : y^2 = x^3 + Ax + B
        
        g = gcd(4*A**3 + 27*B**2,N) #singularity condition
    if g > 1:
        return g # we have got a trivial factor of N
    else:
        logfile = erastosthens_sieve(bound)
        for p in logfile: # for p prime less than bound
            p1 = p
            while p1 < bound :
                P = multiply_ec(p,P,A,B,N) # compute Q = [p]*P
                if P[2]>1:   # error while computing inverse
                    return gcd(P[2],N)  # factor of N
                p1 = p*p1  # until we reach lcm(1,...bound)
    return("Failed to find divisor please try again")


ti = time.time()
print(lenstras_factor(24962571002895724912726972956970761011,1000000))
tf = time.time()
print("runtime: ",tf-ti)
