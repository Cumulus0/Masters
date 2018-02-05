'''
example 5.1
Upwind/Hybrid/Power-law/QUICK
chapter: 5.6/5.7/5.8/5.9

include L^2 feil
'''
from numpy import *
from pylab import *

# Step 1, mesh
def main(N):
    faces = linspace(0, 1., N+1)
    nodes = 0.5*(faces[1:] + faces[:-1])
    delta = faces[1]-faces[0]

    # constants
    L = 1.0
    phi_a = 1.0
    phi_b = 0.0
    F = 0.1 # rho*u
    D = 0.1 / (1.0/N) # gamma/dx

    # Set up linear system of equations
    A = zeros((N, N))
    b = zeros(N)

    # ap*T_P - aw*T_W - ae*T_E = 0
    for i in range(1, N-1):
        aw = D + F
        ae = D
        sp = 0
        su = 0
        ap = aw + ae - sp
        su = 0
        A[i, i-1] = -aw
        A[i, i]   = ap
        A[i, i+1] = -ae
        b[i] = su

    # Node 1
    ae = D
    aw = 0
    sp = -(2*D + F)
    su = (2*D + F) * phi_a
    ap = aw + ae - sp
    A[0, 0] = ap
    A[0, 1] = -ae
    b[0] = su

    # Node N
    aw = D + F
    ae = 0
    sp = -2*D
    su = 2*D*phi_b
    ap = aw + ae - sp
    A[-1, -1] = ap
    A[-1, -2] = -aw
    b[-1] = su
    return nodes,A,b

def analytical(x_):
    F=0.1
    return (-1 * (exp(F * (x_/0.1)) - 1) / (exp(F * (1.0/0.1)) - 1)) + 1

def error(num, anal, N):
    #return linalg.norm(num - anal)
    return sqrt(1/N * sum((num - anal)**2))

if __name__ == '__main__':
    err = []
    N = list(range(5,21))

    for i in N:
        nodes,A,b = main(i)
        T = solve(A, b)
        exact = analytical(nodes)
        err.append(error(T,exact, i))


    figure(1)
    #subplot(211)
    #plot(nodes, T, 'b', label='numerical')
    #plot(nodes, exact, 'r', label='exact')
    #subplot(212)
    #loglog(1 / asarray(N), err)
    loglog(1 / asarray(N), err, 1 / asarray(N), err[0]*(1 / asarray(N) * N[0]))
    show()