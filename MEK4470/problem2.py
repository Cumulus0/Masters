'''
example 5.1
Upwind/Hybrid/Power-law/QUICK
chapter: 5.6/5.7/5.8/5.9

include L^2 feil
'''
from numpy import *
from pylab import *

# Step 1, mesh
def Upwind(N):
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


def QUICK(N):
    faces = linspace(0, 1., N+1)
    nodes = 0.5*(faces[1:] + faces[:-1])
    delta = faces[1]-faces[0]

    # constants
    L = 1.0
    phi_a = 1.0
    phi_b = 0.0
    F = 0.2 # rho*u
    D = 0.1 / (1.0/N) # gamma/dx

    # Set up linear system of equations
    A = zeros((N, N))
    b = zeros(N)

    for i in range(2, N-1):
        aw = D + (6/8)*F + (1/8)*F
        aww = -(1/8)*F
        ae = D - (3/8)*F
        aee = 0
        ap = aw + ae + aww + aee
        sp = 0
        su = 0
        A[i, i-1] = -aw
        A[i, i-2] = -aww
        A[i, i]   = ap
        A[i, i+1] = -ae
        b[i] = su

    # Node 1
    ae = D + (1/3)*D - (3/8)*F
    aw = 0
    aww = 0
    sp = -((8/3)*D + (2/8)*F + F)
    su = ((8/3)*D + (2/8)*F + F) * phi_a
    ap = aww + aw + ae - sp
    A[0, 0] = ap
    A[0, 1] = -ae
    b[0] = su

    # Node 2
    ae = D - (3/8)*F
    aw = D + (7/8)*F + (1/8)*F
    aww = 0
    sp = 0.25*F
    su = -0.25*F*phi_a
    ap = aww + aw + ae - sp
    A[1, 1] = ap
    A[1, 0] = -aw
    A[1, 2] = -ae
    b[1] = su
    

    # Node N
    aw = D + (1/3)*D + (6/8)*F
    aww = -(1/8)*F
    ae = 0
    sp = -((8/3)*D - F)
    su = ((8/3)*D - F)*phi_b
    ap = aww + aw + ae - sp
    A[-1, -1] = ap
    A[-1, -2] = -aw
    A[-1, -3] = -aww
    b[-1] = su
    return nodes,A,b

def analytical(x_):
    F=0.2
    return (-1 * (exp(F * (x_/0.1)) - 1) / (exp(F * (1.0/0.1)) - 1)) + 1

def error(num, anal, N):
    #return linalg.norm(num - anal)
    return sqrt(1/N * sum((num - anal)**2))

if __name__ == '__main__':
    err = []
    N = list(range(5,21))

    for i in N:
        nodes,A,b = QUICK(i)
        T = solve(A, b)
        exact = analytical(nodes)
        err.append(error(T,exact, i))
    print(A)
    figure(1)
    subplot(211)
    plot(nodes, T, 'b')
    plot(nodes, exact, 'r')
    xlabel('X')
    ylabel('Phi')
    title('Phi distribution')

    subplot(212)
    loglog(1 / asarray(N), err)
    loglog(1 / asarray(N), err, 1 / asarray(N), err[0]*(1 / asarray(N) * N[0]))
    xlabel('Grid Points')
    ylabel('L2 Error')
    title('Error vs. # of nodes')
    show()