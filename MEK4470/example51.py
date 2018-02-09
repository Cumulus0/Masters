import numpy as np 
from scipy.sparse import diags
from matplotlib import pyplot as plt

L = 1 # m

phi_A = 1 #transport property
phi_B = 0
np.
F = 0.1 # rho*u
D = 0.5 #Gamme/delta_x
a_e_factor = D - F/2
a_w_factor = D + F/2

def a_matrix(N):
    a_W =  a_w_factor * np.ones(N)
    a_E = a_e_factor * np.ones(N)
    a_W[0] = 0
    a_E[-1] = 0
    s_p = np.zeros(N)
    s_p[0] = -(2*D + F)
    s_p[-1] = -(2*D - F)
    a_P = a_W + a_E - s_p
    s_u = np.zeros(N)
    s_u[0] = (2*D + F) * phi_A
    s_u[-1] = (2*D - F) * phi_B
    A = diags([-a_E[:-1], -a_W[1:], a_P], [1,-1,0])
    return((A.toarray(), s_u))

def analytical():
    pass

def error():
    return np.linalg.norm()

def L2_order():
    pass

def main():
    A = a_matrix(5)
    return A


if __name__ == '__main__':
   A, s_u = main()
   s_u.reshape((1,5))
   phi = np.linalg.solve(A.toarray(), s_u)
   plt.plot(phi)
   plt.show()
#plot i log-log scale