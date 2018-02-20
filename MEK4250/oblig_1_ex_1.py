from dolfin import *
import numpy as np

class Right(SubDomain):
    def inside(self, x, on_boundary, eps=1.e-14):
        return x[0] > 1 - eps


class Left(SubDomain):
    def inside(self, x, on_boundary, eps=1.e-14):
        return x[0] < eps


def solve_ex_1(N=32, degree=1, k=1):
    mesh = UnitSquareMesh(N, N)
    V = FunctionSpace(mesh, 'P', degree)

    #numerical and analytical expressions
    u_e = Expression('sin(pi*k*x[0])*cos(pi*k*x[1])', degree=degree, k=k)
    u_num = Function(V)

    #Boundary conditions
    bcs = [DirichletBC(V, Constant(0.0), Left()),
           DirichletBC(V, Constant(1.0), Right())]

    u = TrialFunction(V)
    v = TestFunction(V)

    f = Constant(1.0)
    a = dot(grad(u), grad(v)) * dx
    L = f * v * dx

    u = Function(V)
    solve(a == L, u, bcs)

    #plot u and mesh
    plot(u)
    plot(mesh)

    output = File('sol_ex_1.pvd')
    output << u

def prob_1_b():
    pass

def error():
    pass

def main():
    prob_1_b()


if __name__ == '__main__':
    main()

