from sympy import *
from sympy.printing.cxxcode import cxxcode

# definições
a, b, D, L = symbols('a b ft_D ft_L', Real=True)
z = symbols('z', Real=True)

# coloque F(z) aqui
fz = a - b*exp(z)
crescente = True
crescente = False

if (crescente):
    eq1 = -D + fz.subs(z,L)
else:
    eq1 =  fz.subs(z,L)

eq2 = -D/2 + fz.subs(z,0)


# Solve equations for a and b
sol = nonlinsolve([eq1, eq2], [a, b])

print('All solutions')
print(sol)
print("")

sol, = nonlinsolve([eq1, eq2], [a, b])
print('a:')
print(cxxcode(sol[0]))
print("")
print('b:')
print(cxxcode(sol[1]))
print("")

# print("")
# print("a")
# aa = solve(fz.subs(z,0).subs(b, bb) - D/2, a)[0]
# print(cxxcode(aa))

print("")



# Função
print("F(z)")
# print(fz)
print(cxxcode(fz))
print("")

# Derivada
dFdz = diff(fz, z)
print("dF(z)/dz")
# print(dFdz)
print(cxxcode(dFdz))
print("")

# As
# As = 2*pi*integrate(fz*sqrt(1+dFdz**2), (z, 0, z))
# print("As(z)")
# print(As)
# print("")


# dAsDz
# dAsdz = diff(As,z)
# print("dAs(z)/dz")
# # print(dAsdz)
# print(cxxcode(dAsdz))
# print("")

