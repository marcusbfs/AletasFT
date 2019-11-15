from sympy import *
from sympy.printing.cxxcode import cxxcode

# definições
a, b, D, L = symbols('a b ft_D ft_L', Real=True)
z = symbols('z', Real=True)

# coloque F(z) aqui
# fz = a + b * z
fz = a -  b*cosh(z)
crescente = True
crescente = False


if (crescente):
    print("b")
    bb = solve(fz.subs(z,L) - D, b)[0]
    print(cxxcode(bb))
else:
    print("b")
    bb = solve(fz.subs(z,L) - 0, b)[0],
    print(cxxcode(bb))

print("")
print("a")
aa = solve(fz.subs(z,0).subs(b, bb) - D/2, a)[0]
print(cxxcode(aa))



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

