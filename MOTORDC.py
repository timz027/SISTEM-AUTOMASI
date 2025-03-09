import sympy as sp

#defini simbol
L, R, V, K, theta, I, doti, J, b, s, t, omega, omegadot = sp.symbols('L R V K theta I doti J b s t omega omegadot')

dIdt = sp.Function('I')(t).diff(t)
dthetadt = sp.Function('theta')(t).diff(t)

# Persamaan listrik motor DC
# L * (dI/dt) + R * I = V - K * (dtheta/dt)
eq1 = sp.Eq(L * dIdt + R * I, V - K * dthetadt)

# Persamaan mekanik motor DC
# J * (d^2theta/dt^2) + b * (dtheta/dt) = K * I
d2thetadt2 = sp.Function('theta')(t).diff(t, t)
eq2 = sp.Eq(J * d2thetadt2 + b * dthetadt, K * I)

# Transformasi Laplace
I_s = sp.Function('I')(s)
omega_s = sp.Function('omega')(s)
V_s = sp.Function('V')(s)

# Transformasi Laplace dari persamaan listrik
eq1_laplace = sp.Eq((L * s + R) * I_s, V_s - K * s * omega_s)

# Transformasi Laplace dari persamaan mekanik
eq2_laplace = sp.Eq(J * s**2 * omega_s + b * s * omega_s, K * I_s)

# Eliminasi I(s)
I_s_expr = sp.solve(eq1_laplace, I_s)[0]
omega_V_relation = sp.solve(eq2_laplace.subs(I_s, I_s_expr), omega_s)[0]

# Fungsi Alih (omega(s) / V(s))
transfer_function = sp.simplify(omega_V_relation / V_s)

# Output persamaan
print("Persamaan Listrik Motor DC:")
print(eq1)
print("\nPersamaan Mekanik Motor DC:")
print(eq2)
print("\nTransformasi Laplace Persamaan Listrik:")
print(eq1_laplace)
print("\nTransformasi Laplace Persamaan Mekanik:")
print(eq2_laplace)
print("\nFungsi Alih (omega(s)/V(s)):")
print(transfer_function)
