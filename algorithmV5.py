import math
# Define inputs
pi = math.pi
V = float(input("speed of the aircraft in km/h")) / 3.6 #speed of the aircraft
h = float(input("altitude of the aircraft in meters")) #altitude
rho = 0.4 # Air density at altitude in kg/m^3/ would be picked up by a sensor
alpha = float(input("angle of attack of the aircraft in degres")) * pi / 180 #angle of attack of the aircraft
W = 50000 # Weight of the aircraft in N/ specific to the aircraft
S = 100 # Wing area in m^2/ specific to the aircraft
q = 0.5 * rho * V**2 # Dynamic pressure in Pa
T= 224.15 # temperature at 10 300 m as an example, would be picked up by a sensor

# Define objective function
# Minimize drag coefficient CD
def objective(x):
    # x is a vector of wing geometry parameters
    # x[0] = sweep angle Lambda
    # x[1] = span b
    # x[2] = chord c
    # x[3] = camber epsilon

    # Calculate aspect ratio AR
    AR = x[1] / x[2]
    k = 0.05 # k is a constant that depends on the wing planform shape and taper ratio
    e = 1 / (1 + k * AR) # e is the Oswald efficiency factor

    # Calculate Mach number M
    M = V / a # a is the speed of sound at altitude

    # Calculate Reynolds number Re
    Re = rho * V * x[2] / mu # mu is the air viscosity at altitude

    # Calculate lift coefficient CL
    CL = W / (q * S)

    # Calculate drag coefficient CD using empirical formulas
    # CD = CD0 + CDi + CDw + CDv

    # CD0 is the zero-lift drag coefficient
    CD0 = f(x, M, Re) # f is a function based on arbitrary data

    # CDi is the induced drag coefficient
    CDi = CL**2 / (pi * e * AR)

    # CDw is the wave drag coefficient
    CDw = g(x, M) # g is a function based on arbitrary data

    # CDv is the viscous drag coefficient
    CDv = h(x, Re) # h is a function based on arbitrary data

    # Return CD as the objective value
    return CD0 + CDi + CDw + CDv


# Define bounds for wing geometry parameters
# The bounds are based on physical limits or design requirements
sweep_min = 0
sweep_max = pi/4
b_min = 10
b_max = 50
c_min = 1
c_max = 5
epsilon_min = 0
epsilon_max = 0.1

bnds = ((sweep_min, sweep_max), (b_min,b_max), (c_min,c_max), (epsilon_min,epsilon_max))

# Define constraints
# The constraints are based on physical limits or design requirements
constraint1 = {'type': 'ineq', 'fun': lambda x: pi/4 - x[0]}
constraint2 = {'type': 'ineq', 'fun': lambda x: x[1] - b_min}
constraint3 = {'type': 'ineq', 'fun': lambda x: b_max - x[1]}
constraint4 = {'type': 'ineq', 'fun': lambda x: x[2] - c_min}
constraint5 = {'type': 'ineq', 'fun': lambda x: c_max - x[2]}
constraint6 = {'type': 'ineq', 'fun': lambda x: x[3] - epsilon_min}
constraint7 = {'type': 'ineq', 'fun': lambda x: epsilon_max - x[3]}

# Define initial guess for wing geometry parameters
# The initial guess is be based random values

x0 = (pi/8,25,2,0.05)


# Define additional parameters that are needed for the objective function and constraints

# Define the speed of sound at altitude a using a formula
R = 287
a = math.sqrt(1.4 * R * T) # R is the gas constant for air (287 J/(kg.K)), T is the temperature at altitude (K)

# Define the temperature at altitude T using a formula
T_0 = 288.15
L = 0.0065
T = T_0 - L * h # T_0 is the reference temperature at sea level (288.15 K), L is the temperature lapse rate (0.0065 K/m), h is the altitude (m)

# Define the air viscosity at altitude mu using a formula
mu_0 = math.pow(1.7894,-5) # mu_0 is the reference viscosity at sea level (1.7894e-5 kg/(m s))
Sconst = 110.4 # Sconst is a constant (110.4 K)
mu = mu_0 * (T / T_0)**(3/2) * (T_0 + S) / (T + Sconst)  #T is the temperature at altitude (K)

# Define the functions f, g, and h that calculate CD0, CDw, and CDv using formulas
def f(x, M, Re):
    # x[0] = sweep angle Lambda
    # x[2] = chord c
    # x[3] = camber epsilon
    # M = Mach number
    # Re = Reynolds number
    return 0.02 + 0.005 * math.cos(x[0]) + 0.001 * x[2] + 0.002 * x[3] + 0.0001 * M + 0.00001 * Re

def g(x, M):
    # x[0] = sweep angle Lambda
    # x[2] = chord c
    # x[3] = camber epsilon
    # M = Mach number
    return 0.001 * math.sin(x[0]) + 0.0005 * x[2] + 0.0002 * x[3] + 0.01 * (M - 1)**2

def h(x, Re):
    # x[0] = sweep angle Lambda
    # x[2] = chord c
    # x[3] = camber epsilon
    # Re = Reynolds number
    return 0.003 * math.cos(x[0]) + 0.0002 * x[2] + 0.0001 * x[3] + 0.000001 * Re

# Import optimization library
from scipy.optimize import minimize

# Define optimization method and options
method = 'SLSQP' # Sequential Least Squares Programming
options = {'disp': True} # Display convergence information


# Perform optimization
result = minimize(objective, x0, method=method, bounds=bnds,
                  constraints=[constraint1, constraint2, constraint3,
                               constraint4, constraint5, constraint6,
                               constraint7],
                  options=options)

# Output optimal wing geometry parameters and drag coefficient
print(result.x)
print(result.fun)
