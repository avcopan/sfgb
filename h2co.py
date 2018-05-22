import numpy
import sobol


ndimxyz = 12

Vmax = 15000
M = 15000
Nsets = 5
Ng = 15000
LHS = numpy.zeros((Ng, Ng), dtype=numpy.float64)
RHS = numpy.zeros((Ng, Ng), dtype=numpy.float64)
Nsample = 2000000

print(LHS.shape)
print(RHS.shape)

# equilibrium geometry in angstroms
r1eq = 1.10064
r2eq = 1.10064
r3eq = 1.20296
theta1eq = 121.65 * numpy.pi / 180.
theta2eq = 121.65 * numpy.pi / 180.
phieq = numpy.pi

# Cartesian equilibrium coordinates in a.u.
Ceq = numpy.array([0.0, 0.0, 0.0]) * 1.88973
Oeq = numpy.array([0.0, r3eq, 0.0]) * 1.88973
H1eq = numpy.array([-numpy.sin(theta1eq) * r1eq,
                    numpy.cos(theta1eq) * r1eq, 0.0]) * 1.88973
H2eq = numpy.array([numpy.sin(theta2eq) * r2eq,
                    numpy.cos(theta2eq) * r2eq, 0.0]) * 1.88973
xeq = numpy.array([Ceq, Oeq, H1eq, H2eq])
m = numpy.array([12.0 * numpy.ones((1, 3)),  # masses
                 15.99491 * numpy.ones((1, 3)),
                 1.00794 * numpy.ones((1, 3)),
                 1.00794 * numpy.ones((1, 3))])

# Rij equilibrium coordinates in a.u.
COeq = numpy.linalg.norm(Ceq - Oeq)
CH1eq = numpy.linalg.norm(Ceq - H1eq)
CH2eq = numpy.linalg.norm(Ceq - H2eq)
OH1eq = numpy.linalg.norm(Oeq - H1eq)
OH2eq = numpy.linalg.norm(Oeq - H2eq)
H1H2eq = numpy.linalg.norm(H1eq - H2eq)
Req = numpy.array([COeq, CH1eq, CH2eq, OH1eq, OH2eq, H1H2eq])

try:
    print('try')
    internal = numpy.load('sobseq.npy')
except FileNotFoundError:
    print('except')
    internal = sobol.i4_sobol_generate(6, Nsample)
    numpy.save('sobseq', internal)
print(internal)
