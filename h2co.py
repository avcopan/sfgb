import numpy
import sobol


def internal_to_xyz(q):
    Nsample, ndim = numpy.shape(q)
    ndimxyz = 12
    x = numpy.zeros((Nsample, ndimxyz))
    x[:, 4] = internal[:, 0]  # y(O)

    x[:, 6] = -numpy.sin(internal[:, 3]) * internal[:, 1]  # x(H1)
    x[:, 7] = +numpy.cos(internal[:, 3]) * internal[:, 1]  # y(H1)

    x[:, 10] = +numpy.cos(internal[:, 4]) * internal[:, 2]  # y(H2)

    # distance from H2 to the y axis
    ksi = numpy.sin(internal[:, 4]) * internal[:, 2]
    x[:, 9] = ksi * numpy.cos(numpy.pi - internal[:, 5])  # x(H2)
    x[:, 11] = ksi * numpy.sin(numpy.pi - internal[:, 5])  # z(H2)
    value = x

    return value


def xyz_to_rij(x):
    C = x[:, 0:3]
    Ox = x[:, 3:6]
    H1 = x[:, 6:9]
    H2 = x[:, 9:]
    COx = numpy.linalg.norm(C - Ox, axis=1)
    CH1 = numpy.linalg.norm(C - H1, axis=1)
    CH2 = numpy.linalg.norm(C - H2, axis=1)
    OxH1 = numpy.linalg.norm(H1 - Ox, axis=1)
    OxH2 = numpy.linalg.norm(H2 - Ox, axis=1)
    H1H2 = numpy.linalg.norm(H1 - H2, axis=1)
    value = numpy.array([COx, CH1, CH2, OxH1, OxH2, H1H2])
    return value


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

# in ang and degrees, for Vmax = 15000
minmaxinternal = numpy.array([[1.03, 1.50],
                              [0.84, 1.69],
                              [0.84, 1.69],
                              [83, 162],
                              [83, 162],
                              [105, 255]])
# convert angstroms to bohr
minmaxinternal[:3, :] *= 1.88973
# convert degrees to radians
minmaxinternal[3:, :] *= numpy.pi / 180.

for i in range(6):
    internal[:, i] = (numpy.ones((Nsample,)) * minmaxinternal[i, 0] +
                      numpy.ones((Nsample,)) *
                      (minmaxinternal[i, 1] - minmaxinternal[i, 0]) *
                      internal[:, i])

print(internal)
# add equilibrium geometry
internal_eq = numpy.array([r3eq, r1eq, r2eq, theta1eq, theta2eq, phieq])
internal_eq[:3] *= 1.88973

internal = numpy.concatenate(([internal_eq], internal), axis=0)

# convert internals to cartesians
x = internal_to_xyz(internal)
Rij = xyz_to_rij(x)
print(Rij[:, 0] / 1.8873)
