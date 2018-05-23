import numpy


def PESH2CO(R):
    M, ndim = numpy.shape(R)
    r1eq = 1.10064  # Ang
    r2eq = 1.10064
    r3eq = 1.20296
    theta1eq = 121.65 * numpy.pi / 180.
    theta2eq = 121.65 * numpy.pi / 180.
    phieq = numpy.pi

    # convert from Rij to internal coordinates
    # Rij = [CO, CH1, CH2, OH1, OH2, H1H2]
    r1 = R[:, 1]
    r2 = R[:, 2]
    r3 = R[:, 0]
    costheta1 = (r1 ** 2 + r3 ** 2 - R[:, 3] ** 2) / (2 * r1 * r3)
    theta1 = numpy.arccos(costheta1)
    costheta2 = (r2 ** 2 + r3 ** 2 - R[:, 4] ** 2) / (2 * r2 * r3)
    theta2 = numpy.arccos(costheta2)
    costheta3 = (r1 ** 2 + r2 ** 2 - R[:, 5] ** 2) / (2 * r1 * r2)
    sintheta1 = numpy.sqrt(1 - costheta1 ** 2)
    sintheta2 = numpy.sqrt(1 - costheta2 ** 2)
    cosphi = (costheta3 - costheta1 * costheta2) / (sintheta1 * sintheta2)
    phi = numpy.arccos(cosphi)

    # relative displacements in which the PES is given
    r1 = (r1 - r1eq) / r1
    r2 = (r2 - r2eq) / r2
    r3 = (r3 - r3eq) / r3
    r4 = theta1 - theta1eq
    r5 = theta2 - theta2eq
    r6 = phi - phieq

    print(r1)
    print(r2)
    print(r3)
    print(r4)
    print(r5)
    print(r6)
