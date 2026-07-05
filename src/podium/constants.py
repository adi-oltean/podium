"""Physical constants, SI throughout (fermi-style single source of truth)."""

# Earth gravitational parameter [m^3/s^2] (EGM2008)
MU_EARTH = 3.986004418e14

# WGS84 equatorial radius [m]
R_EARTH = 6_378_137.0

# Earth J2 zonal harmonic coefficient [-] (EGM2008, unnormalized)
J2_EARTH = 1.08262668e-3

# Earth rotation rate [rad/s] (WGS84)
OMEGA_EARTH = 7.2921159e-5

# Solar-radiation-pressure constant at 1 AU [N/m^2]: the mean solar
# irradiance (~1361 W/m^2) divided by the speed of light.
SOLAR_PRESSURE = 4.5606e-6

# Earth magnetic dipole moment [A m^2] (IGRF epoch ~2020); with
# mu0/4pi = 1e-7 this gives an equatorial surface field ~3.1e-5 T.
EARTH_DIPOLE_MOMENT = 7.94e22
