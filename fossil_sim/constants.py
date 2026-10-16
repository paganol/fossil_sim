"""Physical and astronomical constants (SI unless the name says otherwise)."""

C_LIGHT_M_S = 299_792_458.0
H_PLANCK_J_S = 6.626_070_15e-34
K_BOLTZMANN_J_K = 1.380_649e-23

# Measured parameters (T_CMB, solar dipole, L2 distance) live in the config.

SECONDS_PER_DAY = 86_400.0
# Julian year, used to convert cycles/yr into rad/s
SECONDS_PER_YEAR = 365.25 * SECONDS_PER_DAY

# 1 MJy/sr in W m^-2 Hz^-1 sr^-1
MJY_SR_IN_SI = 1e-20
