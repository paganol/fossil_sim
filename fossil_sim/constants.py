"""Physical and astronomical constants (SI unless the name says otherwise)."""

C_LIGHT_M_S = 299_792_458.0
H_PLANCK_J_S = 6.626_070_15e-34
K_BOLTZMANN_J_K = 1.380_649e-23

# Fixsen (2009), ApJ 707, 916
T_CMB_K = 2.72548

# Distance between the Earth and the Sun-Earth L2 point; same value as
# litebird_sim (constants.EARTH_L2_DISTANCE_KM)
EARTH_L2_DISTANCE_KM = 1_496_509.30522

# Solar dipole, Planck 2018 I (arXiv:1807.06205)
SOLAR_VELOCITY_KM_S = 369.8160
SOLAR_VELOCITY_GAL_LAT_RAD = 0.842_173_724
SOLAR_VELOCITY_GAL_LON_RAD = 4.608_035_744_4

SECONDS_PER_DAY = 86_400.0
# Julian year, used to convert cycles/yr into rad/s
SECONDS_PER_YEAR = 365.25 * SECONDS_PER_DAY

# 1 MJy/sr in W m^-2 Hz^-1 sr^-1
MJY_SR_IN_SI = 1e-20
