"""
The default cosmological constants in this package are the values derived from 7 year WMAP data. You can define new values using the following functions.
"""

# This file contains cosmological constants, physical constants and conversion factors.
import numpy as np
from astropy.cosmology import Cosmology, FLRW


def set_hubble_h(value):
	"""
	Define new hubble constant value (little h).
	"""
	from tools21cm import cosmology as cm

	global h, H0
	h  = value
	H0 = 100.0*h

	cm.precalc_table_cdist = cm.z_to_cdist(cm.precalc_table_z)
	# cm.recalculate_precalc_cdists()


def set_omega_matter(value):
	"""
	Define new omega matter value.
	"""
	from . import cosmology as cm

	global Omega0, lam
	Omega0 = value
	lam = 1.0 - Omega0

	cm.precalc_table_cdist = cm.z_to_cdist(cm.precalc_table_z)


def set_omega_baryon(value):
	"""
	Define new omega baryon value.
	"""
	from . import cosmology as cm

	global OmegaB
	OmegaB = value

	cm.precalc_table_cdist = cm.z_to_cdist(cm.precalc_table_z)


def set_omega_lambda(value):
	"""
	Define new omega lambda value.
	"""
	from . import cosmology as cm

	global lam, OmegaL
	lam = value
	OmegaL = value

	cm.precalc_table_cdist = cm.z_to_cdist(cm.precalc_table_z)


def set_ns(value):
	"""
	Define new ns value.
	"""
	from . import cosmology as cm

	global n_s
	n_s = value

	cm.precalc_table_cdist = cm.z_to_cdist(cm.precalc_table_z)


def set_sigma_8(value):
	"""
	Define new sigma_8 value.
	"""
	from . import cosmology as cm

	global sigma_8
	sigma_8 = value

	cm.precalc_table_cdist = cm.z_to_cdist(cm.precalc_table_z)

def set_abundance_helium(value):
	global abu_he_mass, abu_h_mass
	abu_he_mass = value
	abu_h_mass = 1.0-abu_he_mass


def set_cosmology(cosmology: FLRW):
	"""
	Sets tools21cm's cosmology parameters to those defined by a cosmology
	instance from astropy
	"""
	set_hubble_h(cosmology.H0.value / 100.)
	set_omega_matter(cosmology.Om0)
	set_omega_baryon(cosmology.Ob0)
	set_omega_lambda(1.0 - cosmology.Om0)
	set_ns(cosmology.meta['n'])
	set_sigma_8(cosmology.meta['sigma8'])


# Various useful physical constants
abu_he = 0.074
abu_h = 1.0-abu_he
c = 3.0e5 # km/s
pc=  3.086e18 #1 pc in cm
Mpc = 1e6*pc
G_grav = 6.6732e-8
m_p = 1.672661e-24 #g
mean_molecular = abu_h+4.0*abu_he
abu_he_mass = 0.2486 
abu_h_mass = 1.0-abu_he_mass

mean_molecular = 1.0/(1.0-abu_he_mass)
solar_masses_per_gram = 5.02785431e-34
kms = 1.e5 #1 km/s in cm/s

# Cosmology
h = 0.7
Omega0 = 0.27
OmegaB = 0.044
lam = 1.0 - Omega0
OmegaL = 1.0 - Omega0
n_s = 0.96
sigma_8 = 0.8

# ES : Update for Planck 2018 for 21cmFast input
# h = 0.6766
# Omega0 = (0.02242 + 0.11933) / 0.6766 ** 2.
# OmegaB = 0.02242 / 0.6766 ** 2.
# lam = 1.0-Omega0; OmegaL = lam
# n_s = 0.9665
# sigma_8 = 0.8102

# Cosmology
H0 = 100.0*h
H0cgs = H0*1e5/Mpc
rho_crit_0 = 3.0*H0cgs*H0cgs/(8.0*np.pi*G_grav)
q0 = 0.5*Omega0- lam
rho_matter = rho_crit_0*Omega0  
Tcmb0=2.725

# Redshift dependent Hubble parameter, km/s/Mpc
Hz = lambda z: H0*np.sqrt(Omega0*(1.0+z)**3.+lam) 

# 21 cm stuff
A10 = 2.85e-15
nu0 = 1.42e3
Tstar = 0.068
lambda0 = c*1.0e5/(nu0*1.0e6) # cm
num_h_0=(1-abu_he_mass)*OmegaB*rho_crit_0/m_p
# meandt = 2.9*0.043/0.04
meandt = 3.0*lambda0**3/(32.*np.pi)*A10*Tstar*num_h_0/(H0cgs/h)*1000.
