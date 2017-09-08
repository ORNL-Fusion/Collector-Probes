# Program: meas_locations_v3.py
# Author:  Shawn Zamperini
# Email:   zamp@utk.edu
# Date:    4/30/17
#
# The variables in the program are pulled from the slides
# named Probe_Geo_Figures. That naming convention is followed here.
#
# New in v3: Cleaned up v2_1 with better comments and docstrings.
#			 Functions now accept the location along probe to be calculated.
#			 Functions now return the radial measurement location.
#            Replaced r_wall with r_offset.
#
# New in v3_1: Was missing the degrees() for f calculation on line 51.
#
# New in v3_2: Changed naming convention from L,R to D,U, respectively.
#			   Made sure no variables were reused. All of them match the slides.
#              Fixed calculation of r_BU and r_CU from n+q+q to 360-(n+q+q).
#
# New in v3_3: Changed beta values.


# Needed for trig functions.
from math import *


def calc_R_measAD(r_probe, location):
	"""Calculate the radial position of a measurement along the left A probe."""

	# Probe arm inserted at angle to centerline of port axis.
	offset = 13.0

	# Inches to centimeters (actually where 13 degrees happens, not the wall)
	r_offset = 113.46 * 2.54

	# Geometry of A probe holder.
	alpha = (5.125 - 4.920) * 2.54
	beta  = (0.443 + 0.5 * 0.256) * 2.54
	delta = sqrt(alpha**2 + beta**2)

	# Angle between r_probe and A probe center axis. 180 degrees since we
	# want the higher value in quadrant II (basic trig stuff).
	c = 180 - degrees(asin(r_offset * sin(radians(offset)) / r_probe))

	# Angle between alpha and beta
	d = degrees(atan(beta / alpha))

	# Angle between r_probe and delta
	e = c - d

	# Radial position of A-D tip using law of cosines
	r_AD = sqrt(delta**2 + r_probe**2 - 2 * r_probe * delta * cos(radians(e)))

	# Angle between tip of left A probe and length of A-D probe. Again want
	# the higher value in quadrant II.
	f = 180 - degrees(asin(r_offset * sin(radians(13.0)) / r_AD))

	# Radial position of measurement location, l, along A-D. l=0 is at
	# the tip of the probe. In cm.
	r_ADmeas = sqrt(location**2 + r_AD**2 - 2 * location * r_AD * cos(radians(f)))

	return r_ADmeas


def calc_R_measAU(r_probe, location):
	"""
	Calculate the radial position of a measurement along the right A probe.
	Inputs need to be in cm? 09/08/2017 EAU
	"""

	# Already defined terms.
	offset = 13
	r_offset = 113.46 * 2.54
	alpha = (5.125 - 4.920) * 2.54
	beta  = (0.443 + 0.5 * 0.256) * 2.54
	delta = sqrt(alpha**2 + beta**2)
	c = 180 - degrees(asin(r_offset * sin(radians(offset)) / r_probe))
	d = degrees(atan(beta / alpha))

	# Angle between r-probe and delta
	k = 360 - d - c

	# Radial position of the tip of the right A probe.
	r_AU = sqrt(delta**2 + r_probe**2 - 2 * r_probe * delta * cos(radians(k)))

	# Angle between length of probe and radial position vector of the tip.
	m = 180 - degrees(asin(r_offset * sin(radians(13)) / r_AU))

	r_AUmeas = sqrt(location**2 + r_AU**2 - 2 * location * r_AU * cos(radians(m)))

	return r_AUmeas

def calc_R_measBD(r_probe, location):
	"""Calculate the radial position of a measurement along the left B probe."""

	# Already defined terms.
	offset = 13.0
	r_offset = 113.46 * 2.54
	alpha = 0.0001 * 2.54
	beta = 0.25 * 2.54
	delta = sqrt(alpha**2 + beta**2)
	c = 180.0 - degrees(asin(r_offset * sin(radians(offset)) / r_probe))

	# Follows same logic as A probes. Refer to slides.
	r_Btip = sqrt(r_probe**2 + beta**2 - 2 * r_probe * beta * cos(radians(c)))
	f = degrees(asin(r_probe * sin(radians(c)) / r_Btip))
	q = degrees(atan(beta / alpha))
	n = 180.0 - f - q
	r_BD = sqrt(delta**2 + r_Btip**2 - 2 * delta * r_Btip * cos(radians(n)))
	v = 180.0 - degrees(asin(r_offset * sin(radians(13)) / r_BD))

	r_BDmeas = sqrt(location**2 + r_BD**2 - 2 * location * r_BD * cos(radians(v)))

	return r_BDmeas


def calc_R_measBU(r_probe, location):
	"""Calculate the radial position of a measurement along the right B probe."""

	# Already defined terms.
	offset = 13.0
	r_offset = 113.46 * 2.54
	alpha = 0.0001 * 2.54
	beta = 0.25 * 2.54
	delta = sqrt(alpha**2 + beta**2)
	c = 180.0 - degrees(asin(r_offset * sin(radians(offset)) / r_probe))
	r_Btip = sqrt(r_probe**2 + beta**2 - 2 * r_probe * beta * cos(radians(c)))
	f = degrees(asin(r_probe * sin(radians(c)) / r_Btip))
	q = degrees(atan(beta / alpha))
	n = 180.0 - f - q

	# Refer to slides for picture.
	r_BU = sqrt(r_Btip**2 + delta**2 - 2 * r_Btip * delta * cos(radians(360 - (n + q + q))))
	p = 180.0 - degrees(asin(r_offset * sin(radians(13)) / r_BU))

	r_BUmeas = sqrt(location**2 + r_BU**2 - 2 * location * r_BU * cos(radians(p)))

	return r_BUmeas


def calc_R_measCD(r_probe, location):
	"""Calculate the radial position of a measurement along the left C probe."""

	# Same exact thing as BL, only different beta value.
	offset = 13.0
	r_offset = 113.46 * 2.54
	alpha = 0.0001 * 2.54
	beta = 0.15 * 2.54
	delta = sqrt(alpha**2 + beta**2)
	c = 180.0 - degrees(asin(r_offset * sin(radians(offset)) / r_probe))
	r_Ctip = sqrt(r_probe**2 + beta**2 - 2 * r_probe * beta * cos(radians(c)))
	f = degrees(asin(r_probe * sin(radians(c)) / r_Ctip))
	q = degrees(atan(beta / alpha))
	n = 180.0 - f - q
	r_CD = sqrt(delta**2 + r_Ctip**2 - 2 * delta * r_Ctip * cos(radians(n)))
	v = 180.0 - degrees(asin(r_offset * sin(radians(13)) / r_CD))

	r_CDmeas = sqrt(location**2 + r_CD**2 - 2 * location * r_CD * cos(radians(v)))

	return r_CDmeas


def calc_R_measCU(r_probe, location):
	"""Calculate the radial position of a measurement along the right C probe."""

	# Same exact thing as BR, only different beta value.
	offset = 13.0
	r_offset = 113.46 * 2.54
	alpha = 0.0001 * 2.54
	beta = 0.15 * 2.54
	delta = sqrt(alpha**2 + beta**2)
	c = 180.0 - degrees(asin(r_offset * sin(radians(offset)) / r_probe))
	r_Ctip = sqrt(r_probe**2 + beta**2 - 2 * r_probe * beta * cos(radians(c)))
	f = degrees(asin(r_probe * sin(radians(c)) / r_Ctip))
	q = degrees(atan(beta / alpha))
	n = 180.0 - f - q
	r_CU = sqrt(r_Ctip**2 + delta**2 - 2 * r_Ctip * delta * cos(radians(360 - (n + q + q))))
	p = 180.0 - degrees(asin(r_offset * sin(radians(13)) / r_CU))

	r_CUmeas = sqrt(location**2 + r_CU**2 - 2 * location * r_CU * cos(radians(p)))

	return r_CUmeas
