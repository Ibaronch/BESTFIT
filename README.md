# BESTFIT
IDL based SED fitting code for stellar templates
V1.1
Ivano Baronchelli 2018
BESTFIT is an IDL based function that can be used to fit photometric data using a collection of
template models of stellar spectra. The method used is based on the χ 2 minimization technique.
Chosing among a collection of models (internal or supplied by the users), BESTFIT computes
the spectral model that better fits the input photometric data. BESTFIT also predicts the
expected flux in a users-selected photometric filter (or at a given wavelength). This can be
useful in all those cases in which a rough photometrical calibration of the fluxes in a given
band is needed or when the user wants to check the photometric calibration in a given band.
The parameters and filters files look similar to those used by the hyperz software (Bolzonella+2000).
This should make the use easier for people that previously used hyperz.

For installation instructions and manual see the BESTFIT pdf manual included in this folder 
