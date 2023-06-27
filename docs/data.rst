=========================
Data and Model Parameters
=========================

The observational data used by the code is as follows. Examples can be
found in the data folder of the repository.

- **time**
The observed times of each burst in your burst train.

- **fluence**
The observed (bolometric) fluence of each burst in your burst train.
Fluence is the integrated burst flux. Expected units are 1e-6 erg/cm^2

- **alpha**
The observed alpha of each burst in your burst train. Alpha is the ratio between fluence and the persistent flux (i.e. nuclear to gravitational energy).

You also need to specify the persistent flux at the time of the bursts.
In "burst train" mode, this information is provided by the file specified
with the ``obsname`` parameter. (The code will interpolate over these
measurements to estimate the appropriate accretion rate to provide to the
ignition code).

In "ensemble" mode the persistent flux at the time of each burst
measurement is also provided in the file specified by the ``burstname``
parameter, and ``obsname`` should be set to ``None``.

Predicted:

These are the model input parameters, for the initial walker positions,
specified by the  *theta* parameter: 

- **X** (or x_0)
Hydrogen mass fraction of the accreted fuel.

- **Z**
CNO metallicity of the accreted fuel.

- **Q_b**
Base heating (MeV/nucleon)

- **f_a** 
Systematic error affecting the alpha values for the likelihood calculation

- **f_E**
Systematic error affecting the fluence values for the likelihood calculation

- **r1** 
Scaling factor relating persistent flux and accretion rate.

- **r2** 
Scaling factor between observed and predicted alphas.

- **r3**
Scaling factor between observed and predicted fluences. 

- **r3**
Scaling factor between observed and predicted fluences. 

- **M**
Neutron star mass, in units of the solar mass

- **R**
Neutron star radius, in km

For equations describing the scaling factors see `Goodwin et al. (2019) <https://doi.org/10.1093/mnras/stz2638>` equations 12, 13, and 14
