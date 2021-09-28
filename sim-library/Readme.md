# pulseq-cest: sim-library
The simulation parameters used in the [pulseq-cest project](https://pulseq-cest.github.io/) are stored in the human
readable [yaml](https://yaml.org/) file format. These yaml/simulation files can be read and written using a
simple text editor. 


## Resources of T1, T2, MTC  and CEST parameters:

|  3T  | T1 [s]      | T2 [ms]     |MTC fraction [%]  |  MTC k [Hz] | MTC T2 [µs]| pub |
| ---- |:-----------:|:-----------:| :---:            |:---:        |:---:         | ---- |
| GM_001| 1.2        |   69        |     5.5/111      |    40       | 9         | [van Zijl 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5650949/)   |
| GM_vanZijl2018| 1.2         |   69        |     5.5/111      |   40        |  9         | [van Zijl 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5650949/)   |
| WM_001| 1.05       |  39.8       |    15.42/111      |    29       | 9/0.23         | [Heo 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4662920/) |
| GLIO_002| 1.5       |  70       |    6/111      |    29       | 9/0.23         | WM_001 adjusted to match [Heo 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4662920/) |



## Detailed information about the .yaml files
Here we give a detailed overview over the structure of the .yaml files and the parameters that must or can be defined. As
an example, the [GM_3T_001_bmsim.yaml](GM_3T_001_bmsim.yaml) file is presented and explained.

### .yaml file header
The header gives an overview over the pools that are defined in the file, some helpful information about the references
for the chosen parameters as well as some general information about the different pool settings that can be adjusted:
```
####################################################################################
# ------ Simulation settings for (py)pulseq-cest BlochMcConnell simulations ------ #
# ------------------------ https://pulseq-cest.github.io/ ------------------------ #
####################################################################################
# Simulation parameters for gray matter at 3T with
# - 4 CEST pools
# - 1 NOE pool      
# - a SuperLorentzian shaped MT pool
#
# sources/references for settings in this file:
# CEST pools: https://doi.org/10.1016/j.neuroimage.2017.04.045

###### Pool settings
#         f: relative pool size fraction (float)
#        t1: longitudinal relaxation time T1 = 1/R1 [s] (float)
#        t2: transversal relaxation time T2 = 1/R2 [s] (float)
#         k: exchange rate [Hz] to water pool (MT and CEST pools only) (float)
#        dw: chemical shift (delta omega) relative to water [ppm] (MT and CEST pools only) (float)
# lineshape: lineshape of the MT Pool ('Lorentzian', 'SuperLorentzian' or 'None') (str)
```

### water pool (mandatory)
The **water pool** is mandatory for all pulseq-cest based simulations. It's parameters and values are usually defined 
directly under the header:
```
water_pool: {
  f: 1.0,
  t1: 1.2,
  t2: 0.069
}
```
The (exactly) three parameters that have to be defined are:
* **f:** the (relative) proton fraction (usually defined as 1 for the water pool) (float)
* **t1:** the longitudinal relaxation time T1 = 1/R1 [s] (float)
* **t2:** the transversal relaxation time T2 = 1/R2 [s] (float)

For relaxation times, it is also possible to use the inverse relaxation rates R1 and R2. The same **water pool** could therefore be defined with: 

```
water_pool: {
  f: 1.0,
  r1: 0.8333,
  r2: 14.4928
}
```

### CEST pool(s) (optional)
An arbitrary number of pools can be defined in the *cest_pool* section/dictionary. For every **CEST pool**, a unique 
name and its parameters have to be defined: 
```
cest_pool: {
  'amide': {
    f: 0.00064865, # 72e-3 / 111,
    t1: 1.3,
    t2: 0.1,
    k: 30,
    dw: 3.5
  },
  'amine': {
    f: 0.00018018, # 20e-3 / 111,
    t1: 1.3,
    t2: 0.2,
    k: 5500,
    dw: 3
  },
  'guanidinium': {
    f: 0.00018018, # 20e-3 / 111,
    t1: 1.3,
    t2: 0.17,
    k: 1100,
    dw: 2
  },
  'hydroxyl': {
    f: 0.00040541, #45e-3 / 111,
    t1: 1.3,
    t2: 0.055,
    k: 2000,
    dw: 0.9
  },
  'NOE': { #average of 5 NOE pools from https://doi.org/10.1016/j.neuroimage.2017.04.045
    f: 0.0045, #500e-3 / 111,
    t1: 1.3,
    t2: 0.005,
    k: 16,
    dw: -2.75
  }
}
```
Additionally to the three parameters (*f*, *r1*, *r2*), which have to be defined for the **water pool**, every 
**CEST pool** requires two more parameters. These are:
* dw: the chemical shift relative to water [ppm] (float)
* k: the exchange rate to the water pool [Hz] (float)

### MT pool (optional)
The (macromolecular) Magnetization Transfer (MT) pool is an optional pool:
```
mt_pool: {
  f: 0.0495, # 5.5 / 111,
  t1: 1.3,
  t2: 9.0e-6,
  k: 40,
  dw: 0,
  lineshape: 'SuperLorentzian'
}
```
Additionally to the five parameters (*f*, *t1*, *t2*, *k*, *dw*), which have to be defined for every 
**CEST pool**, the **MT pool** further requires the *lineshape* parameter:
* lineshape: string with MT lineshape. Can be either *'Lorentzian'* or *'SuperLorentzian'*.

### Scanner/field properties (partially mandatory)
Additionally to the pool settings, some scanner/field properties have to be defined:
```
b0: 3
gamma: 267.5153
b0_inhom: 0.0
rel_b1: 1
```
Please note the default values for *gamma*, *b0_inhom* and *rel_b1* that are used if you don't define them:
* b0: magnetic field strength [T] (float)
* gamma: gyromagnetic ratio [rad / µT], default is 267.5153 (float)
* b0_inhom: B0 field inhomogeneity [ppm], default is 0.0 (float)
* rel_b1: relative B1 inhomogeneity/factor, default is 1.0 (float)


### Further optional parameters
```
verbose: False
```
* verbose: True to enable additional output from the C++ code. Default is False (bool)
```
reset_init_mag: true
```
* reset_init_mag: True if magnetization should be reset after each ADC. Default is True (bool)
```
scale: 0.5
```
* scale: relative value (range [0, 1]) to scale the magnetization if `reset_init_mag: True`. Default is 1 (float)
```
max_pulse_samples: 200
```
* max_pulse_samples: defines the number of samples used for the simulation of shaped pulses. Default is 500 (int)
