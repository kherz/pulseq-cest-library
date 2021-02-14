# sim-library
This folder cointains defined settings for simulations in [yaml](https://yaml.org/)

All relevant parameters are directly obtainable and ce be defines with a simple text editor.


Resources of T1, T2, MTC  and CEST parameters:

|  3T  | T1 [s]      | T2 [ms]     |MTC fraction [%]  |  MTC k [Hz] | MTC T2 [µs]| pub |
| ---- |:-----------:|:-----------:| :---:            |:---:        |:---:         | ---- |
| GM_001| 1.2         |   69        |     5.5/111      |   40        |  9         | [van Zijl 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5650949/)   |
| GM_002| 1.21 ± 0.13 |  71.3 ± 5.6 | (6.3 ± 0.7)/111  |    40 ± 5   | 67 ± 5      | [Heo 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6422734/)    |
| WM_002| 1.05 ± 0.03 |  39.8 ± 3.8 | (11.2 ± 0.7)/111 |    29 ± 4   | 63 ± 4      | [Heo 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6422734/)    |
| GM_003| 1.31 ±0.11  |   71  ± 5   |                  |             |             | [Zhu 2014](https://cds.ismrm.org/protected/14MPresentations/abstracts/3208.pdf)     |
| WM_003| 0.939±0.068 |  62 ± 2     |                  |             |             | [Zhu 2014](https://cds.ismrm.org/protected/14MPresentations/abstracts/3208.pdf)     |
| GM_004| 1.82        |   99        |     0.05         |   40        |  9.1        | [Stanisz 2005](https://doi.org/10.1002/mrm.20605)    |
| WM_004| 1.084       |   69        |     0.139        |   23        |  10         | [Stanisz 2005](https://doi.org/10.1002/mrm.20605)    |


|  7T  | T1 [s]      | T2 [ms]     |MTC fraction [%]  |  MTC k [Hz] | MTC T2 [µs]| pub |
| ---- |:-----------:|:-----------:| :---:            |:---:        |:---:         | ---- |
| GM_001| 1.670 ±0.070 |   43  ± 3   |                  |             |             | [Zhu 2014](https://cds.ismrm.org/protected/14MPresentations/abstracts/3208.pdf)     |
| WM_001| 1.222 ± 0.058 |  37 ± 3     |                  |             |             | [Zhu 2014](https://cds.ismrm.org/protected/14MPresentations/abstracts/3208.pdf)     |


The MT pool is at -2.6 in WM, and at 0 ppm in GM? However unclear how to model, as SuperLorentzian has a pole.
Deviations between MTR T2 can be due to Superlorentzian and Lorentzian linehsape, heursitically we found that when using a Lorentzian R2mt must be 0.2-0.25 smaller to have similar shape.

