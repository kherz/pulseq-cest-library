# sim-library
This folder contains defined settings for simulations in [yaml](https://yaml.org/)

All relevant parameters are directly obtainable and can be defines with a simple text editor.


Resources of T1, T2, MTC  and CEST parameters:

|  3T  | T1 [s]      | T2 [ms]     |MTC fraction [%]  |  MTC k [Hz] | MTC T2 [µs]| pub |
| ---- |:-----------:|:-----------:| :---:            |:---:        |:---:         | ---- |
| GM_001| 1.2        |   69        |     5.5/111      |    40       | 9         | [van Zijl 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5650949/)   |
| GM_vanZijl2018| 1.2         |   69        |     5.5/111      |   40        |  9         | [van Zijl 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5650949/)   |
| WM_001| 1.05       |  39.8       |    15.42/111      |    29       | 9/0.23         | [Heo 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4662920/) |
| GLIO_002| 1.5       |  70       |    6/111      |    29       | 9/0.23         | WM_001 adjusted to match [Heo 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4662920/) |
