# sim-library
This folder cointains defined settings for simulations in [yaml](https://yaml.org/)

All relevant parameters are directly obtainable and ce be defines with a simple text editor.


Resources of T1, T2, MTC  and CEST parameters:

|  3T  | T1 [s]      | T2 [ms]     |MTC fraction [%]  |  MTC k [Hz] | MTC T2 [µs]| pub |
| ---- |:-----------:|:-----------:| :---:            |:---:        |:---:         | ---- |
| GM   | 1.21 ± 0.13 |  71.3 ± 5.6 | (6.3 ± 0.7)/111  |    40 ± 5   | 67 ± 5      | [Heo 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6422734/)    |
| WM   | 1.05 ± 0.03 |  39.8 ± 3.8 | (11.2 ± 0.7)/111 |    29 ± 4   | 63 ± 4      | [Heo 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6422734/)    |
| GM   | 1.31 ±0.11  |   71  ± 5   |                  |             |             | [Zhu 2014](https://cds.ismrm.org/protected/14MPresentations/abstracts/3208.pdf)     |
| WM   | 0.939±0.068 |  62 ± 2     |                  |             |             | [Zhu 2014](https://cds.ismrm.org/protected/14MPresentations/abstracts/3208.pdf)     |
| GM   | 1.82        |   99        |     0.05         |   40        |  9.1        | [Stanisz 2005](https://doi.org/10.1002/mrm.20605)    |
| WM   | 1.084       |   69        |     0.139        |   23        |  10         | [Stanisz 2005](https://doi.org/10.1002/mrm.20605)    |
| GM   | XX          |   XX        |     XX        |   XX        |  XX         | van Zijl Neuroimage   |

 case 'WM'  % stanisz 3T
        T1=1.084;
        T2=0.069;
        Sim.dwA=0.0;
        Sim.R1A=1/T1;
        Sim.R2A=1/T2;
        
        %MT
        Sim.R1C=1/T1;
        Sim.fC=0.139;
        Sim.R2C=100000;
        Sim.kCA=23;
        Sim.kAC=Sim.kCA*Sim.fC;
        Sim.dwC=0;
        
    case 'GM' % stanisz 3T
        T1=1.820;
        T2=0.099;
        Sim.dwA=0;
        Sim.R1A=1/T1;
        Sim.R2A=1/T2;
        
        %MT
        Sim.R1C=1/T1;
        Sim.fC=0.05;
        Sim.R2C=109890;  % 1/9.1µs
        Sim.kCA=40;
        Sim.kAC=Sim.kCA*Sim.fC;
        Sim.dwC=-2.6;             %asymmetric MT ?
