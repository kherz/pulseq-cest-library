%Install pulseq CEST (READ.me) and add path
addpath(genpath('...\pulseq_cest'))

%% MRF_CEST_3T_Gaussian_004
cd('....\pulseq-cest-library\MRF_CEST_3T_Gaussian_004')
Ma=simulate_pulseqcest('.....\pulseq-cest-library\MRF_CEST_3T_Gaussian_004\MRF_CEST_3T_Gaussian_004.seq','...\pulseq-cest-library\sim-library\WM_3T_default_7pool_bmsim.yaml')

offsets_ppm=[3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5] 

figure;
plot(1:numel(offsets_ppm),Ma)

%% MRF_CEST_3T_LOAS_MTC_005
cd('.....\pulseq-cest-library\MRF_CEST_3T_LOAS_MTC_005')
Mb1=simulate_pulseqcest('.....\pulseq-cest-library\MRF_CEST_3T_LOAS_MTC_005\MRF_CEST_3T_LOAS_MTC_005_10.seq','...\pulseq-cest-library\sim-library\WM_3T_default_7pool_bmsim.yaml')
Mb2=simulate_pulseqcest('.....\pulseq-cest-library\MRF_CEST_3T_LOAS_MTC_005\MRF_CEST_3T_LOAS_MTC_005_40.seq','...\pulseq-cest-library\sim-library\WM_3T_default_7pool_bmsim.yaml')

offsets_ppm1=[9.1 8.9 50 50 11.5 8.6 20.7 50 50 50]
offsets_ppm2=[50 10.2 23.6 15.6 11.3 15.2 8.9 23.2 28.4 50 28.6 19.9 11.9 11.7 9.7 8.7 9 50 49.9 35.6 50 50 10.1 50 10.1 8.7 35.1 9.4 49.8 34.8 10.4 50 50 18.9 8.1 24.2 10.8 9.6 50 50]

figure;
plot(1:numel(offsets_ppm1),Mb1)
figure;
plot(1:numel(offsets_ppm2),Mb2)


%% MRF_CEST_MTC_3T_UNSUPER_006
cd('.....\pulseq-cest-library\MRF_CEST_MTC_3T_UNSUPER_006')
Mc=simulate_pulseqcest('.....\pulseq-cest-library\MRF_CEST_MTC_3T_UNSUPER_006\MRF_CEST_3T_MTC_UNSUPER_006.seq','...\pulseq-cest-library\sim-library\WM_3T_default_7pool_bmsim.yaml')

offsets_ppm=[8 8 8 8 9 9 9 9 10 10 10 10 11 11 11 11 13 13 13 13 15 15 15 15 20 20 20 20 25 25 25 25 35 35 35 35 50 50 50 50]

figure;
plot(1:numel(offsets_ppm),Mc)

%% MMRF_CEST_3T_SCONE_007
cd('.....\pulseq-cest-library\MRF_CEST_3T_SCONE_007')
Md=simulate_pulseqcest('.....\pulseq-cest-library\MRF_CEST_3T_SCONE_007\MRF_CEST_3T_SCONE_007.seq','...\pulseq-cest-library\sim-library\WM_3T_default_7pool_bmsim.yaml')
offsets_ppm=[3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5] 

figure;
plot(1:numel(offsets_ppm),Md)

