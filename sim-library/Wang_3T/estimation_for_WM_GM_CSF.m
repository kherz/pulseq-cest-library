% estimation of Wang model for GM
%%

T1_zhu_WM= [0.939]      ;                 % [Zhu 2014] http://hdl.handle.net/21.11116/0000-0001-32FE-9
T1_zhu_GM= [1.310]      ;                 % [Zhu 2014] http://hdl.handle.net/21.11116/0000-0001-32FE-9
fmatter=T1_zhu_GM/T1_zhu_WM;

%% WM
R1A=0.4; R1C=12.2./3;
fc=0.289; kac = 1.38; kca= 1.38/fc;
R1obs_wang_approx = (R1A + fc*R1C )/(1 + fc);
R1obs_wang =    0.5*( kac + kca + R1A+ R1C - sqrt(( kac + kca + R1A + R1C ).^2 - 4*( kca*R1A + kac.*R1C + R1A.*R1C )));
T1WM = 1./R1obs_wang
T1WM*fmatter

%% GM
R1A=0.4; R1C=12.2./3;
fc=0.1039; kac = 1.38; kca= 1.38/fc;
R1obs_wang_approx = (R1A + fc*R1C )/(1 + fc);
R1obs_wang =    0.5*( kac + kca + R1A+ R1C - sqrt(( kac + kca + R1A + R1C ).^2 - 4*( kca*R1A + kac.*R1C + R1A.*R1C )));
T1GM = 1./R1obs_wang
