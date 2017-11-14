clear;clc;clf

r_1 = 0.0381; %[m]
r_2 = 0.04; %[m]
k_steel = 14.6; %k of steel
h = 342.1531; %h is a normally a complex function. It is assumed to be a constant here
%h is evaluated at 0.18 kg/s, 80 C
k_insulation = 0.206 + (7.702e-4)* 80; %evauluated at 80 C
r_insulation_thickness = 0.05;
R_insul = (log((r_insulation_thickness+r_2)/r_2))/(k_insulation);
U = R_insul^-1;
V_total = (r_2^2 - r_1^2)*pi*2;
S = 10000/(V_total); %This should be the volumetric source rate (W/m3). 
T_fluid = 80;
T_air = 20;

A = [((1/(r_1*h) + log(r_1))) 1; ((1/(r_2*U)+log(r_2))) 1];
B = [(-S*r_1/2*h) + (-S*(r_1^2)/4*k_steel) - T_fluid;
    (-S*r_2/2*U) + (-S*(r_2^2)/4*k_steel) - T_air
    ];
C = A\-B;
temp_profile = @(r) (-S.*(r.^2)./(4.*k_steel)) + C(1).*log(r) + C(2);
r_profile = linspace(r_1,r_2,20);

plot(r_profile,temp_profile(r_profile))

