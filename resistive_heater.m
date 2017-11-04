%%Calculate 1D axial heater. Each segment acts as a homogeneous control
%%volume. The heater is divided into n_segments, and seperated into the
%%fluid portion and the resistive-heater portion
clear;clc;clf
n_segments = 10;

%input resistive heater dimensions
total_x = 2; %[m]
x_step = total_x/n_segments;
x_profile = linspace(0,total_x,n_segments);
d2T_dx2 = ones(n_segments,1);

%Time step
time_end = 1000; %s
t_segments = 1;
t_step = time_end/t_segments;
t_profile = linspace(0,time_end,t_segments);

%Store all temperatures in both portion in a matrix [resister_T fluid_T];
%input initial temperature profile in both portions
homogeneous_temp_heater_side = 323; %[K]
homogeneous_temp_fluid_side = 323; %[K]
T_heater_initial = ones(n_segments,1).*homogeneous_temp_heater_side;
T_fluid_initial = ones(n_segments,1).*homogeneous_temp_fluid_side;

%input heater power profile. 
%Assume that heater power is a constant 
p_total = 10e3; %[W]
p_profile = ones(n_segments,1).*p_total/n_segments;

%Temp dependent properties of dowtherm A. T in Kelvin
density_oil = @(T) 1078-(0.85*(T-273)); %[kg/m^3]
viscosity_oil = @(T) 0.130/((T-273)^1.072); %[kg/m s]
k_oil = @(T) 0.142 - (0.00016*(T-273)); %[W/m C]
c_p_oil = @(T) 1518 + 2.82.*(T-273); %[J/kg C]

%Temp dependent properties of steel. T in Kelvin
k_steel = @(T) 14.6 + (0.0127*(T-273));

%Temp dependent properties of insulation. T in Kelvin
k_insulation = @(T) 0.206 + (7.702e-4)*(T-273);


%Nusselt correlation for turbulent fluids
Nu = @(Re,Pr) 0.024*(Re^0.8)*(Pr^0.33);

%Pr definition
Re 


%create storer
T_heater{1} = T_heater_initial;
T_fluid{1} = T_fluid_initial;





%heat equation for resister heater
%rate of change = power added - heat loss to coolant + heat diffusion due
%to temp grad before and after 

%heat equation for fluid
%rate of change = heat gain from HS + heat flow due to mass flow rate

for t = 1:t_segments
    %calculate d2T/dx2 using central point difference for the heater side
    d2T_data = [ones(2,1).*T_heater{t}(1); T_heater{t} ;ones(2,1).*T_heater{t}(end)];
    for i = 3:n_segments+2
        %special consideration for the first 2 and last 2 values
        d2T_dx2(i-2) = (-d2T_data(i+2) + 16*d2T_data(i+1)-30*d2T_data(i)+16*d2T_data(i-1)-d2T_data(i-2))/(12*(x_step^2));
    end
    
    %Fluid equation: c_p(T(i))*x_step*pi*r^2*density*T(i)' = m_flow.*(T(i) -
    %T(i-1)
    
    %update Temperature profile at time t
    
end