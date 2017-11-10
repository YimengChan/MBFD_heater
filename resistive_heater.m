%%Calculate 1D axial heater. Each segment acts as a homogeneous control
%%volume. The heater is divided into n_segments, and seperated into the
%%fluid portion and the resistive-heater portion
clear;clc;clf
n_segments = 1000;

%input resistive heater dimensions
total_x = 2; %[m]
x_step = total_x/n_segments;
x_profile = linspace(0,total_x,n_segments);
d2T_dx2 = ones(n_segments,1);
D_hydraulic = 6.6e-3;
r_inner = 0.0381; %[m]
r_outer = 0.04; %[m]
r_insulation_outer = 0.04; %[m]
A_ring = pi*(r_outer^2 - r_inner^2);%Area for the ring section for conductive heat transfer 
volume_heater = A_ring*x_step;
A_HS = 2*pi*r_inner*x_step;

%Inner perforated steel and twisted metal contributes to thermal inertia
inner_assembly_mass = 3.120/n_segments ;%[kg]
vol_fluid = pi*(r_inner^2) - (inner_assembly_mass/7700); %m3


%input fluid flow
mass_flow_fluid = 0.018; %[kg/s]
velocity = @(mass_flow_fluid,density,D_hydraulic) (mass_flow_fluid./density)/(pi*(D_hydraulic^2)/4);

%Time step
time_end = 1; %s
t_segments = 1000;
t_step = time_end/t_segments;
t_profile = linspace(0,time_end,t_segments);

%Store all temperatures in both portion in a matrix [resister_T fluid_T];
%input initial temperature profile in both portions
homogeneous_temp_heater_side = 323; %[K]
homogeneous_temp_fluid_side = 323; %[K]
T_heater_initial = ones(n_segments,1).*homogeneous_temp_heater_side;
T_inlet = 324;
T_fluid_initial = ones(n_segments,1).*homogeneous_temp_fluid_side;

%input heater power profile. 
%Assume that heater power is a constant 
p_total = 0; %[W]
p_profile = ones(n_segments,1).*p_total/n_segments;

%Temp dependent properties of dowtherm A. T in Kelvin
density_oil = @(T) 1078-(0.85.*(T-273)); %[kg/m^3]
viscosity_oil = @(T) 0.130./((T-273).^1.072); %Dynamic viscosity [kg/m s]
k_oil = @(T) 0.142 - (0.00016.*(T-273)); %Thermal conductivity [W/m C]
c_p_oil = @(T) 1518 + 2.82.*(T-273); %[J/kg C]

%Temp dependent properties of steel. T in Kelvin
c_p_steel = @(T)  450 + 0.28.*(T-273);
k_steel = @(T) 14.6 + (0.0127.*(T-273));

%assume k_insulation is a constant. T in Kelvin
k_insulation = @(T) 0.206 + (7.702e-4)* 80;


%Nusselt correlation for turbulent fluids
Nu = @(Re,Pr) 0.024.*(Re.^0.8).*(Pr.^0.33);
%Nu_calc = ones(n_segments,1);

%convective heat transfer
h = @(Nu,k,D_hydraulic) Nu.*k./D_hydraulic;
%h_calc = ones(n_segments,1);

%Re definition
Re = @(density, viscosity, velocity, length) density.*velocity.*length./viscosity;
%Re_calc = ones(n_segments,1);

%Pr definition
Pr = @(heat_capacity, viscosity, thermal_conductivity) heat_capacity.*viscosity./thermal_conductivity;

%Thermal resistance

R_HS = @(k) (log(r_outer/r_inner))./(2*pi*x_step.*k);
R_fluid = @(h) (2*pi*r_inner*x_step*h);
R_fluid_calc = ones(n_segments,1);


R_insul = @(k) (log(r_insulation_outer/r_outer))/(2*pi*k*x_step);
R_insul_calc = ones(n_segments,1);

%create storer
T_heater = T_heater_initial;
T_fluid = T_fluid_initial;
%T_heater_storer{1} = T_heater_initial;
%T_fluid_Storer{1} = T_fluid_initial;

%heat equation for resister heater
%rate of change = power added - heat loss to coolant + heat diffusion due
%to temp grad before and after 

%heat equation for fluid
%rate of change = heat gain from HS + heat flow due to mass flow rate
R_HS_calc = R_HS(k_steel(T_heater));
    Re_calc = Re(density_oil(T_fluid), viscosity_oil(T_fluid), velocity(mass_flow_fluid,density_oil(T_fluid),D_hydraulic), D_hydraulic);
    Pr_calc = Pr(c_p_oil(T_fluid), viscosity_oil(T_fluid), k_oil(T_fluid));
    Nu_calc = Nu(Re_calc,Pr_calc);
    h_calc = h(Nu_calc,k_oil(T_fluid),D_hydraulic);
heat_cap_fluid =  vol_fluid.*c_p_oil(T_fluid) + inner_assembly_mass*(c_p_steel(T_fluid));
T_fluid_calc = [T_inlet; T_fluid(1:end-1)];
    dT_fluid_dt = (1./heat_cap_fluid).*((mass_flow_fluid.*c_p_oil(T_fluid).*(T_fluid-T_fluid_calc)) + h_calc.*A_HS.*(T_heater - T_fluid));
    
    t_store_index = 0;
    
for t = 1:t_segments
    
    %Heater thermal equation: c_p(T(i))*m*dT/dt = Q_heater + Q_Loss to
    %coolant
     
    %calculate updated k values
    R_HS_calc = R_HS(k_steel(T_heater));
    Re_calc = Re(density_oil(T_fluid), viscosity_oil(T_fluid), velocity(mass_flow_fluid,density_oil(T_fluid),D_hydraulic), D_hydraulic);
    Pr_calc = Pr(c_p_oil(T_fluid), viscosity_oil(T_fluid), k_oil(T_fluid));
    Nu_calc = Nu(Re_calc,Pr_calc);
    h_calc = h(Nu_calc,k_oil(T_fluid),D_hydraulic);
    
    %U = h_calc
    
    %calculate d2T/dx2 using central point difference for the heater side
    d2T_data = [ones(2,1).*T_heater(1); T_heater ;ones(2,1).*T_heater(end)];
    for i = 3:n_segments+2
        %special consideration for the first 2 and last 2 values
        d2T_dx2(i-2) = (-d2T_data(i+2) + 16*d2T_data(i+1)-30*d2T_data(i)+16*d2T_data(i-1)-d2T_data(i-2))/(12*(x_step^2));
    end
    
    dT_heater_dt = (1./(volume_heater.*c_p_steel(T_heater))).*((volume_heater.*c_p_steel(T_heater).*k_oil(T_fluid).*d2T_dx2) + p_profile - h_calc.*A_HS.*(T_heater - T_fluid));
    heat_cap_fluid =  vol_fluid.*c_p_oil(T_fluid) + inner_assembly_mass*(c_p_steel(T_fluid));
    dT_fluid_dt = (1./heat_cap_fluid).*((mass_flow_fluid.*c_p_oil(T_fluid).*(T_fluid-T_fluid_calc)) + h_calc.*A_HS.*(T_heater - T_fluid));
    
    T_heater = T_heater + dT_heater_dt.*t_step;
    T_fluid = T_fluid + dT_fluid_dt.*t_step;
    
    if ~mod(t,200)
         t_store_index = t_store_index + 1;
         t_store{t_store_index} = [T_heater T_fluid];
    end
    
    %check the 2nd differiantial 
    %update Temperature profile at time t
    
end