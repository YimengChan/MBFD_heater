%%Calculate 1D axial heater. Each segment acts as a homogeneous control
%%volume. The heater is divided into n_segments, and seperated into the
%%fluid portion and the resistive-heater portion
clear;clc;clf
n_segments = 10;

%input resistive heater dimensions
total_x = 2; %[m] Height of heater
x_step = total_x/n_segments;
x_profile = linspace(0,total_x,n_segments);
d2T_dx2 = ones(n_segments,1);
D_hydraulic = 6.6e-3;
r_inner = 0.0381; %[m]
r_outer = 0.04; %[m]
r_insulation_outer = 0.04; %[m]
A_ring = pi*(r_outer^2 - r_inner^2);%Area for the ring section for conductive heat transfer 
volume_heater = A_ring*x_step;
A_HS = 2*pi*r_inner*x_step; %[m^2] Surface area of contact of Heater and fluid

density_steel = 8030; % treated as constnat [kg/m3]

%Inner perforated steel and twisted metal contributes to thermal inertia
inner_assembly_mass = 3.120/n_segments ;%[kg]
vol_fluid = pi*(r_inner^2) - (inner_assembly_mass/7700); %m3 Difference between inner cylinder vol and the vol of the inner steel assembly


%input fluid flow
mass_flow_fluid = 0.018; %[kg/s]

%Time step
time_end = 5000; %s
t_segments = time_end; %num of time segments
t_step = time_end/t_segments; %s
num_of_snaps = 5; %This variable will tell the program how many snapshots will be taken
t_store_index = 0; %This variable will increase every time a snapshot is taken at 
t_profile = linspace(0,time_end,t_segments);
run_to_end = 0;

%input initial temperature profile in both portions
initial_homogeneous_temp_heater_side = 273+250; %[K]
initial_homogeneous_temp_fluid_side = 273+80; %[K]
T_heater_initial = ones(n_segments,1).*initial_homogeneous_temp_heater_side;
T_fluid_initial = ones(n_segments,1).*initial_homogeneous_temp_fluid_side;
T = [T_heater_initial T_fluid_initial];


%input heater power profile. 
%Assume that heater power is a constant 
p_total = 10000; %[W]
p_profile = ones(n_segments,1).*p_total/n_segments;

%input inlet temperature of fluid. Assume that this is a constant
T_inlet = 273+80; %[K] 

%calculate steady state values
T_steady = lsqnonlin(@(T) dT_dt(T,T_inlet,p_profile,D_hydraulic,r_inner,x_step,volume_heater,density_steel, A_HS,vol_fluid, inner_assembly_mass,mass_flow_fluid,n_segments),T);

for t = 1:t_segments
    
    
    %implement an fsolve
    dT_dt_calc = dT_dt(T,T_inlet,p_profile,D_hydraulic,r_inner,x_step,volume_heater,density_steel, A_HS,vol_fluid, inner_assembly_mass,mass_flow_fluid,n_segments);
    
    
    %Step forward in time Euler method
    T = T + dT_dt_calc.*t_step;
    
    %Taking a snapshot of T
    if ~mod(t,(t_segments/num_of_snaps))  
         t_store_index = t_store_index + 1;
         T_store{t_store_index} = T;
    end
    
    %Check if steady state is reached. 
    if run_to_end == 0
        if all((T_steady - T)./T_steady < 0.001)
        time_to_steady = t;
        prompt = sprintf('Steady state has been reached at %.1f seconds. Do you want to stop the simulation? Y/N ',time_to_steady);
        str = input(prompt,'s');
            if isempty(str)
            str = 'Y';
            end
            if strcmp(str,'Y')
            break
            else
            run_to_end = 1;
            end
        end
    end
 
end

hold on
for i = 1:t_store_index
plot(x_profile, T_store{i}(:,2))
end





%Thermal resistance
R_HS = @(k) (log(r_outer/r_inner))./(2*pi*x_step.*k); %not used
R_fluid = @(h) (2*pi*r_inner*x_step*h);
R_fluid_calc = ones(n_segments,1);
R_insul = @(k) (log(r_insulation_outer/r_outer))/(2*pi*k*x_step); %not used
