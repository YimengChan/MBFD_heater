%%Calculate 1D axial heater. Each segment acts as a homogeneous control
%%volume. The heater is divided into n_segments, and seperated into the
%%fluid portion and the resistive-heater portion
clear;clc;clf
n_segments = 10;

%% input resistive heater dimensions
total_x = 1.924; %[m] Height of heater
x_step = total_x/n_segments;
x_profile = linspace(0,total_x,n_segments);
d2T_dx2 = ones(n_segments,1);
%D_hydraulic = 6.6e-3; %[m]
D_hydraulic = 2.725e-2;
r_inner = 0.0381/2; %[m]
r_outer = 0.04/2; %[m]
r_insulation_thickness = 0.05; %[m]
A_ring = pi*(r_outer^2 - r_inner^2);%Area for the ring section for conductive heat transfer 
volume_heater = A_ring*x_step; %[m^2] %Volume of heater portion in each n segment 
A_HS = 2*pi*r_inner*x_step; %[m^2] Surface area of contact of Heater and fluid of n segment
A_insulation = 2*pi*r_outer*x_step; %[m^2] Surface area of contact of Heater and insulation of n segment
density_steel = 8030; % treated as constnat [kg/m3]

%Inner perforated steel and twisted metal contributes to thermal inertia
inner_assembly_mass = 3.120/n_segments ;%[kg]
vol_fluid = pi*(r_inner^2)*x_step - (inner_assembly_mass/density_steel); %m3 Difference between inner cylinder vol and the vol of the inner steel assembly

%% input fluid flow
mass_flow_fluid = 0.18; %[kg/s]

%Time step
time_end = 150; %s
t_segments = time_end*1000; %num of time segments
t_step = time_end/t_segments; %s
num_of_snaps = 5; %This variable will tell the program how many snapshots will be taken
t_store_index = 0; %This variable will increase every time a snapshot is taken at 
t_profile = linspace(0,time_end,t_segments);
run_to_end = 0;

%input initial temperature profile in both portions and air temperature

T_air = 273 + 18; %[K]
initial_homogeneous_temp_heater_side = 273+80; %[K]
initial_homogeneous_temp_fluid_side = 273+80; %[K]
T_heater_initial = ones(n_segments,1).*initial_homogeneous_temp_heater_side;
T_fluid_initial = ones(n_segments,1).*initial_homogeneous_temp_fluid_side;
T = [T_heater_initial T_fluid_initial];

%Temp indepedent properties of insulation. T in Kelvin
k_insulation = 0.206 + (7.702e-4)* 80; %Thermal conductivity [W/m C]
U_insulation = ((log((r_insulation_thickness+r_outer)/r_outer))/(k_insulation))^-1;

%input heater power profile. 
%Assume that heater power is a constant 
p_total = 9000; %[W]
p_profile = ones(n_segments,1).*p_total/n_segments;

%input inlet temperature of fluid. Assume that this is a constant
T_inlet = 273+80; %[K] 

%calculate steady state values. 
T_steady = lsqnonlin(@(T) dT_dt(T,T_inlet,p_profile,D_hydraulic,r_inner,x_step,volume_heater,density_steel, A_HS,vol_fluid, inner_assembly_mass,mass_flow_fluid,n_segments,A_insulation, U_insulation,T_air),T);
A = dT_dt(T_steady,T_inlet,p_profile,D_hydraulic,r_inner,x_step,volume_heater,density_steel, A_HS,vol_fluid, inner_assembly_mass,mass_flow_fluid,n_segments,A_insulation, U_insulation,T_air);

for t = 1:t_segments
    
    
    %implement an fsolve
    dT_dt_calc = dT_dt(T,T_inlet,p_profile,D_hydraulic,r_inner,x_step,volume_heater,density_steel, A_HS,vol_fluid, inner_assembly_mass,mass_flow_fluid,n_segments,A_insul, U_insulation,T_air);
    
    
    %Step forward in time Euler method
    T = T + dT_dt_calc.*t_step;
    
    %Taking a snapshot of T
    if ~mod(t,(t_segments/num_of_snaps))  
         t_store_index = t_store_index + 1;
         T_store{t_store_index} = T;
    end
    
    %Check if steady state is reached. 
   % if run_to_end == 0
      %  if all((abs(T_steady - T))./T_steady < 0.001)
      %  time_to_steady = t*t_step;
      %  prompt = sprintf('Steady state has been reached at %.1f seconds. Do you want to stop the simulation? Y/N ',time_to_steady);
      %  str = input(prompt,'s');
       %     if isempty(str)
       %     str = 'Y';
       %     end
        %    if strcmp(str,'Y')
       %     break
       %     else
       %     run_to_end = 1;
       %    end
      %  end
 %  end
 
end

plot(x_profile,T(:,1),x_profile,T(:,2))

%figure
hold on
for i = 1:t_store_index
plot(x_profile, T_store{i}(:,1))
end

%figure 
hold on
for j = 1:t_store_index
%plot(x_profile, T_store{j}(:,2))
end




%Thermal resistance
R_HS = @(k) (log(r_outer/r_inner))./(2*pi*x_step.*k); %not used
R_fluid = @(h) (2*pi*r_inner*x_step*h);
R_fluid_calc = ones(n_segments,1);
R_insul = @(k) (log(r_insulation_outer/r_outer))/(2*pi*k*x_step); %not used
