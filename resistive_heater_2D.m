%%Calculate 1D axial heater. Each segment acts as a homogeneous control
%%volume. The heater is divided into n_segments, and seperated into the
%%fluid portion and the resistive-heater portion. Heater wall is assumed to
%%be radially divided 

clear;clc;clf
n_segments = 10;
r_segments = 5; 
%% input resistive heater dimensions
x_total = 1.924; %[m] Height of heater
x_step = x_total/n_segments;
x_profile = linspace(0,x_total,n_segments);
d2T_dx2 = ones(n_segments,1);
%D_hydraulic = 6.6e-3; %[m]
D_hydraulic = 2.725e-2;
r_inner = 0.0381/2; %[m]
r_outer = 0.04/2; %[m]
r_insulation_thickness = 0.05; %[m]

r_step = (r_outer - r_inner)/r_segments;
r_profile = linspace(r_inner,r_outer,r_segments+1);
A_ring = pi*(r_outer^2 - r_inner^2);%Area for the ring section for conductive heat transfer 
Vol_heater = A_ring*x_total;
A_HS = 2*pi*r_inner*x_step; %[m^2] Surface area of contact of Heater and fluid
A_r_normal = 2.*pi.*r_profile(2:end)*x_step; %Surface area between inner and outer rings normal is in radial direction
A_x_normal =  (r_profile(2:end).^2 - r_profile(1:end-1).^2).*pi; %Surface area between top and bottom rings. Normal is in x direction
Vol_r_segment = A_x_normal.*x_step; %Volume of a segment in r direction. Vol. is constant in x direction


density_steel = 8030; % treated as constnat [kg/m3]

dT_dr = zeros(r_segments,n_segments);

%% input fluid flow
mass_flow_fluid = 0.18; %[kg/s]

%Time step
time_end = 1000; %s
t_segments = time_end*100; %num of time segments
t_step = time_end/t_segments; %s
num_of_snaps = 5; %This variable will tell the program how many snapshots will be taken
t_store_index = 0; %This variable will increase every time a snapshot is taken at 
t_profile = linspace(0,time_end,t_segments);
run_to_end = 0;

%input initial temperature profile in both portions
initial_homogeneous_temp_heater_side = 273+80; %[K]
initial_homogeneous_temp_fluid_side = 273+80; %[K]
T_heater_initial = ones(n_segments,r_segments).*initial_homogeneous_temp_heater_side;
T_fluid_initial = ones(n_segments,1).*initial_homogeneous_temp_fluid_side;
T = [T_fluid_initial T_heater_initial];

%input heater power profile. 
%Assume that heater power is a constant 
p_total = 9000*0.8; %[W]
p_density = p_total/Vol_heater;
S_r = p_density.*Vol_r_segment;

%input inlet temperature of fluid. Assume that this is a constant
T_inlet = 273+80; %[K] 

%Inner perforated steel and twisted metal contributes to thermal inertia
inner_assembly_mass = 3.120/n_segments ;%[kg]
Vol_fluid = pi*(r_inner^2)*x_step - (inner_assembly_mass/density_steel); %m3 Difference between inner cylinder vol and the vol of the inner steel assembly

T_steady = lsqnonlin(@(T) dT_dt_2D(T,T_inlet,p_density,D_hydraulic,r_inner,x_step,r_step,Vol_heater,density_steel, A_HS,Vol_fluid, Vol_r_segment, inner_assembly_mass,mass_flow_fluid,n_segments,r_segments,A_r_normal,A_x_normal,S_r),T) ;


for t = 1:t_segments
    
    %implement an fsolve
    dT_dt_calc = dT_dt_2D(T,T_inlet,p_density,D_hydraulic,r_inner,x_step,r_step,Vol_heater,density_steel, A_HS,Vol_fluid, Vol_r_segment, inner_assembly_mass,mass_flow_fluid,n_segments,r_segments,A_r_normal,A_x_normal,S_r);
    
    
    %Step forward in time Euler method
    T = T + dT_dt_calc.*t_step
    
    
    %Taking a snapshot of T
    if ~mod(t,(t_segments/num_of_snaps))  
         t_store_index = t_store_index + 1;
         T_store{t_store_index} = T;
    end
    
 
end




%Thermal resistance
R_HS = @(k) (log(r_outer/r_inner))./(2*pi*x_step.*k); %not used
R_fluid = @(h) (2*pi*r_inner*x_step*h);
R_fluid_calc = ones(n_segments,1);
R_insul = @(k) (log(r_insulation_outer/r_outer))/(2*pi*k*x_step); %not used
