function output = dT_dt_2D(T,T_inlet,p_density,D_hydraulic,r_inner,x_step,r_step,volume_heater,density_steel, A_HS,Vol_fluid, Vol_r_segment, inner_assembly_mass,mass_flow_fluid,n_segments,r_segments,A_r_normal,A_x_normal,S_r) 

%each row is each segment. columns are as follows: [T_fluid T_r1
%T_r2...T_rn]
T_fluid = T(:,1);
T_heater = T(:,2:end);

%%Temp dependent properties
%Temp dependent properties of dowtherm A. T in Kelvin
density_oil = @(T) 1078-(0.85.*(T-273)); %[kg/m^3]
viscosity_oil = @(T) 0.130./((T-273).^1.072); %Dynamic viscosity [kg/m s]
k_oil = @(T) 0.142 - (0.00016.*(T-273)); %Thermal conductivity [W/m C]
c_p_oil = @(T) 1518 + 2.82.*(T-273); %Specific heat capacity [J/kg C]

%Temp dependent properties of steel. T in Kelvin
c_p_steel = @(T)  450 + 0.28.*(T-273); %Specific heat capacity [J/kg C]
k_steel = @(T) 14.6 + (0.0127.*(T-273)); %Thermal conductivity [W/m C]

%Temp depedent properties of insulation. T in Kelvin
k_insulation = @(T) 0.206 + (7.702e-4)* 80; %Thermal conductivity [W/m C]

%% Other variables definition
h = @(Nu,k,D_hydraulic) Nu.*k./D_hydraulic; %[W/m^2 K]
velocity = @(mass_flow_fluid,density,D_hydraulic) (mass_flow_fluid./density)/(pi*(D_hydraulic^2)/4); %[m/s]


%% Dimensionless number definitions
Re = @(density, viscosity, velocity, length) density.*velocity.*length./viscosity;
Pr = @(heat_capacity, viscosity, thermal_conductivity) heat_capacity.*viscosity./thermal_conductivity;
Nu = @(Re,Pr) 0.024.*(Re.^0.8).*(Pr.^0.33); %Nu correlation for turbulent flow
%Nu = @(Re,Pr) 5.44 + 0.034*(Re^0.82)*(Pr^0);
Re_calc = Re(density_oil(T_fluid), viscosity_oil(T_fluid), velocity(mass_flow_fluid,density_oil(T_fluid),D_hydraulic), D_hydraulic);
Pr_calc = Pr(c_p_oil(T_fluid), viscosity_oil(T_fluid), k_oil(T_fluid));
Nu_calc = Nu(Re_calc,Pr_calc);
h_calc = h(Nu_calc,k_oil(T_fluid),D_hydraulic);

%% Calculate output
dT_dt_heater = zeros(n_segments,r_segments);
dT_dx_calc = [T_heater(1,:); T_heater; T_heater(end,:)]; %This portion is to calculate the dT/dx 

for n = 1:n_segments
      dT_dr_calc = [T_heater(n,1) T_heater(n,1:end-1)];
 
      %Heater term
      %calculate Q_r_inner contribution. There is no inner contribution in
      %inner-most section. 
      dT_dr_inner = T_heater(n,:)- dT_dr_calc;
      Q_r_inner = -k_steel(T_heater(n,:)).*dT_dr_inner.*[A_r_normal(1) A_r_normal(1:end-1)]./(r_step);
      %logic array. +ve is heat flowing in positive r direction. This is essentially the heat gain from the rth inner surface in the rth segment is equal to the heat loss in the r-1 segment. 
      Q_r_outer = [Q_r_inner(2:end).*-1 0]; %no heat loss through outer surface
      Q_x = -k_steel(T_heater(n,:)).*A_x_normal.*((-dT_dx_calc(n,:) + 2.*dT_dx_calc(n+1,:) - dT_dx_calc(n+2,:))./x_step);
      
      %Source term due to electrical heating
      %S_r
      
      %Heat exchange between Heater and Fluid
      Q_HS = [h_calc(n).*A_HS.*(T_heater(n,1)-T_fluid(n)) zeros(1,r_segments-1)]; 
      
      %Total heat = heat at inner ring, heat at outer ring, heat difference
      %at top and bottom ring, source  
      Q_total = Q_r_inner + Q_r_outer + S_r + Q_x - Q_HS;
      dT_dt_heater(n,:) = Q_total./(c_p_steel(T_heater(n,:)).*Vol_r_segment.*density_steel);
end

T_fluid_calc = [T_inlet; T_fluid(1:end-1)]; %This line is to calculate the thermal changes due to mass flow into CV
heat_cap_fluid =  Vol_fluid.*density_oil(T_fluid).*c_p_oil(T_fluid) + inner_assembly_mass*(c_p_steel(T_fluid));
dT_dt_fluid = (1./heat_cap_fluid).*((mass_flow_fluid.*c_p_oil(T_fluid).*(T_fluid_calc-T_fluid)) + h_calc.*A_HS.*(T_heater(:,1) - T_fluid));
output = [dT_dt_fluid dT_dt_heater];

end