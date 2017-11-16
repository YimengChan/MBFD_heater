function output = dT_dt(T,T_inlet,p_profile,D_hydraulic,r_inner,x_step,volume_heater,density_steel, A_HS,vol_fluid, inner_assembly_mass,mass_flow_fluid,n_segments,A_insulation, U_insulation,T_air) 

T_heater = T(:,1);
T_fluid = T(:,2);

%%Temp dependent properties
%Temp dependent properties of dowtherm A. T in Kelvin
density_oil = @(T) 1078-(0.85.*(T-273)); %[kg/m^3]
viscosity_oil = @(T) 0.130./((T-273).^1.072); %Dynamic viscosity [kg/m s]
k_oil = @(T) 0.142 - (0.00016.*(T-273)); %Thermal conductivity [W/m C]
c_p_oil = @(T) 1518 + 2.82.*(T-273); %Specific heat capacity [J/kg C]

%Temp dependent properties of steel. T in Kelvin
c_p_steel = @(T)  450 + 0.28.*(T-273); %Specific heat capacity [J/kg C]
k_steel = @(T) 14.6 + (0.0127.*(T-273)); %Thermal conductivity [W/m C]

%% Other variables definition
h = @(Nu,k,D_hydraulic) Nu.*k./D_hydraulic; %[W/m^2 K]
velocity = @(mass_flow_fluid,density,D_hydraulic) (mass_flow_fluid./density)/(pi*(D_hydraulic^2)/4); %[m/s]

%% Using central difference method to calculate d2T/dt2 and central point method requires i-2 and i+2 values. 
d2T_dx2 = ones(n_segments,1);


% d2T_data = [ones(2,1).*T_heater(1); T_heater ;ones(2,1).*T_heater(end)];
 %   for i = 3:n_segments+2
  %      
    %    d2T_dx2(i-2) = (-d2T_data(i+2) + 16*d2T_data(i+1)-30*d2T_data(i)+16*d2T_data(i-1)-d2T_data(i-2))/(12*(x_step^2));
   % end

    %try the other method to see if theres a difference
    d2T_data = [T_heater(1);T_heater;T_heater(end)];
    for i = 2:n_segments+1;
    d2T_dx2 (i-1) = (-d2T_data(i+1)+2.*d2T_data(i)-d2T_data(i-1))/(x_step^2);
    end
    
%% Convective resistance calculation
Re = @(density, viscosity, velocity, length) density.*velocity.*length./viscosity;
Pr = @(heat_capacity, viscosity, thermal_conductivity) heat_capacity.*viscosity./thermal_conductivity;
Nu = @(Re,Pr) 0.024.*(Re.^0.8).*(Pr.^0.33); %Nu correlation for turbulent flow
%Nu = @(Re,Pr) 5.44 + 0.034.*(Re.^0.82).*(Pr.^0);

Re_calc = Re(density_oil(T_fluid), viscosity_oil(T_fluid), velocity(mass_flow_fluid,density_oil(T_fluid),D_hydraulic), D_hydraulic);
Pr_calc = Pr(c_p_oil(T_fluid), viscosity_oil(T_fluid), k_oil(T_fluid));
Nu_calc = Nu(Re_calc,Pr_calc);
h_convec_calc = h(Nu_calc,k_oil(T_fluid),D_hydraulic);

%% Radiative resistance calculation 
emissitivity_steel = 0.80;
h_radiative_calc = emissitivity_steel.*(5.67e-8).*((T_heater).^2 + T_air^2).*(T_heater + T_air);

%% Calculate output
%rate of accumulation of energy within heater = diffusion within heater in
%x-direction + power source + convec heat loss to fluid + conduc heat loss to insulation-
%iar + radiative heat loss
dT_dt_heater = (1./(density_steel.*volume_heater.*c_p_steel(T_heater))).*((-volume_heater.*k_steel(T_heater).*d2T_dx2) + p_profile - h_convec_calc.*A_HS.*(T_heater - T_fluid) - (U_insulation + h_radiative_calc).*A_insulation.*(T_heater - T_air));
T_fluid_calc = [T_inlet; T_fluid(1:end-1)]; %This line is to calculate the thermal changes due to mass flow into CV
heat_cap_fluid =  vol_fluid.*density_oil(T_fluid).*c_p_oil(T_fluid) + inner_assembly_mass*(c_p_steel(T_fluid));
dT_dt_fluid = (1./heat_cap_fluid).*((mass_flow_fluid.*c_p_oil(T_fluid).*(T_fluid_calc-T_fluid)) + h_convec_calc.*A_HS.*(T_heater - T_fluid));
output = [dT_dt_heater dT_dt_fluid];

end