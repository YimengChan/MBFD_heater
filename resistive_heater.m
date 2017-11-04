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
    
    %update Temperature profile at time t
    
end