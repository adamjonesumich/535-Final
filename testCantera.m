clear
clc
% Just tests the cantera() function to make sure it's running correctly

% initialize

rocket.total_delta_v = [8000;9600]; % m/s
rocket.upper_stage_mass = 1000; % kg
rocket.number_launch_stages = 3; 
rocket.fuel_species = "CH4";
rocket.oxidizer_species = "O2";
rocket.mass_fuel_ratio = [4; 4; 4; 4];
rocket.area_throat = [0.05; 0.05; 0.05; 0.05]; % m2
rocket.average_ambient_pressure = [101325; 100; 0; 0]; % Pa
rocket.chamber_pressure = [5; 5.5; 6; 6] * 1e6; % Pa

% Pc  =   rocket.chamber_pressure,...
% MFR =   rocket.mass_fuel_ratio,...
% FS  =   rocket.fuel_species,...
% OS  =   rocket.oxidizer_species);

rocket = cantera(rocket)
