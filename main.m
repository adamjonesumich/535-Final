% initialize

total_delta_v = [8000,9600]; % m/s
upper_stage_mass = 5000; % kg
number_launch_stages = 3; 
%number_engines_per_stage = [1, 1, 1, 1];
number_engines_per_stage = [9, 1, 1, 1];
max_area_ratio = 200;
max_exit_area = 40; % TBD
fuel_species = "CH4";
oxidizer_species = "O2";
%area_throat = [0.5, 0.05, 0.01, 0.01]; % m2
area_throat = [0.05, 0.05, 0.01, 0.001]; % m2
lowest_ambient_pressure = [0, 0, 0, 0]; % Pa

% adjustables
mass_fuel_ratio = [3.8, 3.8, 3.7, 3.8];
highest_ambient_pressure = [101325, 10, 0, 0]; % Pa
ideal_chamber_pressure = [60, 55, 50, 50] * 1e6; % Pa
combustion_area_ratio = [2.5,2,2,1];
        

% struct to hold data between function calls
rocket = Rocket(total_delta_v,upper_stage_mass,number_launch_stages, ...
    number_engines_per_stage,max_area_ratio,max_exit_area,fuel_species, ...
    oxidizer_species,mass_fuel_ratio,area_throat,lowest_ambient_pressure, ...
    highest_ambient_pressure,ideal_chamber_pressure,combustion_area_ratio);

% solve!

rocket = cantera(rocket);
rocket = quasi_1d(rocket);
rocket = combustion_chamber(rocket); % TODO: need to iterate to ensure convergence of gamma and exit area ratio
rocket = quasi_1d(rocket);
rocket = nozzle_adjust(rocket);
rocket = multistage(rocket);

