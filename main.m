% initialize

total_delta_v = [8000;9600]; % m/s
upper_stage_mass = 1000; % kg
number_launch_stages = 3; 
fuel_species = "CH4";
oxidizer_species = "O2";
mass_fuel_ratio = [4; 4; 4; 4];
area_throat = [0.05; 0.05; 0.05; 0.05]; % m2
highest_ambient_pressure = [101325; 10; 0; 0]; % Pa
lowest_ambient_pressure = [0; 0; 0; 0]; % Pa
chamber_pressure = [5; 5.5; 6; 6] * 1e6; % Pa

% struct to hold data between function calls
rocket = Rocket(total_delta_v,upper_stage_mass,number_launch_stages,fuel_species,oxidizer_species, ...
    mass_fuel_ratio,area_throat,lowest_ambient_pressure,highest_ambient_pressure,chamber_pressure);

% solve!

rocket = cantera(rocket);

rocket = quasi_1d(rocket);
rocket = nozzle_adjust(rocket);
rocket = multistage(rocket);

