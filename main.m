% initialize

total_delta_v = [8000,9600]; % m/s
upper_stage_mass = 5000; % kg
number_launch_stages = 3; 
max_area_ratio = [200,200,200,200];
max_exit_area = 40; % TBD
fuel_species = "CH4";
oxidizer_species = "O2";
highest_ambient_pressure = [101325, 100, 0, 0]; % Pa

% adjustables
number_engines_per_stage = [27, 11, 5, 1];
area_throat_launch = 0.01; % m2
area_throat_upper = 0.001; % m2
mass_fuel_ratio_launch = 2.95; 
mass_fuel_ratio_upper = 2.96;
ideal_chamber_pressure_launch = 55e6; % Pa
ideal_chamber_pressure_upper = 55e6; % Pa
combustion_area_ratio_launch = 1.5;
combustion_area_ratio_upper = 1;

area_throat = [area_throat_launch*ones(1,3),area_throat_upper];
mass_fuel_ratio = [mass_fuel_ratio_launch*ones(1,3),mass_fuel_ratio_upper];
ideal_chamber_pressure = [ideal_chamber_pressure_launch*ones(1,3),ideal_chamber_pressure_upper];
combustion_area_ratio = [combustion_area_ratio_launch*ones(1,3),combustion_area_ratio_upper]; 

% struct to hold data between function calls
rocket = Rocket(total_delta_v,upper_stage_mass,number_launch_stages, ...
    number_engines_per_stage,max_area_ratio,max_exit_area,fuel_species, ...
    oxidizer_species,mass_fuel_ratio,area_throat, ...
    highest_ambient_pressure,ideal_chamber_pressure,combustion_area_ratio);

% solve!

rocket = cantera(rocket);
rocket = combustion_chamber(rocket); % TODO: need to iterate to ensure convergence of gamma and exit area ratio
rocket = quasi_1d(rocket,false);
rocket = nozzle_adjust(rocket);
rocket = quasi_1d(rocket,true);
rocket = multistage(rocket);

% rocket.correction_factor = rocket.thrust_sea_level ./ (rocket.stage_masses * 9.81) - 1;
% delta_mass = rocket.stage_masses(1) - prev_mass;
% rocket.number_launch_stages = delta_mass;
% prev_mass = rocket.stage_masses(1);
