% initialize

total_delta_v = [8000;9600]; % m/s
upper_stage_mass = 1000; % kg
number_launch_stages = 3; 
fuel_species = "CH4";
oxidizer_species = "O2";
mass_fuel_ratio = [4; 4; 4; 4];
area_throat = [0.05; 0.05; 0.05; 0.05]; % m2
average_ambient_pressure = [101325; 100; 0; 0]; % Pa
chamber_pressure = [5; 5.5; 6; 6] * 1e6; % Pa

% struct to hold data between function calls
rocket = Rocket(total_delta_v,upper_stage_mass,number_launch_stages,fuel_species,oxidizer_species, ...
    mass_fuel_ratio,area_throat,average_ambient_pressure,chamber_pressure);

% solve!

rocket = cantera(rocket);

    % for testing purposes: fake cantera!
    % rocket.chamber_temperature = [1500; 1600; 1700; 1700];
    % rocket.mixture_gamma = [1.2; 1.15; 1.1; 1.1];
    % rocket.mixture_molecular_weight = [18; 17; 16; 16];

rocket = quasi_1d(rocket);
rocket = nozzle_adjust(rocket);
rocket = multistage(rocket);

function rocket = nozzle_adjust(rocket)
    % valid inputs: all thermodynamic variables

    % expected output: adjusted specific_impulse

    % note: this is a combination of Gamba's 2D nozzle solver, the "nozzle
    % chopper", and the Isp adjuster script
end

