% initialize

total_delta_v = 8000; % m/s
number_stages = 3; 
fuel_species = "CH4";
oxidizer_species = "O2";
mass_fuel_ratio = 4;
area_throat = 0.05; % m2
average_ambient_pressure = 0; % Pa
chamber_pressure = 5 * 1e6;% Pa

rocket = Rocket(total_delta_v,number_stages,fuel_species,oxidizer_species, ...
    mass_fuel_ratio,area_throat,average_ambient_pressure,chamber_pressure);

rocket = cantera(rocket);
% rocket.chamber_temperature = 1500;
% rocket.mixture_gamma = 1.2;
% rocket.mixture_molecular_weight = 16;
rocket = quasi_1d(rocket);
rocket = nozzle_adjust(rocket);
rocket = multistage(rocket);


function rocket = cantera(rocket)
    % valid inputs: chamber_pressure, fuel_species, oxidizer_species,
    % mass_fuel_ratio

    % expected outputs: chamber_temperature, mixture_gamma, 
    % mixture_molecular_weight
end

function rocket = nozzle_adjust(rocket)
    % valid inputs: all thermodynamic variables

    % expected output: adjusted specific_impulse

    % note: this is a combination of Gamba's 2D nozzle solver, the "nozzle
    % chopper", and the Isp adjuster script
end

function rocket = multistage(rocket)
    % valid inputs: all thermodynamic variables

    % expected output: all structural variables
end

