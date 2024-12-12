classdef Rocket
    %RocketParameters Container to pass variables between functions
    %   Detailed explanation goes here

    properties
        % -------------inputs--------------
        % mission parameters - nonadjustable
        total_delta_v % m/s

        % architecture dependent - basically nonadjustable
        number_stages
        fuel_species
        oxidizer_species
        average_ambient_pressure % Pa 
        % could change to be a vector/range; also could change to be 
        % calculated somehow, or be a product of thrust

        % adjustable
        mass_fuel_ratio
        area_throat % m2
        chamber_pressure % Pa
        % ----------------------------------
        
        % --------cantera outputs-----------
        chamber_temperature % K
        mixture_molecular_weight % g/mol
        mixture_gamma
        % ----------------------------------
        
        % --------quasi-1D outputs----------
        area_exit_ratio
        mass_flow_rate % kg/s
        exit_mach_number
        specific_impulse % s
        % ----------------------------------

        % ------multistaging outputs--------
        stage_structural_coefficients
        stage_mass_ratios
        stage_masses % kg
    end

    methods
        function obj = Rocket(total_delta_v,number_stages,fuel_species,oxidizer_species,mass_fuel_ratio,area_throat,average_ambient_pressure,chamber_pressure)
            %Rocket Construct an instance of this class
            %   Detailed explanation goes here
            obj.total_delta_v = total_delta_v;
            obj.number_stages = number_stages;
            obj.fuel_species = fuel_species;
            obj.oxidizer_species = oxidizer_species;
            obj.mass_fuel_ratio = mass_fuel_ratio;
            obj.area_throat = area_throat;
            obj.average_ambient_pressure = average_ambient_pressure;
            obj.chamber_pressure = chamber_pressure;

            obj.chamber_temperature = NaN;
            obj.mixture_molecular_weight = NaN;
            obj.mixture_gamma = NaN;
    
            obj.area_exit_ratio = NaN;
            obj.mass_flow_rate = NaN;
            obj.exit_mach_number = NaN;
            obj.specific_impulse = NaN;
            
            obj.stage_structural_coefficients = NaN;
            obj.stage_mass_ratios = NaN;
            obj.stage_masses = NaN;
        end

    end
end