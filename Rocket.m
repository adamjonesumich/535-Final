classdef Rocket
    %RocketParameters Container to pass variables between functions
    %   Detailed explanation goes here

    properties
        % -------------inputs--------------
        % mission parameters - nonadjustable
        total_delta_v % m/s
        upper_stage_mass; % kg

        % architecture dependent - basically nonadjustable
        number_launch_stages 
        number_engines_per_stage
        max_area_ratio
        max_exit_area
        fuel_species
        oxidizer_species
        highest_ambient_pressure % Pa
        combustion_area_ratio

        % adjustable
        mass_fuel_ratio
        area_throat % m2
        ideal_chamber_pressure % Pa
        effective_chamber_pressure % Pa
        % ----------------------------------
        
        % combustion
        combustion_chamber_mass
        nozzle_mass
        
        % --------cantera outputs-----------
        chamber_temperature % K
        mixture_molecular_weight % g/mol
        mixture_gamma
        % ----------------------------------
        
        % --------quasi-1D outputs----------
        area_exit_ratio
        mass_flow_rate % kg/s
        exit_mach_number
        specific_impulse_vacuum % s
        thrust_vacuum; % N
        % ----------------------------------
        
        % ----------nozzle adjust-----------
        correction_factor
        specific_impulse_sea_level % s
        thrust_sea_level % N
        % ----------------------------------
        
        % ------multistaging outputs--------
        stage_payload_ratios
        stage_mass_ratios
        stage_masses % kg
        stage_delta_vs % m/s
        final_thrust % N
        initial_net_gs
        % ----------------------------------
    end

    methods
        function obj = Rocket(total_delta_v,upper_stage_mass,number_launch_stages, ...
                number_engines_per_stage,max_area_ratio,max_exit_area,fuel_species, ...
                oxidizer_species,mass_fuel_ratio,area_throat, ...
                highest_ambient_pressure,chamber_pressure,combustion_area_ratio)
            %Rocket Construct an instance of this class
            %   Detailed explanation goes here
            obj.total_delta_v = total_delta_v;
            obj.upper_stage_mass = upper_stage_mass;
            obj.number_launch_stages = number_launch_stages;
            obj.number_engines_per_stage = number_engines_per_stage;
            obj.max_area_ratio = max_area_ratio;
            obj.max_exit_area = max_exit_area;
            obj.fuel_species = fuel_species;
            obj.oxidizer_species = oxidizer_species;
            obj.mass_fuel_ratio = mass_fuel_ratio;
            obj.area_throat = area_throat;
            obj.highest_ambient_pressure = highest_ambient_pressure;
            obj.ideal_chamber_pressure = chamber_pressure;
            obj.effective_chamber_pressure = chamber_pressure;
            obj.combustion_area_ratio = combustion_area_ratio;

            obj.chamber_temperature = NaN;
            obj.mixture_molecular_weight = NaN;
            obj.mixture_gamma = NaN;
            
            obj.combustion_chamber_mass = NaN;
            obj.nozzle_mass = NaN;
    
            obj.area_exit_ratio = NaN;
            obj.mass_flow_rate = NaN;
            obj.exit_mach_number = NaN;
            obj.specific_impulse_vacuum = NaN;
            obj.thrust_vacuum = NaN;

            obj.correction_factor = NaN;
            obj.specific_impulse_sea_level = NaN;
            obj.thrust_sea_level = NaN;
            
            obj.stage_payload_ratios = NaN;
            obj.stage_mass_ratios = NaN;
            obj.stage_masses = NaN;
            obj.stage_delta_vs = NaN;
            obj.final_thrust = NaN;
            obj.initial_net_gs = NaN;
        end

    end
end