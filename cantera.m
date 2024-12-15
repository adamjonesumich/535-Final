function rocket = cantera(rocket)
% valid inputs: chamber_pressure, fuel_species, oxidizer_species,
% mass_fuel_ratio

    pycode = [
    "import cantera as ct"
    "def solve_cantera(temp,pressure,fuel_ratio,fuel_species,oxidizer_species):"
        "gas = ct.Solution('gri30.yaml')"
        "initial_mass_fraction = fuel_species + ':1, ' + oxidizer_species + ':' + fuel_ratio"
        "gas.Y = initial_mass_fraction"
        "gas.TP = temp,pressure"
        "gas.equilibrate('HP')"
        "final_temp = gas.T"
        "gamma = gas.cp / gas.cv"
        "MW = gas.mean_molecular_weight"
        "return [final_temp, gamma, MW]"
    "output = solve_cantera(Tc,Pc,MFR,FS,OS)"
                ];
    
    % Inputs are passed in using the A=___, B=___ expressions
    % The variable 'output' is saved as the pyrun output since it is 
    % Designated below
    py_out = pyrun(pycode,'output',...
                    Tc  =   rocket.chamber_temperature,... % NEED 
                    Pc  =   rocket.chamber_pressure,...
                    MFR =   rocket.mass_fuel_ratio,...
                    FS  =   rocket.fuel_species,...
                    OS  =   rocket.oxidizer_species);

    % Assign cantera outputs
    rocket.chamber_temperature = py_out(1);
    rocket.mixture_gamma = py_out(2);
    rocket.mixture_molecular_weight = py_out(3);
end

