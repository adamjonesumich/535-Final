function rocket = cantera(rocket)
% valid inputs: chamber_pressure, fuel_species, oxidizer_species,
% mass_fuel_ratio

    pycode = [
    "import cantera as ct"
    "def solve_cantera(temp,pressure,fuel_ratio,fuel_species,oxidizer_species):"
    "   gas = ct.Solution('gri30.yaml')"
    "   initial_mass_fraction = fuel_species + ':1, ' + oxidizer_species + ':' + str(fuel_ratio)"
    "   gas.Y = initial_mass_fraction"
    "   gas.TP = temp,pressure"
    "   gas.equilibrate('HP')"
    "   final_temp = gas.T"
    "   gamma = gas.cp / gas.cv"
    "   MW = gas.mean_molecular_weight"
    "   return [final_temp, gamma, MW]"
    "output = solve_cantera(T_init,Pc,MFR,FS,OS)"
                ];
    
    % Inputs are passed in using the A=___, B=___ expressions
    % The variable 'output' is saved as the pyrun output since it is 
    % Designated below
    for ii = 1:(rocket.number_launch_stages+1)
        py_out = double(pyrun(pycode,'output',...
                        T_init  =   300,... % K
                        Pc  =   rocket.effective_chamber_pressure(ii),...
                        MFR =   rocket.mass_fuel_ratio(ii),...
                        FS  =   rocket.fuel_species,...
                        OS  =   rocket.oxidizer_species));
    
        % Assign cantera outputs
        rocket.chamber_temperature(ii) = py_out(1);
        rocket.mixture_gamma(ii) = py_out(2);
        rocket.mixture_molecular_weight(ii) = py_out(3);
    end
end

