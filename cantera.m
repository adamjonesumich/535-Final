clear
clc
% Currently set up for testing pyrun

CANTERA_test(20,-2)

%function rocket = cantera(rocket)
function output = CANTERA_test(input,input2)
% valid inputs: chamber_pressure, fuel_species, oxidizer_species,
% mass_fuel_ratio

    pycode = [
                "# Here is where the python code goes"
                "# With each line as a string list element"
                "# Like so"
                "output = 2*A + B"
                "print(output)"
                ];
    
    % Inputs are passed in using the A=___, B=___ expressions
    % The variable output is saved as the pyrun output since it is 
    % Designated below
    out = pyrun(pycode,'output',A=input,B=input2);
    
    %py_out = pyrun(pycode,'outVec',...)

    % Assign cantera outputs
    % rocket.chamber_temperature = py_out(1);
    % rocket.mixture_gamma = py_out(2);
    % rocket.mixture_molecular_weight = py_out(3);
end
