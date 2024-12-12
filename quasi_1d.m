function rocket = quasi_1d(rocket)
    % valid inputs: chamber_temperature, chamber_pressure, 
    % average_ambient_pressure, mixture_gamma, mixture_molecular_weight,
    % area_throat

    % expected outputs: area_exit_ratio, mass_flow_rate, exit_mach_number,
    % preliminary specific_impulse
    
    % TODO: combustor static conditions, plot enhancements,
    % problem 2: combustor properties
    
    % native sim variables
    num_points = 2; % good? 
    T0 = rocket.chamber_temperature; % K
    p0 = rocket.chamber_pressure; % Pa
    p_a = rocket.average_ambient_pressure; % Pa
    MW = rocket.mixture_molecular_weight; % g/mol
    MW = 1/1000 * MW; % kg/mol
    y = rocket.mixture_gamma; % gamma
    A_t = rocket.area_throat; % m2
    
    exit_area_ratio = 15; %
    combustor_area_ratio = 1; % TODO: what should this be?

    params = Quasi1DNozzleParameters(T0,p0,p_a,MW,y,A_t,exit_area_ratio,combustor_area_ratio,num_points);

    % find area ratio based on perfect expansion pressure ratio
    M_e = find_mach_from_pressure_ratio(params.p0/p_a,y);
    params.exit_area_ratio = calculate_area_ratio_from_Mach(M_e,y);
    params = update_area(params);

    results = solve_nozzle(params);

    rocket.area_exit_ratio = params.exit_area_ratio;
    rocket.mass_flow_rate = results.mdot;
    rocket.exit_mach_number = results.M_x(end);
    rocket.specific_impulse = results.Isp;
end


function results = solve_nozzle(params)
    
    % unpack
    T0 = params.T0;
    p0 = params.p0;
    p_a = params.p_a;
    MW = params.MW;
    y = params.y;
    A_t = params.A_t;
    exit_area_ratio = params.exit_area_ratio;
    x_l_c = params.x_l_c;
    x_l_d = params.x_l_d;
    area_ratio_c = params.area_ratio_c;
    area_ratio_d = params.area_ratio_d;
    
    R_u = 8.314; % J/molK
    R = R_u/MW;
    
    % expansion type and shock location
    exit_pressure_ratio = p_a/p0;
    
    M_e = Mach_from_area(exit_area_ratio,y,true);
    upper_bound_pressure_ratio = 1/stagnation_pressure_ratio(M_e,y); % this will give exactly sonic at throat, subsonic exit
    
    M_e = Mach_from_area(exit_area_ratio,y,false);
    lower_bound_pressure_ratio = 1/stagnation_pressure_ratio(M_e,y); % this will give perfect expansion
    
    M_e_2 = ns_jump_M2(M_e,y);
    p02_p01 = ns_jump_pressure(M_e,y);
    p02 = p0 * p02_p01;
    p2 = p02/stagnation_pressure_ratio(M_e_2,y);
    lower_bound_pressure_ratio_with_ns_at_exit = p2/p0;
    
    is_diverging_subsonic = false;
    normal_shock_exists = false;
    shock_location_index = -1;
    is_valid = true;

    tol = 5e-5;
    if abs(exit_pressure_ratio-lower_bound_pressure_ratio) < tol
        disp("Perfectly expanded flow.")
    elseif lower_bound_pressure_ratio > exit_pressure_ratio
        disp("Underexpanded flow.");
    elseif lower_bound_pressure_ratio_with_ns_at_exit > exit_pressure_ratio
        disp("Overexpanded flow with oblique shocks at exit.");
    elseif abs(exit_pressure_ratio-lower_bound_pressure_ratio_with_ns_at_exit) < 1e-5
        disp("Overexpanded flow with normal shock at exit.");
        is_valid = false;
        shock_location_index = length(x_l_d);
        normal_shock_exists = true;
    elseif upper_bound_pressure_ratio > exit_pressure_ratio
        disp("Overexpanded flow with normal shock within nozzle.");
        shock_location_index = find_normal_shock(area_ratio_d,y,exit_pressure_ratio,exit_area_ratio);
        normal_shock_exists = true;
        is_valid = false;
    else
        disp("Insufficient pressure differential to generate supersonic flow.");
        is_diverging_subsonic = true;
        is_valid = false;
    end
    
    % calculate distributions

    % independent variable
    x_l = [x_l_c,x_l_d];
    
    % calculate Mach distribution based on expansion type
    M_x_c = Mach_from_area(area_ratio_c,y,true);
    p_r_x_c = 1./stagnation_pressure_ratio(M_x_c,y);
    if(normal_shock_exists)
        M1 = Mach_from_area(area_ratio_d(shock_location_index),y,false);
        M2 = ns_jump_M2(M1,y);
        Astar1_Astar2 = calculate_area_ratio_from_Mach(M2,y)/calculate_area_ratio_from_Mach(M1,y);
    
        area_ratio_d1 = area_ratio_d(1:shock_location_index-1);
        area_ratio_d2 = area_ratio_d(shock_location_index:end) * Astar1_Astar2;
    
        M_x_d1 = Mach_from_area(area_ratio_d1,y,false);
        M_x_d2 = Mach_from_area(area_ratio_d2,y,true);
        M_x_d = [M_x_d1,M_x_d2];

        % finding pressure ratio post shock

        p02_p01 = ns_jump_pressure(M1,y);
        p_r_x_d1 = 1./stagnation_pressure_ratio(M_x_d1,y);
        p_r_x_d2 = 1./stagnation_pressure_ratio(M_x_d2,y) * p02_p01;
        p_r_x_d = [p_r_x_d1,p_r_x_d2];


    else
        M_x_d = Mach_from_area(area_ratio_d,y,is_diverging_subsonic);
        % calculate stagnation ratios
        p_r_x_d = 1./stagnation_pressure_ratio(M_x_d,y);
    end
    
    M_x = [M_x_c,M_x_d];
    p_ratio_x = [p_r_x_c,p_r_x_d];
    
    
    
    
    T_ratio_x = 1./stagnation_temperature_ratio(M_x,y);

    % calculate true temperature distribution
    T_x = T_ratio_x * T0;
    
    % calculate speed of sound
    a_x = sqrt(y*R*T_x);
    
    % calculate velocity
    u_x = M_x .* a_x;
    
    % calculate thrust
    u_e = u_x(end);
    A_e = A_t * exit_area_ratio;
    
    % % need sonic density
    % rho_0 = p0/(R*T0);
    % rho_s = rho_0 / stagnation_density_ratio(M_x(num_points),y);
    % 
    % % calculate mass flow rate
    % mdot = rho_s * A_t * u_x(num_points);

    p0e = p0;
    if normal_shock_exists
        p0e_p01 = ns_jump_pressure(M1,y);
        p0e = p0 * p0e_p01;
    end
    
    p_e = p0e / stagnation_pressure_ratio(M_x(end),y);
    rho_e = p_e / (R*T_x(end));
    
    % calculate mass flow rate
    mdot = rho_e * A_e * u_e;
    
    % thrust equation
    thrust = mdot * u_e + (p_e - p_a) * A_e;
    
    g0 = 9.81;
    
    Isp = thrust / (mdot*g0);
    
    coeff_thrust = thrust / (A_t * p0);

    
    results = Quasi1DNozzleResults(x_l,p_ratio_x,T_ratio_x,M_x,u_x,normal_shock_exists,...
        shock_location_index,p0e,mdot,thrust,Isp,coeff_thrust,is_valid);
    
end

%% functions

% stagnation relations
function p0_p = stagnation_pressure_ratio(M,y)
    p0_p = (1+(y-1)/2 * M.^2).^(y/(y-1));
end

function z0_z = stagnation_density_ratio(M,y)
    z0_z = (1+(y-1)/2 * M.^2).^(1/(y-1));
end

function T0_T = stagnation_temperature_ratio(M,y)
    T0_T = (1+(y-1)/2 * M.^2);
end

function M = find_mach_from_pressure_ratio(p0_p,y)
    tol = 1e-5;
    M_high = 10;
    M_low = 1;
    M_guess = (M_high + M_low)/2;
    p0_p_guess = stagnation_pressure_ratio(M_guess,y);
    while(abs(M_high - M_low) > tol) 
        if p0_p_guess > p0_p % guess is too high
            M_high = M_guess;
        else
            M_low = M_guess; 
        end
        M_guess = (M_low + M_high)/2;
        p0_p_guess = stagnation_pressure_ratio(M_guess,y);
    end
    M = M_guess;
end
% end stagnation relations

% calculating Mach number from area ratio
function M = Mach_from_area(A_At_x,y,isSubsonic)
    tol = 1e-5;
    M = zeros(1,numel(A_At_x));
    if(isSubsonic)
        for i = 1:numel(A_At_x)
            A_At = A_At_x(i);
            M_guess = 0.5;
            M_high = 1;
            M_low = 0;
            A_At_guess = calculate_area_ratio_from_Mach(M_guess,y);
            while(abs(M_high - M_low) > tol)
                if(A_At_guess > A_At) % then M is too small
                    M_low = M_guess;         
                else % M is too big
                    M_high = M_guess;
                end
                M_guess = (M_low + M_high)/2;
                A_At_guess = calculate_area_ratio_from_Mach(M_guess,y);
            end
            M(i) = M_guess;
        end
    else
        for i = 1:numel(A_At_x)
            A_At = A_At_x(i);
            M_guess = 3.5;
            M_high = 10;
            M_low = 1;
            A_At_guess = calculate_area_ratio_from_Mach(M_guess,y);
            while(abs(A_At_guess - A_At) > tol)
                if(A_At_guess > A_At) % then M is too big
                    M_high = M_guess;         
                else % M is too small
                    M_low = M_guess;
                end
                M_guess = (M_low + M_high)/2;
                A_At_guess = calculate_area_ratio_from_Mach(M_guess,y);
            end
            M(i) = M_guess;
        end
    end
end

% calculating area ratio from Mach number - used for iteration in above
% function
function A_At = calculate_area_ratio_from_Mach(M,y)
    A_At = sqrt(   1./M.^2 * ( 2/(y+1) * (1+(y-1)/2.*M^2) ).^( (y+1)/(y-1) )  );
end

% locating a normal shock, assuming it exists
function index = find_normal_shock(A_At_x,y,exit_pressure_ratio,exit_area_ratio)
    q = exit_pressure_ratio^2 * exit_area_ratio^2;
    M_e = sqrt( -1/(y-1) + sqrt( 1/(y-1)^2 + (2/(y-1)) * (2/(y+1))^((y+1)/(y-1)) * 1/q) );

    p0e_pe = stagnation_pressure_ratio(M_e,y);
    

    p02_p01 = p0e_pe * exit_pressure_ratio;
    M1 = ns_jump_M1_from_pressure(p02_p01,y);
    
    A_At = calculate_area_ratio_from_Mach(M1,y);
    index = length(A_At_x);
    for i = 1:length(A_At_x)
        if A_At_x(i) > A_At
            index = i;
            break
        end
    end

end

% calculating M1 based on pressure ratio - used to find ns
function M1 = ns_jump_M1_from_pressure(p02_p01,y)
    tol = 1e-5;
    M1_low = 1;
    M1_high = 10;
    M1_guess = (M1_low + M1_high)/2;
    p02_p01_guess = ns_jump_pressure(M1_guess,y);
    

    while(abs(M1_high - M1_low) > tol)
        if(p02_p01_guess > p02_p01) % then M is too small, because the pressure didn't drop enough
            M1_low = M1_guess;         
        else % M is too big
            M1_high = M1_guess;
        end
        M1_guess = (M1_low + M1_high)/2;
        p02_p01_guess = ns_jump_pressure(M1_guess,y);
    end

    M1 = M1_guess;
end

% ns jump relations
function M2 = ns_jump_M2(M1,y)
    M2 = sqrt( (1 + ( (y-1)/2 * M1^2) ) / ( y*M1^2 - (y-1)/2) );
end

function p02_p01 = ns_jump_pressure(M1,y)
    M2 = ns_jump_M2(M1,y);
    p02_p01 = M1/M2 * ( (2+(y-1)*M2^2) / (2+(y-1)*M1^2) ) ^ ( (y+1) / (2*y-2) );
end
%end ns jump relations

% thrust
function thrust = calc_thrust(p0,A_t,y,p_e,varargin)
    pressure_force = 0;
    if ~isempty(varargin)
        p_a = varargin{1};
        A_e = varargin{2};
        pressure_force = p_a*(p_e/p_a-1)*A_e;
    end
    
    reaction_force = p0*A_t*sqrt( 2*y*y/(y-1) * (2/(y+1))^((y+1)/(y-1)) * (1 - (p_e/p0)^((y-1)/(y)) ) );

    thrust = reaction_force + pressure_force;
end

% isp
function isp = calc_isp(p0,T0,y,R,p_e,varargin)
    g0 = 9.81;
    pressure_contribution = 0;
    if ~isempty(varargin)
        p_a = varargin{1};
        exit_area_ratio = varargin{2};
        pressure_contribution = sqrt(T0)/g0 * p_a/p0 * (p_e/p_a-1) * exit_area_ratio * sqrt(R/y * ((y+1)/2) ^ ((y+1)/(y-1)));
    end
    
    reaction_contribution = 1/g0*sqrt( 2*y*R/(y-1) * T0 * (1 - (p_e/p0)^((y-1)/y) ) );

    isp = reaction_contribution + pressure_contribution;
end

% mdot

% isp
function mdot = calc_mdot(p0,A_t,T0,y,R)
    mdot = p0*A_t/sqrt(T0) * sqrt(y/R * (2/(y+1))^ ((y+1)/(y-1)));
end
