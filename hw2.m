close all
clearvars;
% TODO: combustor static conditions, plot enhancements,
% problem 2: combustor properties
% sim variables
num_points = 10000;
% Problem 1a inputs
T0 = 3200; % K
p0 = 65; % atm
p0 = 101325 * p0; % Pa
MW = 13; % g/mol
MW = 1/1000 * MW; % kg/mol
y = 1.2; % gamma
p_a = .485; % atm
p_a = 101325 * p_a; % Pa
A_t = 1; %TBDTBDTBDTBD
exit_area_ratio = 15;
 % A_e = A_t * area_expansion_ratio;
combustor_area_ratio = 5;
 % A_c = A_t * area_compression_ratio;
R_u = 8.314; % J/molK
R = R_u/MW;

% define area at a given x/l value


params = NozzleParameters(T0,p0,p_a,MW,y,A_t,exit_area_ratio,combustor_area_ratio,num_points);

results = solve_nozzle(params);

plot_distributions(results.x_l,results.p_ratio_x,results.T_ratio_x,results.M_x,results.u_x);

%% Parametric inputs

%% a

exit_area_ratios = logspace(0,2,20);
p0_pa = [10,25,50,100,200,300,500,1000];
thrust_coeff_matrix = zeros(length(exit_area_ratios),length(p0_pa));
for i = 1:length(p0_pa)
    params.p_a = params.p0/p0_pa(i);
    for j = 1:length(exit_area_ratios)
        params.exit_area_ratio = exit_area_ratios(j);
        params = update_area(params);
        results = solve_nozzle(params);
        thrust_coeff_matrix(j,i) = results.coeff_thrust;
    end
end

figure;
hold on;
Legend = cell(1,length(p0_pa));
for i = 1:length(p0_pa)
    semilogx(exit_area_ratios,thrust_coeff_matrix(:,i));
    Legend{i} = string(p0_pa(i));
end
legend(Legend,"Interpreter","latex","FontSize",12);
xlim([1,1000]);
xscale(gca,"log");
ylim([0.6,2]);
xlabel("Exit Area Ratio, [-]","Interpreter","latex","FontSize",12);
ylabel("Thrust Coefficient, [-]","Interpreter","latex","FontSize",12);


%% b

params.exit_area_ratio = 15;
params = update_area(params);
altitudes = linspace(0,50000,101);
thrust_coeff_matrix = zeros(1,length(altitudes));
for i = 1:length(altitudes)
    [T,a,P,rho] = atmosisa(altitudes(i),extended=true);
    params.p_a = P;
    results = solve_nozzle(params);
    thrust_coeff_matrix(i) = results.coeff_thrust;
end

figure;
hold on;
plot(altitudes/1000,thrust_coeff_matrix);

for i = 1:length(altitudes)
    [T,a,P,rho] = atmosisa(altitudes(i),extended=true);
    params.p_a = P;
    
    % find area ratio based on pressure ratio
    M_e = find_mach_from_pressure_ratio(params.p0/params.p_a,y);
    params.exit_area_ratio = calculate_area_ratio_from_Mach(M_e,y);
    params = update_area(params);

    results = solve_nozzle(params);
    thrust_coeff_matrix(i) = results.coeff_thrust;
end

plot(altitudes/1000,thrust_coeff_matrix);
xlabel("Altitude, [km]","Interpreter","latex","FontSize",12);
ylabel("Thrust Coefficient, [-]","Interpreter","latex","FontSize",12);
legend("Area Ratio = 15","Variable Area Ratio","Interpreter","latex","FontSize",12);

%% d

p0_pes = [10,20,30,50,100];
pe_pas = [100,10,1,0.1];

T0_MWs = linspace(1,50000,100);

isp_matrix = zeros(length(p0_pes)*length(pe_pas),length(T0_MWs));

for i = 1:length(p0_pes)
    for j = 1:length(pe_pas)
        for k = 1:length(T0_MWs)
            params.T0 = 3200;
            params.MW = params.T0/T0_MWs(k);
            params.p0 = 1;
            p_e = params.p0/p0_pes(i);
            params.p_a = p_e/pe_pas(j);
    
            % find area ratio based on pressure ratio
            M_e = find_mach_from_pressure_ratio(params.p0/p_e,y);
            params.exit_area_ratio = calculate_area_ratio_from_Mach(M_e,y);
            params = update_area(params);
    
            results = solve_nozzle(params);
            if results.is_valid
                isp_matrix((i-1)*length(pe_pas)+j,k) = results.Isp;
            % else
            %     isp_matrix((i-1)*length(pe_pas)+j,k) = NaN;
            end
        end
    end
end

figure;
hold on;
Legend = cell(1,length(p0_pes)*length(pe_pas));
colors = {'red';[0 0.7 0];'blue';'magenta'};
styles = {'-';'--';':';'-o';'-*'};
for i = 1:length(p0_pes)
    style = styles{i};
    for j = 1:length(pe_pas)
        % if ~isnan(isp_matrix((i-1)*length(pe_pas)+j,1))
            color = colors{j};
            plot(T0_MWs/1000,isp_matrix((i-1)*length(pe_pas)+j,:),style,'Color',color);
            Legend{(i-1)*length(pe_pas)+j} = "p0/pe: " + string(p0_pes(i)) + ", pe/pa: " + string(pe_pas(j));
        % end
    end
end
legend(Legend,"Interpreter","latex","FontSize",12);
xlabel("$\frac{T_0}{MW}$, [K mol / g]","Interpreter","latex","FontSize",12);
ylabel("Isp, [s]","Interpreter","latex","FontSize",12);

%% 2
clearvars;
pressures = {"1 atm", "10 atm", "100 atm"};
figure;

for i = 1:length(pressures)
    j = 1;
    filename = 'cantera1_' + string(i-1) + ".csv";
    T = readtable(filename);
    
    subplot(2,2,j);
    hold on;
    plot(T.r___,T.T_K_);
    j = j+1;
    xlabel("Mixture ratio [-]","Interpreter","latex","FontSize",12);
    ylabel("Adiabatic flame temperature $T_0$ [K]","Interpreter","latex","FontSize",12);
    titlename = "Adiabatic flame temperature as a function of mixture ratio";
    title(titlename, 'Units', 'normalized', 'Position', [0.5 -0.2 0],"Interpreter","latex","FontSize",12);

    subplot(2,2,j);
    hold on;
    plot(T.r___,T.MolarFraction_H2O____);
    j = j+1;
    xlabel("Mixture ratio [-]","Interpreter","latex","FontSize",12);
    ylabel("Molar fraction of water [-]","Interpreter","latex","FontSize",12);
    titlename = "Molar fraction of water as a function of mixture ratio";
    title(titlename, 'Units', 'normalized', 'Position', [0.5 -0.2 0],"Interpreter","latex","FontSize",12);

    subplot(2,2,j);
    hold on;
    plot(T.r___,T.Gamma___);
    j = j+1;
    xlabel("Mixture ratio [-]","Interpreter","latex","FontSize",12);
    ylabel("Product mixture gamma [-]","Interpreter","latex","FontSize",12);
    titlename = "Product mixture gamma as a function of mixture ratio";
    title(titlename, 'Units', 'normalized', 'Position', [0.5 -0.2 0],"Interpreter","latex","FontSize",12);

    subplot(2,2,j);
    hold on;
    plot(T.r___,T.MW_g_mol_);
    xlabel("Mixture ratio [-]","Interpreter","latex","FontSize",12);
    ylabel("Product mixture molecular weight $MW$ [g/mol]","Interpreter","latex","FontSize",12);
    titlename = "Product mixture molecular weight as a function of mixture ratio";
    title(titlename, 'Units', 'normalized', 'Position', [0.5 -0.2 0],"Interpreter","latex","FontSize",12);
end

for j = 1:4
    subplot(2,2,j);
    hold on;
    legend(pressures,"Interpreter","latex","FontSize",12);
end

%% 2g-o

clearvars;
rs = {"2","3","4","5","6","7"};
A_t = 0.01;
p_a = 0.26; % atm
p_e = p_a * 101325;
R_u = 8.134;
pressure_ratios = [0.1,1,10];

f_thrust = figure;
f_coeff = figure;
f_isp = figure;
f_area = figure;
f_mach = figure;
f_mdot = figure;

fs = [f_thrust,f_coeff,f_isp,f_area,f_mach,f_mdot];

colors = [1 0 0;
    0.7 0.4 0;
    0.1 0.8 0;
    0 .5 0.5;
    0.1 0.1 0.8;
    0.5 0 0.5];

colors = [0.8 0 0.1;
    0.5 0 0.3;
    0.3 0 0.7;
    0 0 0.8;
    0 0.3 0.8;
    0 0.4 0.8];

for i = 1:length(rs)
    color = colors(i,:);
    filename = 'cantera2_' + string(i-1) + ".csv";
    T = readtable(filename);
    T = table2array(T);
    
    pressures = T(:,1)*101325;
    temps = T(:,2);
    MWs = T(:,3)/1000;
    gammas = T(:,4);
    for k = 1:3
        
        p_a = p_e / pressure_ratios(k);

        thrusts = zeros(1,length(pressures));
        thrust_coeffs = zeros(1,length(pressures));
        isps = zeros(1,length(pressures));
        exit_area_ratios = zeros(1,length(pressures));
        machs = zeros(1,length(pressures));
        mdots = zeros(1,length(pressures));
        
        % thrust coefficient
        for j = 1:length(pressures)
            
            p0 = pressures(j);
            T0 = temps(j);
            y = gammas(j);
            MW = MWs(j);


            R = R_u / MW;
            
            mach_equal = find_mach_from_pressure_ratio(p0/p_e,y); % for exit area ratio
            exit_area_ratios(j) = calculate_area_ratio_from_Mach(mach_equal,y);
            A_e = exit_area_ratios(j) * A_t;
            
            

            thrusts(j) = calc_thrust(p0,A_t,y,p_e,p_a,A_e);
            thrust_coeffs(j) = thrusts(j)/(p0*A_t);
            isps(j) = calc_isp(p0,T0,y,R,p_e,p_a,exit_area_ratios(j));
            machs(j) = find_mach_from_pressure_ratio(p0/p_e,y);
            mdots(j) = calc_mdot(p0,A_t,T0,y,R);

            % combustor_area_ratio = 5;
            % exit_area_ratio = exit_area_ratios(j);
            % num_points = 1000;
            % 
            % params = NozzleParameters(T0,p0,p_a,MW,y,A_t,exit_area_ratio,combustor_area_ratio,num_points);
            % results = solve_nozzle(params);
            % thrusts(j) = results.thrust;
            % thrust_coeffs(j) = results.coeff_thrust;
            % isps(j) = results.Isp;
            % machs(j) = results.M_x(end);
            % mdots(j) = results.mdot;

            

        end
    
        set(0,'CurrentFigure',f_thrust);
        subplot(1,3,k)
        hold on;

        plot(pressures/101325,thrusts/1e6,'Color',color);
        xlabel("Combustion pressure [atm]","Interpreter","latex","FontSize",12);
        ylabel("Thrust [MN]","Interpreter","latex","FontSize",12);
        ylim([0,0.5])
        
        
        set(0,'CurrentFigure',f_coeff);
        subplot(1,3,k)
        hold on;

        plot(pressures/101325,thrust_coeffs,'Color',color);
        xlabel("Combustion pressure [atm]","Interpreter","latex","FontSize",12);
        ylabel("Thrust coefficient [-]","Interpreter","latex","FontSize",12);
        ylim([0,2])
    
        set(0,'CurrentFigure',f_isp);
        subplot(1,3,k)
        hold on;

        plot(pressures/101325,isps,'Color',color);
        xlabel("Combustion pressure [atm]","Interpreter","latex","FontSize",12);
        ylabel("Isp [s]","Interpreter","latex","FontSize",12);
        ylim([0,380])
    
        set(0,'CurrentFigure',f_area);
        subplot(1,3,k)
        hold on;

        plot(pressures/101325,exit_area_ratios,'Color',color);
        xlabel("Combustion pressure [atm]","Interpreter","latex","FontSize",12);
        ylabel("Exit area ratio [-]","Interpreter","latex","FontSize",12);
        ylim([0,60])
    
        set(0,'CurrentFigure',f_mach);
        subplot(1,3,k)
        hold on;
        

        plot(pressures/101325,machs,'Color',color);
        xlabel("Combustion pressure [atm]","Interpreter","latex","FontSize",12);
        ylabel("Exit mach number [-]","Interpreter","latex","FontSize",12);
        ylim([0,4.8])

        set(0,'CurrentFigure',f_mdot);
        subplot(1,3,k)
        hold on;
        

        plot(pressures/101325,mdots,'Color',color);
        xlabel("Combustion pressure [atm]","Interpreter","latex","FontSize",12);
        ylabel("Mass flow rate [kg/s]","Interpreter","latex","FontSize",12);
        ylim([0,140])

    end 
    
end

for j = 1:6
    set(0,'CurrentFigure',fs(j));
    for k = 1:3
        subplot(1,3,k);
        hold on;
        title("Exit-Ambient Pressure Ratio: " + string(pressure_ratios(k)), 'Units', 'normalized', 'Position', [0.5 -0.1 0],"Interpreter","latex","FontSize",12);
        legend(rs,'location','southeast');
    end
end
%% functions
function results = solve_nozzle(params)
    
    % unpack
    T0 = params.T0;
    p0 = params.p0;
    p_a = params.p_a;
    MW = params.MW;
    y = params.y;
    A_t = params.A_t;
    exit_area_ratio = params.exit_area_ratio;
    combustor_area_ratio = params.combustor_area_ratio;
    num_points = params.num_points;
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

    
    results = NozzleResults(x_l,p_ratio_x,T_ratio_x,M_x,u_x,normal_shock_exists,...
        shock_location_index,p0e,mdot,thrust,Isp,coeff_thrust,is_valid);
    
end

function plot_distributions(x_l,p_ratio_x,T_ratio_x,M_x,u_x)

    figure(1);
    hold on;

    plot(x_l,p_ratio_x);
    % title("Pressure Ratio Distribution");
    xlabel("Normalized Axial Location, [-]","Interpreter","latex","FontSize",12);
    ylabel("Pressure ratio (\textit{p/p0}), [-]","Interpreter","latex","FontSize",12);

    figure(2);
    hold on;
    
    plot(x_l,T_ratio_x);
    % title("Temperature Ratio Distribution");
    xlabel("Normalized Axial Location, [-]","Interpreter","latex","FontSize",12);
    ylabel("Temperature Ratio (\textit{T/T0}), [-]","Interpreter","latex","FontSize",12);

    figure(3);
    hold on;
    
    plot(x_l,M_x);
    % title("Mach Distribution");
    xlabel("Normalized Axial Location, [-]","Interpreter","latex","FontSize",12);
    ylabel("Mach Number (M), [-]","Interpreter","latex","FontSize",12);
    
    figure(4);
    hold on;
    
    plot(x_l,u_x);
    % title("Velocity Distribution");
    xlabel("Normalized Axial Location, [-]","Interpreter","latex","FontSize",12);
    ylabel("Axial Velocity (u), [m/s]","Interpreter","latex","FontSize",12);

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


