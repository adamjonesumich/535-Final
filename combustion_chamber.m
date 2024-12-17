function rocket = combustion_chamber(rocket)
    
    for i = 1:rocket.number_launch_stages+1
        % account for Rayleigh losses
        y = rocket.mixture_gamma(i);
        Ac_At = rocket.combustion_area_ratio(i);
        Ae_At = rocket.area_exit_ratio(i);
        A_t = rocket.area_throat(i);
        M_c = Mach_from_area(Ac_At,y,true);
        p0_p0_in = rayleigh_loss(M_c,y);
        p_c = rocket.ideal_chamber_pressure(i);

        rocket.effective_chamber_pressure(i) = p_c * p0_p0_in;
        % calculate engine mass
        rocket.combustion_chamber_mass(i) = calculate_combustion_chamber_mass(A_t,p_c);
        rocket.nozzle_mass(i) = calculate_nozzle_mass(Ae_At,A_t,p_c,Ac_At);
    end

end

function m = calculate_combustion_chamber_mass(A_t,p_c)
    rho_c = 6000; %kg/m3
    o_h = 50e6; %Pa
    L_star = .115; 
    m = A_t * p_c * 2 * rho_c * L_star / o_h;
end

function m = calculate_nozzle_mass(area_exit_ratio,A_t,p_c,Ac_At)
    rho_c = 6000; %kg/m3
    o_h = 50e6; %Pa
    r_t = sqrt(A_t/pi);
    theta = 15 * pi/180;
    A_c = Ac_At * A_t;

    t = p_c * sqrt(A_c/pi) / o_h;
    L15 = calculate_L15(area_exit_ratio,A_t);

    m = rho_c * ((2*pi*t*r_t-pi*t^2)*L15 * pi/2*tan(theta)*L15^2);
    
end

function L = calculate_L15(area_exit_ratio,A_t)
    r_t = sqrt(A_t/pi);
    r_e = 0.4 * r_t; % rule of thumb
    theta = 15 * pi/180;
    L = ( r_t * ( area_exit_ratio ^ 0.5 - 1) + r_e * (1 / cos(theta) - 1) ) / tan(theta);

end

function p0_p0_in = rayleigh_loss(M,y)
    p0_p0_in = (1+(y-1)/2 * M.^2).^(y/(y-1)) / (1 + y * M.^2);
end

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

function A_At = calculate_area_ratio_from_Mach(M,y)
    A_At = sqrt(   1./M.^2 * ( 2/(y+1) * (1+(y-1)/2.*M^2) ).^( (y+1)/(y-1) )  );
end