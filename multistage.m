function rocket = multistage(rocket)
    
    % sim values
    dv_tot = rocket.total_delta_v(1);
    g0 = 9.807;

    
    Isp = rocket.specific_impulse_sea_level(1:end);
    ISP = Isp(1:rocket.number_launch_stages); % lower stages Isp

    % here, I just manually inputs what the structural coefficients are;
    % this needs to be replaced with a loop that looks up e_i based on the
    % value of R_i

    % We do need an initial guess based on how this is set up, 
    % so should I leave the preset epsilons here?
    
    e_1 = 0.032;
    e_2 = 0.035;
    e_3 = 0.055;
    e = [e_1;e_2;e_3];

    e_residual = e;
    tol = 1e-4;

    %do_loop = true;
    while(norm(e_residual,Inf) > tol) % TODO: remove "do_loop"
        %do_loop = false;
        
        % -----------lagrange multiplier math--------------
        % the Lagrange system forms a set of 4 equations in 4 unknowns; I
        % solved the system for R1, thus the loop is able to calculate the
        % constraint in terms of only R1. The residual should be 0 (based
        % on the constraint equality) so the loop searches until the R1
        % value basically returns 0. It then calculates the other mass
        % ratios based on R1, and returns them. This is only for ONE set of
        % structural coefficients; the outer loop updates the structural
        % coefficients and runs again until they converge.

        R_residual = 1;
        R_max = 5;
        R_min = 1;
        R_guess = R_min + (R_max-R_min)/2;
        
        
        while abs(R_residual) > tol
            R_residual = multistage_residual(Isp,e,1,dv_tot,g0,R_guess);
            if(R_residual > 0)
                R_min = R_guess;
            else
                R_max = R_guess;
            end
            R_guess = R_min + (R_max-R_min)/2;
        end
        
        Rs = [0;0;0];
        Rs(1) = R_guess;
        
        lambda = @(n) Rs(n)*(1 + e(n)) / (Isp(n) * g0 * ( 2 + Rs(n)*(1-3*e(n)) - Rs(n)^2*(e(n)^2-e(n)) ) );
        lambda_to_use = lambda(1);
        
        A = @(n) e(n)^2-e(n);
        B = @(n) 1-3*e(n)-(1+e(n))/(lambda_to_use*Isp(n)*g0);
        
        Rlow = @(n) (-B(n) + sqrt(B(n)^2-8*A(n)))/(2*A(n));
        Rhigh = @(n) (-B(n) - sqrt(B(n)^2-8*A(n)))/(2*A(n));
        R = @(n) max(Rlow(n),Rhigh(n));
        
        Rs(2) = R(2);
        Rs(3) = R(3);


        e_prev = e;
        e = fetch_e(Rs,ISP);

        e_residual = e - e_prev;

    end
    
    % calculate parameters based on finalized mass ratios and structural
    % coefficients

    dvs = [0;0;0];
    for i = 1:3
        dvs(i) = Isp(i)*g0*log(Rs(i));
    end
    
    prs = [0;0;0];
    for i = 1:3
        prs(i) = (1-Rs(i)*e(i))/(Rs(i)-1);
    end
   




    % upper stage math
    
    dv_upper = rocket.total_delta_v(end);
    Isp_upper = rocket.specific_impulse_sea_level(end);
    m_pl_upper = rocket.upper_stage_mass + rocket.combustion_chamber_mass(end) + rocket.nozzle_mass(end);
    e_upper = fetch_e_from_dv(dv_upper); % TODO: fix with real value
    
    C = exp(-dv_upper/(Isp_upper*g0));
    
    m_p_upper = m_pl_upper*(e_upper/(1-e_upper)*(C-1)+C)^(-1)*(1-C);
    m_s_upper = m_p_upper*(e_upper/(1-e_upper));
    
    m_upper = m_pl_upper+m_p_upper+m_s_upper;
    ms = [0;0;0;m_upper];
    for i = 3:-1:1
        ms(i) = ms(i+1)*(prs(i)+1)/prs(i) + rocket.combustion_chamber_mass(i) + rocket.nozzle_mass(i);
    end

    rocket.stage_payload_ratios = prs;
    rocket.stage_mass_ratios = Rs;
    rocket.stage_masses = ms;
    rocket.stage_delta_vs = dvs;

    function e = fetch_e(R,ISP)
        % Inputs:
        %    R = Initial to burn-out mass ratio of each stage, vector
        %    ISP = Specific impulse of each stage, vector
    
        % Outputs:
        %    e = structural coefficients for each stage, vector
    
        % Need to convert R1,R2,R3 into delta v
        dv1 = 9.807*ISP(1)*log(R(1));
        dv2 = 9.807*ISP(2)*log(R(2));
        dv3 = 9.807*ISP(3)*log(R(3));
        
        % Curve fits for epsilon as a function of x, where x = delta v
        %H2fit = [91.3825   -1.0164    0.0586];
        HCfit = [76.7767   -1.1191    0.0266];
        %epsfunc_H2 = @(x)H2fit(1)*x.^H2fit(2)+H2fit(3); % H2 fuel
        epsfunc_HC = @(x)HCfit(1)*x.^HCfit(2)+HCfit(3); % Hydrocarbon fuel
    
        e1 = epsfunc_HC(dv1);
        e2 = epsfunc_HC(dv2);
        e3 = epsfunc_HC(dv3);
        e = [e1;e2;e3];
        
        %e = zeros(size(R)) + 0.04; % ?
    end

    function e = fetch_e_from_dv(x)
        HCfit = [76.7767   -1.1191    0.0266];
        %epsfunc_H2 = @(x)H2fit(1)*x.^H2fit(2)+H2fit(3); % H2 fuel
        e = HCfit(1)*x.^HCfit(2)+HCfit(3);
    end
end