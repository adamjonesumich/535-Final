function rocket = nozzle_adjust(rocket)
% Get needed vals from object 
g0 = 9.81;
R = 8.314;
m_dot = rocket.mass_flow_rate;
p0 = rocket.effective_chamber_pressure;
Me = rocket.exit_mach_number;
gamma = rocket.mixture_gamma;
pa_vacuum = rocket.lowest_ambient_pressure;
pa_sea_level = rocket.highest_ambient_pressure;
MW = rocket.mixture_molecular_weight;
T0 = rocket.chamber_temperature;
Ae_At = rocket.area_exit_ratio;
At = rocket.area_throat;

r_t = sqrt(At/pi); % TODO: verify




for i = 1:rocket.number_launch_stages+1
    % full-length
    nlines = 100;
    thi = 0.1 * pi/180;

    % [xw,yw,xcl,Mcl] = MinLenNozDes(r_t(i),Me(i),gamma(i),nlines,thi);
    %L_full(i) = xcl(end);
    L_full(i) = -1;

    % L15

    L15(i) = calculate_L15(Ae_At(i),r_t(i));
    % input correction factor
    fprintf("Enter correction factor below for engine %d." + ...
    " Expansion ratio: %3.1f, L15: %3.1f, Full: %3.1f\n", ...
    i,rocket.area_exit_ratio(i), L15(i), L_full(i));
    
    %rocket.correction_factor(i) = input("");
    rocket.correction_factor(i) = .974;

end







% Calculate Isp 
T = T0.*(1+(gamma-1)/2.*Me.^2).^(-1);
a = sqrt(gamma.*R.*T./MW*1000);
ue = Me.*a;
pe = p0.*(1+(gamma-1)/2.*Me.^2).^(-gamma./(gamma-1));
Ae = Ae_At.*At;
correction_factor = rocket.correction_factor;

FT_vacuum = (m_dot.*ue + (pe-pa_vacuum).*Ae) .* correction_factor;
FT_sea_level = (m_dot.*ue + (pe-pa_sea_level).*Ae) .* correction_factor;
Isp_vacuum = FT_vacuum./(m_dot.*g0);
Isp_sea_level = FT_sea_level./(m_dot.*g0);

% Set rocket Isp
rocket.thrust_vacuum = FT_vacuum;
rocket.thrust_sea_level = FT_sea_level;
rocket.specific_impulse_vacuum = Isp_vacuum;
rocket.specific_impulse_sea_level = Isp_sea_level;

end


function L = calculate_L15(area_exit_ratio,r_t)
    
    r_e = 0.4 * r_t; % rule of thumb
    theta = 15 * pi/180;
    L = ( r_t * ( area_exit_ratio ^ 0.5 - 1) + r_e * (1 / cos(theta) - 1) ) / tan(theta);

end