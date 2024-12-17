function rocket = nozzle_adjust(rocket)
% Get needed vals from object 
g0 = 9.81;
R = 8.314;
m_dot = rocket.mass_flow_rate;
p0 = rocket.chamber_pressure;
Me = rocket.exit_mach_number;
gamma = rocket.mixture_gamma;
pa = rocket.average_ambient_pressure;
MW = rocket.mixture_molecular_weight;
T0 = rocket.chamber_temperature;
Ae_At = rocket.area_exit_ratio;
At = rocket.area_throat;

r_t = sqrt(At/pi); % TODO: verify

% full-length

nlines = 100;
thi = 0.1 * pi/180;
[xw,yw,xcl,Mcl] = MinLenNozDes(r_t,Me,gamma,nlines,thi);
L_full = xcl(end);

% L15

L15 = calculate_L15(area_exit_ratio,r_t);

% input correction factor
fprintf("Enter correction factor below." + ...
    " Expansion ratio: %3.1f, L15: %3.1f, Full: %3.1f\n", ...
    rocket.area_exit_ratio, L15, L_full);

rocket.correction_factor = input("");
correction_factor = rocket.correction_factor;

% Calculate Isp 
T = T0*(1+(gamma-1)/2*Me^2)^(-1);
a = sqrt(gamma*R*T/MW);
ue = Me*a;
pe = p0*(1+(gamma-1)/2*Me^2)^(-gamma/(gamma-1));
Ae = Ae_At*At;

FT = (m_dot*ue - (pe-pa)*Ae) * correction_factor;
Isp = FT/(m_dot*g0);

% Set rocket Isp
rocket.specific_impulse = Isp;

end


function L = calculate_L15(area_exit_ratio,r_t)
    
    r_e = 0.4 * r_t; % rule of thumb
    theta = 15 * pi/180;
    L = ( r_t * ( area_exit_ratio ^ 0.5 - 1) + r_e (1 / cos(theta) - 1) ) / tan(theta);

end