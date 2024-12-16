function rocket = post_nozzle_analysis(rocket, correction_factor)
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