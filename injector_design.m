% Injector design script 
% NOTE: This is separate from iterative loop code 
% Assume impinging injector with ring manifold and like-like impingement 

function m_dot_fi, m_dot_oxi, N_f, N_ox = injector_design(pc, m_dot, F_OX_ratio)
    % Geometry 
    D = 1e-3; % m
    Cd = 0.88; % for short tube with rounded entrance injectors 
    delta_p = 0.2*pc;

    % Mass flow rate of fuel and ox
    m_dot_ox = m_dot/4.8;
    m_dot_f = m_dot - m_dot_ox;

    % Mass flow rate through each hole 
    m_dot_oxi = Cd*pi*(D/2)^2*sqrt(2*rho*delta_p);
    m_dot_fi = m_dot_oxi;

    % Number of each injector
    N_f = m_dot_f/m_dot_fi;
    N_ox = m_dpt_ox/m_dot_oxi;

end