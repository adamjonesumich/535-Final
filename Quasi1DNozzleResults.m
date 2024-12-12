classdef Quasi1DNozzleResults
    %NOZZLERESULTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x_l
        p_ratio_x
        T_ratio_x
        M_x
        u_x
        normal_shock_exists
        shock_location_index
        p0e
        mdot
        thrust
        Isp
        coeff_thrust
        is_valid
    end
    
    methods
        function obj = Quasi1DNozzleResults(x_l,p_ratio_x,T_ratio_x,M_x,u_x,...
                normal_shock_exists,shock_location_index,p0e,mdot,thrust,Isp,coeff_thrust, ...
                is_valid)
            %NOZZLERESULTS Construct an instance of this class
            %   Detailed explanation goes here
            obj.x_l = x_l;
            obj.p_ratio_x = p_ratio_x;
            obj.T_ratio_x = T_ratio_x;
            obj.M_x = M_x;
            obj.u_x = u_x;
            obj.normal_shock_exists = normal_shock_exists;
            obj.shock_location_index = shock_location_index;
            obj.p0e = p0e;
            obj.mdot = mdot;
            obj.thrust = thrust;
            obj.Isp = Isp;
            obj.coeff_thrust = coeff_thrust;
            obj.is_valid = is_valid;
        end
    end
end

