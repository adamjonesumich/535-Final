classdef Quasi1DNozzleParameters
    %SOLVERPARAMETERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        T0
        p0
        p_a
        MW
        y
        A_t
        exit_area_ratio
        combustor_area_ratio
        num_points;
        x_l_c
        x_l_d
        area_ratio_c
        area_ratio_d
    end
    
    methods
        function obj = Quasi1DNozzleParameters(T0,p0,p_a,MW,y,A_t,exit_area_ratio,...
                combustor_area_ratio,num_points)
            %SOLVERPARAMETERS Construct an instance of this class
            %   Detailed explanation goes here
            obj.T0 = T0;
            obj.p0 = p0;
            obj.p_a = p_a;
            obj.MW = MW;
            obj.y = y;
            obj.A_t = A_t;
            obj.exit_area_ratio = exit_area_ratio;
            obj.combustor_area_ratio = combustor_area_ratio;
            obj.num_points = num_points;

            x_l_c = linspace(-1,0,num_points);
            x_l_d = linspace(0,1,num_points);
           
            obj.x_l_c = x_l_c;
            obj.x_l_d = x_l_d;
            obj.area_ratio_c = ( (1-(obj.combustor_area_ratio-1)*obj.x_l_c) );
            obj.area_ratio_d = ( (1+(obj.exit_area_ratio-1)*obj.x_l_d) );
        end
        
        function obj = update_area(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.area_ratio_c = ( (1-(obj.combustor_area_ratio-1)*obj.x_l_c) );
            obj.area_ratio_d = ( (1+(obj.exit_area_ratio-1)*obj.x_l_d) );
        end
    end
end

