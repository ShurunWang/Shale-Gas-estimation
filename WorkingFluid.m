classdef WorkingFluid < handle
    % WorkingFluid 
    
    properties
        name char
        % oC
        temperature double
        % kPa
        pressure double
        massfraction double
    end
    
    properties (Constant)
      T0 = 273.15;
   end
    
    methods
        function obj = WorkingFluid(fluid_name,t,p,x)
            % 构造此类的实例
            obj.name = fluid_name;
            % oC
            obj.temperature = t;
            % kPa
            obj.pressure = p;
            if nargin == 4
                obj.massfraction = x;
            else
                obj.massfraction = 1;
            end
        end
        
        function nu = kineticViscosity(obj)
            % 计算动力粘度，单位[m^2/s]
            T = obj.temperature + obj.T0;
            nu = refproparray('$','T',T,'P',obj.pressure,obj.name);
            nu = nu/10000;
        end
        
        function miu = dynamicViscosity(obj)
            % [Pa*s]
            T = obj.temperature + obj.T0;
            miu = refproparray('V','T',T,'P',obj.pressure,obj.name);
        end
        
        function cg = isothrmCompress(obj)
            % [1/kPa]
            T = obj.temperature + obj.T0;
            cg = refproparray('=','T',T,'P',obj.pressure,obj.name);
        end
        
        function Z = compFactor(obj)
            % [1]
            T = obj.temperature + obj.T0;
            Z = refproparray('Z','T',T,'P',obj.pressure,obj.name);
        end
        
        function rho = density(obj)
            % [kg/m^3]
            T = obj.temperature + obj.T0;
            rho = refproparray('D','T',T,'P',obj.pressure,obj.name);
        end
            
    end
end

