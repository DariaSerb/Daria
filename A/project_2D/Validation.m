classdef Validation
    % validation of the initial parameters
      
    properties
      P = Initialize_Parameters_2D();
    end
    
    methods 
      function ValidMethod(obj,arg1)
       obj = pi*obj.R1^2;    
    
      end
    end
end

