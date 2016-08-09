classdef Validation
    % validation of the initial parameters
      
    properties
      P = Initialize_Parameters_2D();
    end
    
    methods (Static)
      function out = ValidMethod()
        if (P.radiusout < 0.5*P.lx)&&(P.radiusout < 0.5*P.ly)   
        out = true;
        else
        out = false;    
        end
      end
    end
end
