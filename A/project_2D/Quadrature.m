function [zgp,wgp] = Quadrature(elem, ngaus) 
% [zgp,wgp] = Quadrature(elem, ngaus)
% Input:    
% elem:   type of element = 1 for triangles
% ngaus:  number of Gauss points in each element
% Output:   
% zgp, wgp: Gauss points and weights on the reference element

% Symmetric quadrature for the unit triangle
if elem == 1
    if ngaus == 1          % degree 1
        pos1  = 1/3;
        zgp   = [pos1   pos1]; 
        wgp   = 1/2; 
    elseif ngaus == 3      % degree 2
        pos1  = 1/2;
        zgp   = [pos1   pos1
                 0      pos1
                 pos1   0];
        pes1  = 1/6;
        wgp   = [pes1   pes1   pes1];
    elseif ngaus == 4      % degree 3
        zgp   = [1/3    1/3
                 0.6    0.2
                 0.2    0.6
                 0.2    0.2];
       wgp    = [-27/96   25/96   25/96   25/96]; 
    elseif ngaus == 7      % degree 5
        a     = 0.101286507323456338800987361915123;
        b     = 0.470142064105115089770441209513447;
        P1    = 0.0629695902724135762978419727500906;
        P2    = 0.0661970763942530903688246939165759;
        zgp   = [a        a
                 a        1 - 2*a
                 1 - 2*a  a
                 b        b
                 b        1 - 2*b
                 1 - 2*b  b
                 1/3      1/3];
        wgp   = [P1, P1, P1, P2, P2, P2, 0.1125]; 
    elseif ngaus == 12   % degree 6
        a     = 0.0630890144915022283403316028708191;
        b     = 0.249286745170910421291638553107019;
        c     = 0.0531450498448169473532496716313981; 
        d     = 0.310352451033784405416607733956552;
        P1    = 0.0254224531851034084604684045534344;
        P2    = 0.0583931378631896830126448056927897;
        P3    = 0.0414255378091867875967767282102212;
        zgp   = [a          a  
                 a          1 - 2*a
                 1 - 2*a    a
                 b          b  
                 b          1 - 2*b
                 1 - 2*b    b
                 c          d
                 c          1 - c - d
                 1 - c - d  c
                 1 - c - d  d
                 d          1 - c - d
                 d          c]; 
        wgp   = [P1, P1, P1, P2, P2, P2, P3, P3, P3, P3, P3, P3];
    else
    error('Unavailable quadrature') 
   end
end
