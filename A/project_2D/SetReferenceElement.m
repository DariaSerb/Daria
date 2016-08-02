function Element = SetReferenceElement()
% Element = SetReferenceElement()
% The properties of Reference element:
% the discription of the type of element, degree, number of nodes, nodal coordinates
% integration gauss points and gauss weights, shape functions N, Ns, Nt


    elem  = 1;      % type of element = 1 for triangles
    p     = 1;      % degree 
    ngaus = 7;      % number of Gauss points in each element
    nen   = (p + 1)*(p + 2)/2; 
    
    Element.elem   = elem; 
    Element.nen    = nen; 
    Element.degree = p; 
%   Element.Xe_ref = [1, 0; 0, 1; 0, 0]; % nodal coordinates in (s,t)
    Element.Xe_ref = [0, 0; 1, 0; 0, 1]; % nodal coordinates in (s,t)

    Element.ngaus = ngaus; 
    [zgp,wgp]     = Quadrature(elem, ngaus); 
    [N, Ns, Nt]   = ShapeFunction_T3(zgp); 
    Element.zgp  = zgp; 
    Element.wgp  = wgp; 
    Element.N    = N; 
    Element.Ns   = Ns; 
    Element.Nt   = Nt; 
end
