function [zgp,wgp] = ModifyQuadrature(LSe,ReferenceElement,zgp_tri,wgp_tri)

% [zgp,wgp] = ModifyQuadrature(LSe,RefElement,zgp_tri,wgp_tri)
% Integration points and weights for an element cut by the interface
% Input: 
% LSe:              leves set values at the element's nodes
% ReferenceElement: reference element properties
% zgp_tri,wgp_tri:  quadrature in the reference element

elem   = ReferenceElement.elem;
degree = ReferenceElement.degree;
Xe_ref = ReferenceElement.Xe_ref; 


tol = 1e-15;

if elem == 1 && degree == 1
    node0 = find(abs(LSe) < tol); LSe(node0) = 0;
    if LSe >= 0
        zgp = zgp_tri; wgp = wgp_tri;
    elseif LSe <= 0
        zgp = zgp_tri; wgp = 0 * wgp_tri;
    else
     node0 = find(LSe <= 0);   
        if length(node0) == 1
            node1 = node0+1; if(node1) == 4, node1 = 1; end
            node2 = node0-1; if(node2) == 0, node2 = 3; end
            P0 = Xe_ref(node0,:);  d0 = LSe(node0); 
            P1 = Xe_ref(node1,:);  d1 = LSe(node1); 
            P2 = Xe_ref(node2,:);  d2 = LSe(node2); 
            
            PInt1 = FindPointLS0(P0,d0,P1,d1); 
            PInt2 = FindPointLS0(P0,d0,P2,d2); 
            
%           [zgp1,wgp1] = ModQuadTriangle([P1;P2;PInt1],zgp_tri,wgp_tri);
%           [zgp2,wgp2] = ModQuadTriangle([P2;PInt1;PInt2],zgp_tri,wgp_tri);
%           [zgp3,wgp3] = ModQuadTriangle([PInt2;PInt1;P0],zgp_tri,wgp_tri);

            [zgp1,wgp1] = ModQuadTriangle([PInt1;P1;      P2],zgp_tri,wgp_tri);
            [zgp2,wgp2] = ModQuadTriangle([PInt2;PInt1;   P2],zgp_tri,wgp_tri);
            [zgp3,wgp3] = ModQuadTriangle([P0;   PInt1;PInt2],zgp_tri,wgp_tri);
            
            zgp = [zgp1; zgp2; zgp3];
            wgp = [wgp1; wgp2; wgp3];
         
        elseif length(node0) == 2
            node1 = setdiff([1 2 3],node0);
            P0 = Xe_ref(node1,:);  d0 = LSe(node1); 
            P1 = Xe_ref(node0(1),:);  d1 = LSe(node0(1)); 
            P2 = Xe_ref(node0(2),:);  d2 = LSe(node0(2)); 
            
            PInt1 = FindPointLS0(P0,d0,P1,d1); 
            PInt2 = FindPointLS0(P0,d0,P2,d2); 
            
            [zgp,wgp] = ModQuadTriangle([P0;PInt1;PInt2],zgp_tri,wgp_tri);
        end
     end    
end

    figure(3);
    hold on;
    plot(Xe_ref(:,1), Xe_ref(:,2), 'm.-','LineWidth', 3) 
    plot(PInt1(1), PInt1(2), 'r.', PInt2(1), PInt2(2),'b.','LineWidth', 3)
    plot(zgp1(:,1), zgp1(:,2), 'k.', 'LineWidth', 4)
    plot(zgp2(:,1), zgp2(:,2), 'r.', 'LineWidth', 4)
    plot(zgp3(:,1), zgp3(:,2), 'b.', 'LineWidth', 4)
    hold off;

end