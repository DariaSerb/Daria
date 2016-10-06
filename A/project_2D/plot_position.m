function plot_position(Nodes,Elements)
% the definition of position

for in=1:size(Elements,1)
 if Elements(in,5) == 1
    Te     = Elements(in,2:4);
    Xe     = Nodes(Te,2:3);
    LSe    = Nodes(Te,4); 
    [p1 p2] = FindPointLS0test(Xe,LSe);
    plot(p1(1),p1(2),'m.',p2(1),p2(2),'m.','LineWidth',8)
    node0 = find(LSe < 0);  
    if length(node0) == 1
            node1 = node0+1; if(node1) == 4, node1 = 1; end
            node2 = node0-1; if(node2) == 0, node2 = 3; end
            P0tr = Te(node0);  
            P1tr = Te(node1);  
            P2tr = Te(node2); 
            
        tab(1,1:2) = [p1(1,:) p1(2,:)];
        tab(2,1:2) = [Nodes(P1tr,2) Nodes(P1tr,3)];
        tab(3,1:2) = [Nodes(P2tr,2) Nodes(P2tr,3)];
        tab(4,1:2) = [p2(1,:) p2(2,:)];
        tab(5,1:2) = [p1(1,:) p1(2,:)];
        plot(tab(:,1),tab(:,2),'c.-') 
    elseif length(node0) == 2 
        node1 = setdiff([1 2 3],node0);
            P0tr = Te(node1); 
            P1tr = Te(node0(1)); 
            P2tr = Te(node0(2));
        
        tab(1,1:2) = [p2(1,:) p2(2,:)];
        tab(2,1:2) = [Nodes(P0tr,2) Nodes(P0tr,3)];
        tab(3,1:2) = [p1(1,:) p1(2,:)];
        tab(4,1:2) = [p2(1,:) p2(2,:)];
        plot(tab(:,1),tab(:,2),'c.-') 
    end
 end    
end
end
