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
        tab(1,1:2) = [p1(1,:) p1(2,:)];
        tab(2,1:2) = [Nodes(Elements(in,3),2) Nodes(Elements(in,3),3)];
        tab(3,1:2) = [Nodes(Elements(in,4),2) Nodes(Elements(in,4),3)];
        tab(4,1:2) = [p2(1,:) p2(2,:)];
        tab(5,1:2) = [p1(1,:) p1(2,:)];
        plot(tab(:,1),tab(:,2),'c.-') 
    elseif length(node0) == 2 
        tab(1,1:2) = [p1test(1,:) p1test(2,:)];
        tab(2,1:2) = [Nodes(Elements(in,3),2) Nodes(Elements(in,3),3)];
        tab(3,1:2) = [p2(1,:) p2(2,:)];
        tab(4,1:2) = [p1(1,:) p1(2,:)];
        plot(tab(:,1),tab(:,2),'c.-') 
    end
 end    
end
end
