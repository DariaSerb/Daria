function [Nodes,Elements] = ReadGMSH(Name_GMSH,Graphic_display)

fid =fopen(Name_GMSH,'rt');

while ~feof(fid)
    
    line = fgetl(fid);
    if strfind(line,'$Nodes')
        Nb_nodes = str2num(fgetl(fid));
        Nodes    = zeros(Nb_nodes,3);
        for i=1:Nb_nodes
            Tmp          = str2num(fgetl(fid));
            Nodes(i,1:3) = Tmp(1:end-1);
        end
    end

    if strfind(line,'$Elements')
        Nb_elements_2_read = str2num(fgetl(fid));
        Nb_elements        = 1;
        for i=1:Nb_elements_2_read
            Tmp          = str2num(fgetl(fid));
            if Tmp(1,2)==2
                Elements(Nb_elements,1)   = Nb_elements;
                Elements(Nb_elements,2:4) = Tmp(end-2:end);
                Nb_elements               = Nb_elements + 1;
            end
        end
    end
    
end

fclose(fid);

xmin = min(Nodes(:,2));
xmax = max(Nodes(:,2));
ymin = min(Nodes(:,3));
ymax = max(Nodes(:,3));

if strcmp(upper(Graphic_display),'YES')
    
    figure(1);
    hold on;
    for j = 1:size(Elements,1)
        tab(1,1:2) = [Nodes(Elements(j,2),2) Nodes(Elements(j,2),3)];
        tab(2,1:2) = [Nodes(Elements(j,3),2) Nodes(Elements(j,3),3)];
        tab(3,1:2) = [Nodes(Elements(j,4),2) Nodes(Elements(j,4),3)];
        tab(4,1:2) = [Nodes(Elements(j,2),2) Nodes(Elements(j,2),3)];
        plot(tab(:,1),tab(:,2),'g.-') 
        axis([xmin xmax ymin ymax]) 
    end
    hold off;
end