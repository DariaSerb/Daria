function [freqNum,eigvec_XFEM] = EigenSolutionNum(Nodes,Elements,radius)

% Init parameters in meters
P = Initialize_Parameters_2D();

% Read the material file
Data_LS    = P.Data_LS; 
Data_LS(3) = radius; 

Type_LS          = 'Circle';
Graphic_display  = 'NO';

% Initialization and computing Level Set
[Nodes] = ComputeLS(Nodes, Type_LS, Data_LS, Graphic_display);

% Separation the elements in three sets
if strcmp(Type_LS, 'Circle')
    [Elements] = SeparateElements(Elements,Nodes,Graphic_display);

    % Separation the elements inside the void
    ind_In  = 1;
    for ii = 1:size(Elements,1)
        if Elements(ii,5) == -1;
            Elems_In(ind_In) = ii;
            ind_In = ind_In + 1;
        end
    end
else
    Elements(:,5) = zeros(size(Elements,1),1);
end

[freq,eigvec,mui] = eig_val(Nodes,Elements);

freqNum     = freq(1:P.ModeEst);
eigvec_XFEM = eigvec(1:P.ModeEst,:);
