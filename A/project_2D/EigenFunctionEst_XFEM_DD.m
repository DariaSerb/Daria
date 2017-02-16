function [uest,NodesCurrent,Elements,eigvecnorm,G,Knorm,Mnorm,ratio_norm] = EigenFunctionEst_XFEM_DD(Nodes,Elements,Type_LS,dTau)
global P;
P = Initialize_Parameters_2D();

% Read the material file
ModeEst = P.ModeEst;
Data_LS = P.Data_LS; 
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

% u = InitialApproximation;
[lambda,u,DerivLambda,DerivShapes,dof_out_BC_hole,C] = DirectDeriv(Nodes,Elements);

lambda = lambda';
LambdaEst = zeros(1,ModeEst);
  for n = 1:ModeEst
   LambdaEst(1,n) = dTau*DerivLambda(1,n) + lambda(1,n);    
  end

% u = u + alp * u * dTau; the estimated eigenshapes 
uesttemp = u(:,1:ModeEst) + dTau*DerivShapes(:,1:ModeEst);

[eigvecnorm,G,Knorm,Mnorm,ratio_norm,valid_F] = calc_norm(Nodes,Elements,LambdaEst,uesttemp,dof_out_BC_hole,dTau);
maxtmp = cellfun(@(x) max(x),valid_F,'Uniform',1)';

% to reestablish the ovelall eihenvector
uesttemp(dof_out_BC_hole,1:ModeEst) = uesttemp;  
uesttemp(C,:) = uesttemp(C,:) * 0;

% normalization of eigenshapes 
for in = 1:ModeEst
uest(:,in) = uesttemp(:,in)/max(uesttemp(:,in));
end

NodesCurrent = change_configuration(Nodes,dTau);

flag = 0;
if flag == 1 
 % plot of the initial configuration and current one
    figure(2);
    hold on;
    for j = 1:size(Elements,1)
        % Plot the FEM elements
        if Elements(j,5) == 0
        % plot of the current configuration    
        tab(1,1:2) = [NodesCurrent(Elements(j,2),2) NodesCurrent(Elements(j,2),3)];
        tab(2,1:2) = [NodesCurrent(Elements(j,3),2) NodesCurrent(Elements(j,3),3)];
        tab(3,1:2) = [NodesCurrent(Elements(j,4),2) NodesCurrent(Elements(j,4),3)];
        tab(4,1:2) = [NodesCurrent(Elements(j,2),2) NodesCurrent(Elements(j,2),3)];
        plot(tab(:,1),tab(:,2),'b.-') 
        end
        % Plot the XFEM elements
        if Elements(j,5) == 1
       % plot of the current configuration
        tab(1,1:2) = [NodesCurrent(Elements(j,2),2) NodesCurrent(Elements(j,2),3)];
        tab(2,1:2) = [NodesCurrent(Elements(j,3),2) NodesCurrent(Elements(j,3),3)];
        tab(3,1:2) = [NodesCurrent(Elements(j,4),2) NodesCurrent(Elements(j,4),3)];
        tab(4,1:2) = [NodesCurrent(Elements(j,2),2) NodesCurrent(Elements(j,2),3)];
        plot(tab(:,1),tab(:,2),'m.-') 
        end
        if Elements(j,5) == -1
       % plot of the current configuration
        tab(1,1:2) = [NodesCurrent(Elements(j,2),2) NodesCurrent(Elements(j,2),3)];
        tab(2,1:2) = [NodesCurrent(Elements(j,3),2) NodesCurrent(Elements(j,3),3)];
        tab(3,1:2) = [NodesCurrent(Elements(j,4),2) NodesCurrent(Elements(j,4),3)];
        tab(4,1:2) = [NodesCurrent(Elements(j,2),2) NodesCurrent(Elements(j,2),3)];
        plot(tab(:,1),tab(:,2),'c.-') 
        end
%       hold off;      
        axis([min(NodesCurrent(:,2)) max(NodesCurrent(:,2)) min(NodesCurrent(:,3)) max(NodesCurrent(:,3))]) 
        title(strcat('The current configuration for tau = ',  num2str(dTau)),'FontSize', 14);

   end
end

Data_LS(3) = P.r0 + dTau; 
% [Nodes] = ComputeLS(Nodes,Type_LS,Data_LS,'NO');
[NodesTemp] = ComputeLS(NodesCurrent,Type_LS,Data_LS,'YES');
if flag == 1
    figure(2);
    % plot real LS for total radius    
    hold on;
    if strcmp(Type_LS, 'Circle')
        Angle = [0:2*pi/100:2*pi];
%         xold  = P.r0*cos(Angle) + P.pos_x_center;
%         yold  = P.r0*sin(Angle) + P.pos_y_center;
        x     = Data_LS(3)*cos(Angle) + P.pos_x_center;
        y     = Data_LS(3)*sin(Angle) + P.pos_y_center;
%         plot(xold,yold,'r.-');
        plot(x,y,'k.-');
        axis([min(Nodes(:,2)) max(Nodes(:,2)) min(Nodes(:,3)) max(Nodes(:,3))]) 
    end
end 
% Separation the elements in three sets
NodesCurrent(:,4) = NodesTemp(:,4);
   
% the definition of the first five estimated modes
if flag == 1 
  for in = 1:ModeEst;
    plot_displacementLS(in,uesttemp(:,in),dof_out_BC_hole,C,NodesCurrent,Elements);
  end
end

end
