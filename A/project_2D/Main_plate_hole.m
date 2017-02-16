clear all;
close all;
clc;
    disp('****************************************************************************************************'); 
    disp('Calculation of eigenvalue problem in 2D using XFEM and DD for plate with hole /modelling hole by LS/')
    disp('The transformation between RC and CC is determined by different choice of q(X)/init parameters in m/') 
    disp('the mesh which used is from GMSH')
    disp('****************************************************************************************************'); 
       
% Calculation of eigenvalue problem in 2D using XFEM for plate with hole /modelling hole by LS/
% The transformation between RC and CC is determined by different choice of q(X)
% Integration of [K],[Kdq],[M] and [Mdq] by gauss points using the modification of gauss points for cut elements
% 16/02/2017
% Begin Timer
tic

% to create a file with the mesh
% N = 86;  W = 0.6; H = 0.4;
  N = 321;  W = 0.6; H = 0.4;
% N = 891;  W = 0.6; H = 0.4;
% N = 1211; W = 0.6; H = 0.4;
% N = 1669; W = 0.6; H = 0.4;
% N = 2035; W = 0.6; H = 0.4;
% N = 2098; W = 0.6; H = 0.4;
% N = 3273; W = 0.6; H = 0.4;
% N = 3888; W = 0.6; H = 0.4;
% N = 4237; W = 0.6; H = 0.4;
% N = 4634; W = 0.6; H = 0.4;
% N = 5969; W = 0.6; H = 0.4;
% N = 7562; W = 0.6; H = 0.4; Out of memory. Type HELP MEMORY for your options.

filename = create_name_mesh(N,W,H);

Graphic_display  = 'YES';
  Name_GMSH  = strcat('C:\Users\Dasha\Documents\MATLAB\myTest_2D\new\meshes_XFEM\',filename,'.msh');
% Name_GMSH  = strcat('D:\PhD_thesis\MATLAB\myTest_2D\new\meshes_XFEM\',filename,'.msh');
Type_LS      = 'Circle';

% Init parameters in meters
P = Initialize_Parameters_2D();

% Read the material file
ModeCnt  = P.ModeCnt;
ModeEst  = P.ModeEst;
ModeEst2 = P.ModeEst2;
r0      = P.r0;
dr      = P.dr;
dTau    = P.dTau; 
N_tau   = P.N_tau;
N_dr    = P.N_dr;

% Validation of initial parameters
out = Validation(P);
if out == 0
 return
end

% Read of GMSH information
  [Nodes,Elements] = ReadGMSH(Name_GMSH,Graphic_display);
% [Nodes,Elements] = ReadINP(Name_INP,Graphic_display);

% *************************************************************************
% test experiment for error estimator
% load('NodesCurrent.mat','NodesCurrent') 
% load('ElementsCurrent.mat','ElementsCurrent') 
% plot_mesh(NodesCurrent,ElementsCurrent)
% freqEstNew = EigenValueEst_CurrentConf(NodesCurrent,ElementsCurrent);
  error_indicator(Nodes,Elements,Type_LS) 
% *************************************************************************

% *************************************************************************
% The calculation of estimated frequencies using DD
  [freqEst,LambdaEst] = EigenValueEst_XFEM_DD(Nodes,Elements,Type_LS,Graphic_display);
%   [uEst,NodesCurrent,ElementsCurrent,~,~,~,~] = EigenFunctionEst_XFEM_DD(Nodes,Elements,Type_LS,0.0010);
%   plot_mesh(NodesCurrent,ElementsCurrent)
%   calc_mesh_quality(NodesCurrent,ElementsCurrent)   
% save('NodesCurrent.mat','NodesCurrent') 
% save('ElementsCurrent.mat','ElementsCurrent') 
  
% *************************************************************************

% *************************************************************************
% Calculation of eigenvalue problem numerically (using XFEM) 
% [freqNumCur,uNumCur] = EigenSolutionNum_XFEM(Nodes,Elements,Type_LS,0.0400);
% [freqNumCur,uNumCur] = eig_val(NodesCurrent,ElementsCurrent);
% [freqNum,uNum] = EigenSolutionNum_XFEM(Nodes,Elements,Type_LS,0.0933);
% [freqNum,uNum] = EigenSolutionNum_XFEM(NodesCurrent,Elements,Type_LS,0.0410);
% save('freqNum.mat','freqNum')
% load('NodesCurrent.mat','NodesCurrent') 
% load('ElementsCurrent.mat','Elements') 
% plot_mesh(NodesCurrent,ElementsCurrent)
% [freqNumCur,uNumCur] = eig_val(NodesCurrent,ElementsCurrent);
% *************************************************************************

% [errorL2,errorH1,nelements]=error_convergence(Nodes,Elements,uNumCur,Type_LS);     

flag = 0;
if flag == 1
% The calculation of numerical frequencies 
radius   = zeros(length(dr),1);
freqXFEM = zeros(length(dr),ModeEst);
   for n = 1:N_dr
      radius(n) = r0 + dr(n);   
      [freqNum,uNum] = EigenSolutionNum_XFEM(Nodes,Elements,Type_LS,radius(n));
%     [freqNum,uNum] = EigenSolutionNum_XFEM(NodesCurrent,ElementsCurrent,Type_LS,radius(n));
      freqXFEM(n,:)  = freqNum;
   end
end

plot_freq(radius,freqXFEM,freqEst);
plot_MAC(Nodes,Elements,Type_LS);

% End Timer
computation_time = toc;

disp(['Total computation time: ',num2str(computation_time),' seconds']); 
disp('**************************************************'); 