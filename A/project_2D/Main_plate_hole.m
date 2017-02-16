clear all;
close all;
clc;
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
  error_indicator(Nodes,Elements,Type_LS) 
% *************************************************************************
% The calculation of estimated frequencies using DD
  [freqEst,LambdaEst] = EigenValueEst_XFEM_DD(Nodes,Elements,Type_LS,Graphic_display);

flag = 0;
if flag == 1
% The calculation of numerical frequencies 
radius   = zeros(length(dr),1);
freqXFEM = zeros(length(dr),ModeEst);
   for n = 1:N_dr
      radius(n) = r0 + dr(n);   
      [freqNum,uNum] = EigenSolutionNum_XFEM(Nodes,Elements,Type_LS,radius(n));
      freqXFEM(n,:)  = freqNum;
   end
end

plot_freq(radius,freqXFEM,freqEst);

% End Timer
computation_time = toc;

disp(['Total computation time: ',num2str(computation_time),' seconds']); 
disp('**************************************************'); 
