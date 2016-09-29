clear all;
close all;
clc;
    
% Begin Timer
tic

% to create a file with the mesh
N = 321; W = 0.6; H = 0.4;
filename = create_name_mesh(N,W,H);

Graphic_display  = 'YES';
Name_GMSH        = strcat('C:\Users\Dasha\Documents\MATLAB\myTest_2D\new\meshes_XFEM\',filename,'.msh');
Type_LS          = 'Circle';

% Init parameters in meters
P = Initialize_Parameters_2D();

% Read the material file
ModeCnt = P.ModeCnt;
ModeEst = P.ModeEst;
r0      = P.r0;
dr      = P.dr;
dTau    = P.dTau; 
N_tau   = P.N_tau;

% Validation of initial parameters
out = Validation(P);
if out == 0
 return
end

% Read of GMSH information
[Nodes,Elements]  = ReadGMSH(Name_GMSH,Graphic_display);

[freqNum, uNum] = EigenSolutionNum_XFEM(Nodes,Elements,0.07);

for n = 1:N_dr
   radius(n) = r0 + dr(n);   
   [freqNum, uNum] = EigenSolutionNum_XFEM(Nodes,Elements,Type_LS,radius(n));
   freqNum(n,:) = freqNum;
end

for n = 1:ModeEst
   figure(n + 1);
   plot(radius, freqNum(:,n), radius, freqEst(:,n), 'LineWidth', 3);
   legend('freqNum','freqEst',3)
   grid on;
   title(strcat('Natural frequency # ', num2str(n), ' for the position of the initial radius', ' r_0 = ',  num2str(r0), ' [m]'),'FontSize', 12);
   xlabel('r_0 [m]','FontSize', 12);
   ylabel('Natural frequency, f [Hz]','FontSize', 12);
end
