clear all;
close all;
clc;
    
% Begin Timer
tic

Graphic_display  = 'yes';
Name_GMSH        = '2Dmesh_0_6_0_4_N__321.msh';
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
