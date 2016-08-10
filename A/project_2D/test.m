clear all;
close all;
clc;
  
% Begin Timer
tic

Graphic_display  = 'yes';
Name_GMSH        = '2Dmesh_0_6_0_4 - 2.msh';
Type_LS          = 'Circle';

P                = Initialize_Parameters_2D();
ReferenceElement = SetReferenceElement();

out = Validation(P);
