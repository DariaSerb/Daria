clear all;
close all;
clc;
  
% Begin Timer
tic

Graphic_display  = 'yes';
% Name_GMSH     = 'plate_hole.msh';
% Name_GMSH     = 'plate_hole_190216.msh';
% Name_GMSH     = 'plate_hole_250416.msh';
% Name_GMSH     = 'plate_hole_270416.msh';
% Name_GMSH     = 'plate_hole_270416_2.msh';
% Name_GMSH     = 'plate_hole_060516.msh';
% Name_GMSH     = 'plate_hole_140616.msh';
% Name_GMSH     = '2Dmesh_hole_0_4_0_2_N_293_r_0_040.msh';
  Name_GMSH     = 'C:\Users\Dasha\Documents\MATLAB\myTest_2D\new\meshes_XFEM\2Dmesh_0_4_0_2_N_258.msh';
% Name_GMSH     = 'C:\Users\Dasha\Documents\MATLAB\myTest_2D\new\meshes_FEM\2Dmesh_hole_0_4_0_2_N_907_r_0_070.msh';
% Name_GMSH     = '2Dmesh_hole_0_4_0_2_N_279_r_0_045.msh';
% Name_GMSH     = '2Dmesh_hole_0_4_0_2_N_261_r_0_050.msh';
Type_LS          = 'Circle';

P                = Initialize_Parameters_2D();
ReferenceElement = SetReferenceElement();

out = Validation(P);
if out == 0
return
end
