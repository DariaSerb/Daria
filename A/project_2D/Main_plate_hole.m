clear all;
close all;
clc;
    disp('****************************************************************************************************'); 
    disp('Calculation of eigenvalue problem in 2D using XFEM and DD for plate with hole /modelling hole by LS/')
    disp('The transformation between RC and CC is determined by different choice of q(X)/init parameters in m/') 
    disp('****************************************************************************************************'); 
       
% Calculation of eigenvalue problem in 2D using XFEM for plate with hole /modelling hole by LS/
% The transformation between RC and CC is determined by different choice of q(X)
% Integration of [K],[Kdq],[M] and [Mdq] by gauss points using the modification of gauss points for cut elements
% Begin Timer
tic

Graphic_display  = 'yes';
Name_GMSH        = 'C:\Users\Dasha\Documents\MATLAB\myTest_2D\new\meshes_XFEM\2Dmesh_0_6_0_4_N__321.msh';
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

freqEst = EigenValueEst_XFEM_DD(Nodes,Elements,Graphic_display);
% [freqNum, uNum] = EigenSolutionNum_XFEM(Nodes,Elements,0.07);

% uest  = EigenFunctionEst_XFEM_DD(Nodes,Elements,dTau);

% mac = zeros(ModeEst, N_tau);
ModeEst2 = ModeEst*ModeEst;
macfin   = zeros(ModeEst2,N_tau-3);
radius   = zeros(length(dr),1);

for n = 1:N_tau - 3
  radius(n) = r0 + dr(n);  
  uEst = EigenFunctionEst_XFEM_DD(Nodes,Elements,dTau(n));
  [freqNum, uNum] = EigenSolutionNum_XFEM(Nodes,Elements,radius(n));
  mac  = MAC(uNum, uEst);   
  macfin(:,n) = reshape(mac',1,ModeEst2);
end

figure(3);
plot(radius(1:7), macfin(1:5,:), 'LineWidth', 2);
legend('Uappr1, Unum1','Uappr1, Unum2','Uappr1, Unum3','Uappr1, Unum4','Uappr1, Unum5',3)
grid on;
axis([-inf inf -inf 1]);
title(['Modal Assurance Criterium (Uappr1, Unum)']);
xlabel('radius');
ylabel('log(MAC)');

figure(4);
plot(radius(1:7), macfin(6:10,:), 'LineWidth', 2);
legend('Uappr2, Unum1','Uappr2, Unum2','Uappr2, Unum3','Uappr2, Unum4','Uappr2, Unum5',3)
grid on;
axis([-inf inf -inf 1]);
title(['Modal Assurance Criterium (Uappr2, Unum)']);
xlabel('radius');
ylabel('log(MAC)');

figure(5);
plot(radius(1:7), macfin(11:15,:), 'LineWidth', 2);
legend('Uappr3, Unum1','Uappr3, Unum2','Uappr3, Unum3','Uappr3, Unum4','Uappr3, Unum5',3)
grid on;
axis([-inf inf -inf 1]);
title(['Modal Assurance Criterium (Uappr3, Unum)']);
xlabel('radius');
ylabel('log(MAC)');


% End Timer
computation_time = toc;

disp(['Total computation time: ',num2str(computation_time),' seconds']); 
disp('**************************************************'); 

















% % mac = zeros(ModeEst, N_tau);
% ModeEst2 = ModeEst*ModeEst;
% macfin   = zeros(ModeEst2, N_tau);
% 
% for n = 1:N_tau
% %   uest = EigenFunctionEst(dTau(n));
% %   unum = NumCalculation(radius(n));
%     mac = MAC(u_num, u_ext);   
%     macfin(:,n) = reshape(mac', 1, ModeEst2);
% end
% 
% figure(3);
% plot(dTau, macfin(1:3,:), 'LineWidth', 2);
% legend('Uappr1, Unum1','Uappr1, Unum2','Uappr1, Unum3',3)
% grid on;
% axis([-inf inf -inf 1]);
% title(['Modal Assurance Criterium (Uappr1, Unum)']);
% xlabel('dTau');
% ylabel('log(MAC)');
% 
% figure(4);
% plot(dTau, macfin(4:6,:), 'LineWidth', 2);
% legend('Uappr2, Unum1','Uappr2, Unum2','Uappr2, Unum3',3)
% grid on;
% axis([-inf inf -inf 1]);
% title(['Modal Assurance Criterium (Uappr2, Unum)']);
% xlabel('dTau');
% ylabel('log(MAC)');
% 
% figure(5);
% plot(dTau, macfin(7:9,:), 'LineWidth', 2);
% legend('Uappr3, Unum1','Uappr3, Unum2','Uappr3, Unum3',3)
% grid on;
% axis([-inf inf -inf 1]);
% title(['Modal Assurance Criterium (Uappr3, Unum)']);
% xlabel('dTau');
% ylabel('log(MAC)');

