function [zgp,wgp] = ModQuadTriangle(PtsTri,zgp_tri,wgp_tri)
% [zgp,wgp] = ModQuadTriangle(PtsTri,zgp_tri,wgp_tri)
% Integration points and weights in a triangle with vertices PtsTri
% zgp_tri, wgp_tri are the integration points and weights in the Reference element

v1 = [PtsTri(1,:)-PtsTri(3,:),0];
v2 = [PtsTri(2,:)-PtsTri(3,:),0];
AreaTri = norm(cross(v1,v2),2)/2;

zgp = reftri2tri(zgp_tri, PtsTri); 
wgp = wgp_tri*AreaTri/0.5; 



function Pts_x = reftri2tri(Pts_xi,PtsTri)
% reftri2tri (s,t,X,Y)

s = Pts_xi(:,1); t = Pts_xi(:,2); 
X = PtsTri(:,1); Y = PtsTri(:,2); 

X1 = X(1); X2 = X(2); X3 = X(3);
Y1 = Y(1); Y2 = Y(2); Y3 = Y(3);

a = X1-X3;
b = X2-X3;
c = X3;
x = a*s + b*t + c;

A = Y1-Y3;
B = Y2-Y3;
C = Y3;
y = A*s + B*t + C;

Pts_x = [x,y]; 