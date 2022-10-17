clear all; close all; clc


L1 = 0.2;
L2 = 0.06;
L = L1;
[w_1, x_1] = solid_1(L1);
[w_2, x_2] = solid_1(L2);
plot(x_1,w_1)
hold on
plot(x_2,w_2)
function [w, x] = solid_1(L)
E = 200; %[GPa] or Pa? hooke does not specify
v = 0.3; %[poisson]
F = 2000; %[N]
h = 0.02; % [m]
b = h;
p1 = [0,0];
p2 = [L, h];
nelx = 100;
nely = 5;
[D]=hooke(1,E,v);
ndofs = 2;
[Edof,Ex,Ey,B1,B2,B3,B4]=quad_mesh(p1,p2,nelx,nely,ndofs);
n = length(Edof);
%eldraw2(Ex,Ey) % draws mesh 
ep = [1,1,2];

K = zeros(max(max(Edof)));
f = zeros(max(max(Edof)),1);
for i = 1:n
    ex = Ex(i,:);
    ey = Ey(i,:);
    
    if i == nelx
        eq = [F; 0];
        [Ke,fe] = plani4e(ex,ey,ep,D,eq);
    else
        eq = [0; 0];
        [Ke,fe] = plani4e(ex,ey,ep,D,eq);
    end
    [K,f] = assem(Edof(i,:),K,Ke,f,fe);
end
B4_nodes = [B4(:,1)' , B4(:,2)']';
bc = [B4_nodes, zeros(length(B4_nodes), 1)];
aa = solveq(K,f,bc);
ed = extract(Edof, aa);
[sfac] = eldisp2(Ex,Ey,ed);
w = aa(4:2:202);
x = linspace(L/nelx,L,100);
for i = 1:n
    [es(i,:),~] = plani4s(Ex(i,:),Ey(i,:),ep,D,ed(i,:));
end
end