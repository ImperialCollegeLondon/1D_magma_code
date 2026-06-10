%% a trial of the new sill intrusion function file
clear all
clc
%% Inputs needed
depthN = 7;
depthCN = 7;
N = 20;
dz = [10*ones(1,5), 5*ones(1,10), 10*ones(1,5)];
dzF = 5;
nodez=zeros(1,N+1);
for i = 2:N+1
    nodez(i)=nodez(i-1)+dz(i-1);
end
cellz = (nodez(1:end-1)+nodez(2:end))/2;
Cb = ones(N,1);
H = ones(N,1);
T = ones(N,1);
C_all = [ones(N,1);zeros(N,1)];
phi = [ones(N,1);zeros(N,1)];
u_all = [ones(N+1,1);zeros(N+1,1)];
SillNodez = 10;
Sillcomp = 2;
SillH = 2;
SillTemp = 2; 
SillMF = 2;
SillCf = 2;
SillCm = 2;

[dz1, nodez1, N1, cellz1, Cb1, H1, T1, phi1, C_all1, Cphi_all1, u_all1] = sill_intrusion(N,depthN,depthCN,dz,dzF, Cb, H, T, C_all,phi, u_all, SillNodez, Sillcomp, SillH, SillTemp, SillMF, SillCf, SillCm);
