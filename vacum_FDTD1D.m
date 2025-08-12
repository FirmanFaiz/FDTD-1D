close all; clc; clear all; 
%% Konstanta
mu0 = 2.013354451e-3;  % permeability of free space in x0.1(V fs^2/e nm)
ep0 = 55.26349597e-2;  % permittivity of free space in x0.1(e / V nm)
c0  = 29.9792458;      % light speed in x10(nm/fs)

% Parameter FDTD
dx = 2.65;      % 26,5 nm
Nx = 300;       % 3 mikrometer
dt = dx/(2*c0); % 0,0442 fs 
Nt = 300;       % total time step

x = (0:Nx-1)*dx;        % sumbu x
t = (0:Nt-1)*dt;        % time-step

% Elemen permitivitas dan permeabilitas (awal) untuk seluruh ruang
ep = ones(1,Nx)*ep0;
mu = ones(1,Nx)*mu0;

% Elemen Ex dan Hy (awal) untuk seluruh ruang
Hy = zeros(1,Nx);
Ez = zeros(1,Nx);
mHx = (dt)./mu; mEy = (dt)./ep; 

%% Looping 
for T = 1 : Nt
% Update H dari E
   for nx = 1 : Nx-1 
     Hy(nx) = Hy(nx) + mHx(nx)*(Ez(nx+1) - Ez(nx))/dx;
   end
    
% Update E dari H
   for nx = 2 : Nx
     Ez(nx) = Ez(nx) + mEy(nx)*(Hy(nx) - Hy(nx-1))/dx;
   end

plot (Ez,'-b');
axis ([0 Nx -1 1]);
xlabel ('x (x10 nm)','FontSize', 15);
ylabel ('E (V/nm)','FontSize', 15);
set(gca, 'FontSize', 13);
title(['Simulasi pada Iterasi ', num2str(T), ' (', sprintf('%.4f', t(T)), ' femtosekon)']);
%legend('Medan Listrik');
%set(Legend, 'Position', [0.2 0.7 0.2 0.1]); % Posisi legenda
getframe;
end