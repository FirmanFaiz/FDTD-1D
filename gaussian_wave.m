close all; clc; clear all; 
%% Konstanta
mu0 = 2.013354451e-3;  % permeability of free space in x0.1(V fs^2/e nm)
ep0 = 55.26349597e-2;  % permittivity of free space in x0.1(e / V nm)
c0  = 29.9792458;      % light speed in x10(nm/fs)

% Parameter FDTD
dx = 2.65;      % 15,6 nm
Nx = 300;       % 2 mikrometer
dt = dx/(2*c0); % 0,0432 fs
Nt = 800;       % total time step

% Elemen permitivitas dan permeabilitas (awal) untuk seluruh ruang
ep = ones(1,Nx)*ep0;
mu = ones(1,Nx)*mu0;

mHx = (dt)./mu; mEy = (dt)./ep; 

x = (0:Nx-1)*dx; %sumbu x
t = (0:Nt-1)*dt; %time-step

% Elemen Ex dan Hy (awal) untuk seluruh ruang
Hy = zeros(1,Nx);
Ez = zeros(1,Nx);

%% Parameter Gelombang
nxsrc = round(Nx/2);    % titik sumber
nbc=1;                  % Indeks refraktif ruang (awal)
tau = (sqrt(2))*10*dt;  % lebar pulsa
t0 = 4*tau;

src = exp(-((t-t0)/tau).^2); %Persamaan gelombang

%% Looping 
for T = 1 : Nt
    % Update H from E
    for i = 1 : Nx-1 
        Hy(i) = Hy(i) + mHx(i)*(Ez(i+1) - Ez(i))/dx;
    end

    % Update E from H
    for i = 2 : Nx
        Ez(i) = Ez(i) + mEy(i)*(Hy(i) - Hy(i-1))/dx;
    end

    % Wave source for Ez
    Ez(nxsrc) = src(T);

plot(Ez, '-k','LineWidth', 2);
axis([0 Nx -1.5 1.5]);
xlabel('x (x10 nm)', 'FontSize', 25);
ylabel('E (V/nm)', 'FontSize', 25);
set(gca, 'FontSize', 25);
title(['Simulasi pada Iterasi ', num2str(T), ' (', sprintf('%.4f', t(T)), ' femtosekon)']);
% L=legend('Medan Listrik', 'Location', 'southeast', 'FontSize', 20);
% L.Position = [0.65, 0.25, 0.3, 0.1]; % [x, y, width, height] dalam normalisasi
getframe;

end