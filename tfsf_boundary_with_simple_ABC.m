close all; clc; clear all; 

%% Konstanta
mu0 = 2.013354451e-3;   % permeability of free space in x0.1(V fs^2/e nm)
ep0 = 55.26349597e-2;   % permittivity of free space in x0.1(e / V nm)
c0  = 29.9792458;       %light speed in x10(nm/fs)

% Parameter FDTD
dx = 2.65;           % 26,5 nm
Nx = 300;            % 3 mikrometer
dt = dx/(2*c0);      % 0,0442 fs
Nt = 650;            % total time step

% Elemen permitivitas dan permeabilitas (awal) untuk seluruh ruang
ep = ones(1,Nx)*ep0;
mu = ones(1,Nx)*mu0;

x = (0:Nx-1)*dx;        % sumbu x
t = (0:Nt-1)*dt;        % time-step

% Elemen Ex dan Hy (awal) untuk seluruh ruang
Hy = zeros(1,Nx);
Ez = zeros(1,Nx);
mHx = (dt)./mu; mEy = (dt)./ep; 

%% Parameter Gelombang
nxsrc = round(Nx/4);    % titik sumber
nbc=1;                  % Indeks refraktif ruang (awal)
tau = (sqrt(2))*10*dt;  % lebar pulsa
t0 = 4*tau;

% Formulasi batas TFSF
A = -sqrt(ep(nxsrc)/mu(nxsrc));     % Amplitudo untuk H
st = ((nbc*dx)/(2*c0))+ dt/2;       % 1,5 karena H iterasi +0,5. jadi iterasi selanjutnya adalah 1,5
Esrc = exp(-((t-t0)/tau).^2);
Hsrc = A*exp(-((t-t0+st)/tau).^2);

%% Nilai tetapan batas Dirichlet

E_low_m2=0;
E_low_m1=0;
E_high_m2=0;
E_high_m1=0;

%% Looping
for T = 1 : Nt
% Update H dari E
   for nx = 1 : Nx-1 
     Hy(nx) = Hy(nx) + mHx(nx)*(Ez(nx+1) - Ez(nx))/dx;
   end
    
% H Field Source
    Hy(nxsrc-1) = Hy(nxsrc-1) - mHx(nxsrc-1)*Esrc(T)/dx;


% Update E dari H
   for nx = 2 : Nx
    Ez(nx) = Ez(nx) + mEy(nx)*(Hy(nx) - Hy(nx-1))/dx;
   end
 
% Boundary Condition
    Ez(1) =E_low_m2;
    E_low_m2=E_low_m1;
    E_low_m1= Ez(2);
    Ez(Nx)= E_high_m2;
    E_high_m2 = E_high_m1;
    E_high_m1 = Ez(Nx - 1);
    
%E Field Source  
    Ez(nxsrc) = Ez(nxsrc) - mEy(nxsrc)*Hsrc(T)/dx;

plot(Ez,'-k');
axis([0 Nx -1.5 1.5]);
hold on
tfsf = line([nxsrc nxsrc], [-1.5 1.5], 'Color', 'b', 'LineStyle', '--', 'LineWidth', 1); % Garis vertikal kiri
%legend([h1,tfsf],'Medan Listrik','Batas TFSF','FontSize', 20, 'Location', 'southeast');
xlabel('x (x10 nm)', 'FontSize', 25);
ylabel('E (V/nm)', 'FontSize', 25);
set(gca, 'FontSize', 25);
title(['Simulasi pada Iterasi ', num2str(T), ' (', sprintf('%.4f', t(T)), ' femtosekon)']);
hold off
getframe;

end