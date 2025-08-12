close all; clc; clear all; 

%% Konstanta
mu0 = 2.013354451e-3;   % permeability of free space in x0.1(V fs^2/e nm)
ep0 = 55.26349597e-2;   % permittivity of free space in x0.1(e / V nm)
c0  = 29.9792458;       % light speed in x10(nm/fs)

% Parameter FDTD
dx = 2.65;              % 25,9 nm
Nx = 300;               % 2,1 mikrometer
dt = dx/(2*c0);         % 0,0442 fs
Nt = 1700;              % total time step

% Elemen permitivitas dan permeabilitas (awal) untuk seluruh ruang
ep = ones(1,Nx)*ep0;
mu = ones(1,Nx)*mu0;
er = 13;
% 2.1; 2.4; 2.6; 3.0; 3.3; 3.5; 3.9; 7.0; 7.4; 9.8

% Medium dielektrik
d1 = Nx/2;  %titik awal dielektrik
d2 = Nx;    %tiitk akhir dielektrik
ep(d1:d2) = er*ep0;

disp(ep); disp(mu);
mHx = (dt)./mu; mEy = (dt)./ep; 

x = (0:Nx-1)*dx;        % sumbu x
t = (0:Nt-1)*dt;        % time-step

% Elemen Ex dan Hy (awal) untuk seluruh ruang
Hy = zeros(1,Nx);
Ez = zeros(1,Nx);

%% Parameter Gelombang
nxsrc = round(Nx/4);    % titik sumber
nbc = 1;                % Indeks refraktif ruang (awal)
tau = (sqrt(2))*10*dt;  % lebar pulsa
t0 = 4*tau;

disp (dt)
disp (tau)

% Formulasi batas TFSF
A = -sqrt(ep(nxsrc)/mu(nxsrc));     % Amplitudo untuk H
st = ((nbc*dx)/(2*c0))+ dt/2;       % 1,5 karena H iterasi +0,5. jadi iterasi selanjutnya adalah 1,5
Esrc = exp(-((t-t0)/tau).^2);
Hsrc = A*exp(-((t-t0+st)/tau).^2);

%% Formulasi PML
sigma = zeros(1, Nx); % initialize conductivity array
d = 20*dx; % jumlah grid untuk lapisan PML
m = 2; % polynomial order for grading sigma array
R = 1e-5; % required reflectivity

% Left PML (untuk area 1 sampai d)
for nz = 1:d
    sigma_max = (-(m+1)) * log10(R) * ep0 * c0 / (2 * d);
    sigma(nz) = sigma_max * ((1 - nz/d).^m);
end

% Right PML (untuk area Nx-d sampai Nx)
for nz = round(Nx-d+1):Nx % Use round() or floor() to ensure integer indices
    sigma_max = (-(m+1)) * log10(R) * ep0*er * c0 / (2 * d);
    sigma(nz) = sigma_max * (((nz - (Nx-d))/d).^m);
end

% Recalculate PML constants for each point in the grid
sigma_star = sigma .* mu ./ ep; % Update sigma_star based on local ep

% Konstanta Update PML
A=((ep-0.5*dt*sigma)./(ep+0.5*dt*sigma));
B=(dt/dx)./(ep+0.5*dt*sigma);
C=((mu-0.5*dt*sigma_star)./(mu+0.5*dt*sigma_star)); 
D=(dt/dx)./(mu+0.5*dt*sigma_star);                          

%% Looping
for T = 1 : Nt    
% Update H dari E
   for nx = 1 : Nx-1 
     Hy(nx) = C(nx)*Hy(nx) + D(nx)*(Ez(nx+1) - Ez(nx));
   end
    
%H Field Source
    Hy(nxsrc-1) = Hy(nxsrc-1) - D(nxsrc-1)*Esrc(T);

% Update E dari H
   for nx = 2 : Nx
    Ez(nx) = A(nx)*Ez(nx) + B(nx)*(Hy(nx) - Hy(nx-1));
   end
 
%E Field Source  
    Ez(nxsrc) = Ez(nxsrc) - B(nxsrc)*Hsrc(T);

h1 = plot(Ez,'-k','LineWidth',2);
axis([0 Nx -1.5 1.5]);
hold on
batas_1=line([d d], [-1.5 1.5], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2); % Garis vertikal kiri
line([Nx-d Nx-d], [-1.5 1.5], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2); % Garis vertikal kanan
tfsf = line([nxsrc nxsrc], [-1.5 1.5], 'Color', 'b', 'LineStyle', '--', 'LineWidth', 1); % Garis vertikal kiri
dielektrik=line([round(Nx/2) round(Nx/2)], [-1.5 1.5], 'Color', 'g', 'LineStyle', '--', 'LineWidth', 2); % Garis vertikal kiri
%legend([h1,tfsf,dielektrik,batas_1],'Medan Listrik','Batas TFSF','Batas Dielektrik','Batas PML','FontSize', 20, 'Location', 'southeast');
xlabel('x (x10 nm)', 'FontSize', 25);
ylabel('E (V/nm)', 'FontSize', 25);
set(gca, 'FontSize', 25);
title(['Simulasi pada Iterasi ', num2str(T), ' (', sprintf('%.4f', t(T)), ' femtosekon)']);
hold off
getframe;

end