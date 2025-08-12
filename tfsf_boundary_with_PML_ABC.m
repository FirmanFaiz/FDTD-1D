close all; clc; clear all; 

%% Konstanta
mu0 = 2.013354451e-3;   % permeability of free space in x0.1(V fs^2/e nm)
ep0 = 55.26349597e-2;   % permittivity of free space in x0.1(e / V nm)
c0  = 29.9792458;       %light speed in x10(nm/fs)

% Parameter FDTD
dx = 2.65;              % 26,5 nm
Nx = 300;               % 2,1 mikrometer
dt = dx/(2*c0);         % 0,0442 fs
Nt = 1000;               % total time step

% Elemen permitivitas dan permeabilitas (awal) untuk seluruh ruang
ep = ones(1,Nx)*ep0;
mu = ones(1,Nx)*mu0;

disp(ep); disp(mu);
mHx = (dt)./mu; mEy = (dt)./ep; 

x = (0:Nx-1)*dx;        % sumbu x
t = (0:Nt-1)*dt;        % time-step

% Elemen Ex dan Hy (awal) untuk seluruh ruang
Hy = zeros(1,Nx);
Ez = zeros(1,Nx);

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

%% Formulasi PML
sigma(1:Nx)=0; % initialize conductivity array
d=20*dx; % width of PML layer
m=2; % polynomial order for grading sigma array (1-5)
neta=sqrt(mu0/ep0); 
R=1e-5; % required reflectivity
sigma_max=((-1)*(m+1)*log(R))/(2*neta*d*dx);
rho=((1:d+1)./d).^m*sigma_max;
sigma(Nx-d:Nx)=rho; % lossy conductivity profile
sigma(1:d+1)=fliplr(rho);
sigma_star(1:Nx)=sigma.*mu./ep;

disp(sigma_max)

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
%legend([h1,tfsf,batas_1],'Medan Listrik','Batas TFSF','Batas PML','FontSize', 20, 'Location', 'southeast');
xlabel('x (x10 nm)', 'FontSize', 25);
ylabel('E (V/nm)', 'FontSize', 25);
set(gca, 'FontSize', 25);
title(['Simulasi pada Iterasi ', num2str(T), ' (', sprintf('%.4f', t(T)), ' femtosekon)']);
hold off
getframe;

end