% this program calculate the dispersion relation of a bonded section
% longitudinal mode
% First created by Prof. Haiying Huang at University of Texas Arlington
% on April 5th, 2019
%
close all;clear; clc;

tstart = tic;
% Default Material properties
Lbond=1e-3; % no used
Hbond=100e-6;
Ebond=3e9; 
Segment=fSMeshBond(Lbond, Hbond, Ebond);

% construct the structure
S=fSectionS(Segment);

% Set up frequency sampling points
fStart=1e3;
fEnd=1000e3;
fk=[fStart:10e3:fEnd].';
    Amat = [];
    Bmat = [];
    ABmat =[];
% calculate Wavenumber 
for ii=1:length(fk)  % loop over non-zero frequency
    Omega=2*pi*fk(ii);  % angular frequency

%% parameters for different type of sections  
    [S,A,B]=fBeta(Segment, S, Omega);
    Beta(ii,:)=S(1).beta;
    Amat = [Amat,A];
    Bmat = [Bmat,B];
                      
end

%plot the dispersion curve
figure(10);plot(fk/1e3, imag(Beta(:,1)), fk/1e3, imag(Beta(:,2)), 'r-.', 'linewidth', 2); grid on;
title(['Ebond=', num2str(Ebond/1e9), ' GPa']);
xlabel('Freq (kHz)');ylabel('Beta-Imaginary');legend('beta_1', 'beta_2');

% plot 3D dispersion curve
figure(50); plot3(real(Beta(:,1)), imag(Beta(:,1)), fk/1e3, real(Beta(:,2)), imag(Beta(:,2)), fk/1e3, 'r','linewidth', 2);grid on;
xlabel('real(beta)');ylabel('imag(beta)');zlabel('Freq. (kHz)');legend('beta_1', 'beta_2');

% calculate group velocity
K=imag(Beta);  % wave number, only take the imaginary part
dK=diff(K);
df=diff(fk);
Cg(:,1)=2*pi*df./dK(:,1); % group velocity
Cg(:,2)=2*pi*df./dK(:,2); 
figure(20);plot(fk(2:end)/1e3, Cg(:,1), fk(2:end)/1e3, Cg(:,2), 'r-.', 'linewidth', 2); grid on;
title(['Group Velocity, Ebond=', num2str(Ebond/1e9), ' GPa']);
xlabel('Freq (kHz)');ylabel('Group Velocity (m/s)');legend('Cg_1', 'Cg_2');

% plot wavenumber
lambda0=sqrt(68.95e9/2715)./fk; % E/density 
figure(30); plot(fk/1e3, lambda0, 'k-. ',fk(2:end)/1e3, Cg(:,1)./fk(2:end), fk(2:end)/1e3, Cg(:,2)./fk(2:end), 'r-.', 'linewidth', 2); grid on;
title(['Wavelength, Ebond=', num2str(Ebond/1e9), ' GPa']);
xlabel('Freq (kHz)');ylabel('Wavelength (m)');legend('Al plate', 'lambda_1', 'lambda_2');
