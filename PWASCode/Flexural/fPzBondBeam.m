function [SP]=fPzBondBeam(SP, Gamma)
% this function set up the material properties etc. for the piezo
% bonded section
% Inputs:
%  SP: struct holds the info for the piezo bonded section
%  Gamma: shear lag parameter
% Outputs:
%  SP: modified input

%piezo geometry and material properties from Table 1 of Yan et al
SP.lP=10e-3; % piezo length
SP.aP=SP.lP/2;
SP.wP=10e-3; % piezo width
SP.hP=0.34e-3; % piezo thickness
SP.itaP=0.03;   % mechanical loss factor
SP.RhoP=7800;  % density
SP.d31P=-2.10e-10;  % charge constant
SP.e33P=2.14e-8;  % dielectric constant
SP.deltaP=0.0185;  % dielectric loss factor
SP.EP=9.1e10; %6.66e10;  % Young's modulus
SP.mu=0.3;    % Poisson's ratio, assumed for all materials

%construct complex parameters
SP.e33Pc=SP.e33P*(1-1i*SP.deltaP);
SP.EPc=SP.EP*(1+1i*SP.itaP);

%beam info
SP.hS=6e-3;  % beam thickness
SP.bS=10e-3; % beam width should be the same as piezo width
SP.Es=6.69e10;  %beam Young's modulus
SP.RhoS=2715; % beam density
SP.itaS=0.01;
SP.ESc=SP.Es*(1+1i*SP.itaS);  % complex stiffness, including losses
SP.GSc=SP.ESc/2/(1+SP.mu);
SP.IS=1/12*SP.hS^3*SP.bS;
SP.AS=SP.hS*SP.bS;

% glue information from Huang's paper
SP.EA=1.09e9;  % Young's modulus of glue
SP.GA=SP.EA/2/(1+SP.mu);  % shear modulus of adhesive
SP.ks=pi^2/12;  % shear correction factor; Mindlin plate theory
SP.thetaA=2; %0.17/0.2;  %Ha/hp
SP.hA=SP.thetaA*SP.hP;

% shear lag model
SP.psi=SP.ESc*SP.hS/SP.EPc/SP.hP;  % eq. (6) in Yan 2007
SP.Gbar=SP.GA/SP.EPc;
SP.ha_bar=SP.hA/SP.aP;
SP.Gamma=Gamma; %47.1;
% SP.Gamma= sqrt((SP.Gbar*SP.thetaA/(SP.ha_bar^2))*()/ ); % eq. (6) in Yan 2007

SP.alpha=SP.psi*(SP.ha_bar^2*SP.Gamma^2/SP.Gbar/SP.thetaA-1);
SP.alphaB=SP.alpha/2; %  divide by 2 to account for piezo on one side

% constants see derivation p3
SP.C1=SP.ESc*SP.IS*SP.hA/SP.GA/SP.hS/SP.alphaB;
SP.C2=(SP.ks*SP.GSc*SP.AS*SP.hA/SP.GA/SP.hS/SP.alphaB-SP.hS/2);
return
