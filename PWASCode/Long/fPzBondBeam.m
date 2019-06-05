function [SP]=fPzBondBeam(SP, alpha)
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
SP.hP=0.2e-3; % piezo thickness
SP.itaP=0.03;   % mechanical loss factor
SP.RhoP=7800;  % density
SP.d31P=-2.10e-10;  % charge constant
SP.e33P=2.14e-8;  % dielectric constant
SP.deltaP=0.0185;  % dielectric loss factor
SP.EP=7.1e10; %6.66e10;  % Young's modulus
SP.muP=0.3;    % Poisson's ratio, assumed for all materials
SP.muB=0.3;% Beam
SP.muA=0.3; % Adhesive
%construct complex parameters
SP.e33Pc=SP.e33P*(1-1i*SP.deltaP);
SP.EPc=SP.EP*(1+1i*SP.itaP);

%beam info
SP.hS=1.52e-3;  % beam thickness
SP.bS=10e-3; % beam width should be the same as piezo width
SP.Es=6.895e10;  %beam Young's modulus
SP.RhoS=2715; % beam density
SP.itaS=0.015;
SP.ESc=SP.Es*(1+1i*SP.itaS);  % complex stiffness, including losses
SP.GSc=SP.ESc/2/(1+SP.muB);%% was it used in the program ever?????????;see fBeamBeta.m
SP.IS=1/12*SP.hS^3*SP.bS;% useless for extensional mode!!1
SP.AS=SP.hS*SP.bS;

% glue information from Huang's paper
SP.GA=0.6666e9;  % shear modulus of adhesive


SP.ks=pi^2/12;  % shear correction factor; Minline plate theory
SP.thetaA=1.0; %0.17/0.2;  %Ha/hp
SP.hA=SP.thetaA*SP.hP;

% shear lag model
SP.psi=SP.Es*SP.hS/SP.EP/SP.hP;  % eq. (6)
SP.Gbar=SP.GA/SP.EP;
SP.ha_bar=SP.hA/SP.aP;
SP.Gamma= sqrt((SP.Gbar*SP.thetaA/(SP.ha_bar^2))*(SP.psi+alpha)/SP.psi) %47.1 ideally;
SP.alpha=alpha;%SP.psi*(SP.ha_bar^2*SP.Gamma^2/SP.Gbar/SP.thetaA-1);
SP.alphaB=SP.alpha/2; %  divide by 2 to account for piezo on one side

% constants see derivation p3
SP.C1=-SP.ESc*SP.AS*SP.hA/SP.GA/SP.alphaB/SP.bS;
SP.C2=0;
%SP.C3=1-SP.RhoS*SP.AS*SP.hA*(2*pi*f).^2/SP.GA/SP.alphaB;% No f incoming so uncomment this line
%%SP.C1=SP.ESc*SP.IS*SP.hA/SP.GA/SP.hS/SP.alphaB; %flexure mode
%%SP.C2=(SP.ks*SP.GSc*SP.AS*SP.hA/SP.GA/SP.hS/SP.alphaB-SP.hS/2);%% flexure mode
return
