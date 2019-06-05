function [S]=fBeamBeta(S, Omega, Di)
% This function calculates the dynamic parameters for the beam sections 
% Inputs:
%   S: structure holding info for the beam section
%   Omega: angular velocity
%   Di: section number
% Outputs:
%   S: updated S

A=S.Ec*S.I;
B=S.Rho*Omega.^2*S.I*(1+S.Ec/S.k(Di)/S.Gc);
C=S.Rho*Omega^2*S.I/S.k(Di)/S.Gc-S.Rho*Omega^2*A;

S.betaI(Di, :)=sqrt(roots([A,B,C])).';  % row vector

S.KG2(Di)=S.Rho.*Omega^2/S.k(Di)/S.Gc;
S.gammaI(Di, :)=S.betaI(Di, :)+S.KG2(Di)./S.betaI(Di, :);
S.muI(Di, :)=-S.Ec*S.I*S.gammaI(Di, :).*S.betaI(Di, :);
S.qI(Di, :)=-S.k(Di)*S.A*S.Gc.*S.KG2(Di)./S.betaI(Di, :);
return


