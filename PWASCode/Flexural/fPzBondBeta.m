function [SP]=fPzBondBeta(SP, Omega)
% This function calculates the parameters for a piezo bonded to a beam of
% the same size as the piezo
% Inputs:
%   SP: structure for the PZ bonded section
%   Omega: angular velocity
% Outputs:
%   SP: updated SP

SP.C3=(SP.RhoS*Omega.^2*SP.IS*SP.hA/SP.GA/SP.hS-SP.ks*SP.GSc*SP.AS*SP.hA/SP.GA/SP.hS)/SP.alphaB;

% constants see derivation P5
SP.A=(SP.C1.*SP.RhoS.*Omega.^2/SP.ks/SP.GSc+SP.C2+SP.C3)./SP.C1-...
    (-SP.RhoP*Omega.^2+SP.GA/SP.hA/SP.hP)./SP.EPc;
SP.B=-SP.GA*SP.hS/2/SP.hA/SP.hP./(SP.EPc*SP.C1)-(-SP.RhoP.*Omega.^2+SP.GA/SP.hA/SP.hP).*...
        (SP.C1.*SP.RhoS.*Omega.^2/SP.ks/SP.GSc+SP.C2+SP.C3)./(SP.EPc*SP.C1)+...
        SP.RhoS*Omega.^2.*SP.C3/SP.ks/SP.GSc./SP.C1;
SP.C=-SP.RhoS*Omega.^2.*SP.C3/SP.ks/SP.GSc.*(-SP.RhoP.*Omega.^2+SP.GA/SP.hA/SP.hP)./(SP.EPc*SP.C1);

% solve for betaI see derivation P6
SP.betaI=sqrt(roots([1,SP.A,SP.B,SP.C])).';  % reduce to 3rd order differential Eq.

SP.KG2=SP.RhoS.*Omega.^2/SP.ks/SP.GSc;
SP.gammaI=SP.betaI+SP.KG2./SP.betaI;
SP.lambdaI=SP.C1.*SP.gammaI.*SP.betaI.^2+SP.C2.*SP.betaI+SP.C3.*SP.gammaI;
SP.kxiI=SP.lambdaI.*SP.betaI;
SP.itaI=SP.gammaI.*SP.betaI;
SP.chiI=SP.betaI-SP.gammaI;
SP.muI=-SP.ESc*SP.IS*SP.itaI;
SP.qI=-SP.ks*SP.AS*SP.GSc*SP.KG2./SP.betaI;
return