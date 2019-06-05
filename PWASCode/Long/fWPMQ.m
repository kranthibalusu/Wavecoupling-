function [S, Mu,Ms]=fWPMQ(S, x, jj);
% this function calculates the displacement, slope, moment, and shear
% forces in the beam at x
% based on the solution stored in S
% Inputs:
%   S: structure stored information
%   x: coordinate
%   jj: section of the beam
% Outputs:
%   S: updated S with displacement
%   

N=length(S.betaI(jj, :));

 Mu=[];Ms=[];Mup=[];%Mp=[]; Mm=[]; 
for ii=1:N
    Mu=[Mu, exp(S.betaI(jj, ii)*x)];
   % Mup=[Mup, S.gammaI(jj,ii)*exp(S.betaI(jj, ii)*x)]; % 
%     Mp=[Mp, S.gammaI(jj, ii)*exp(S.betaI(jj, ii)*x)];
%     Mm=[Mm, S.muI(jj,ii)*exp(S.betaI(jj,ii)*x)];
    Ms=[Ms, S.KI(jj,ii)*exp(S.betaI(jj,ii)*x)];
end

for ii=1:N
    Mu=[Mu, exp(-S.betaI(jj, ii)*x)];
   % Mup=[Mup, S.gammaI(jj,ii)*exp(-S.betaI(jj, ii)*x)];
%     Mp=[Mp, -S.gammaI(jj,ii)*exp(-S.betaI(jj,ii)*x)];
%     Mm=[Mm, S.muI(jj,ii)*exp(-S.betaI(jj,ii)*x)];
    Ms=[Ms, -S.KI(jj,ii)*exp(-S.betaI(jj,ii)*x)];
end
S.u(:,jj)=Mu*S.adI(:,jj);
S.Mu=Mu;
%S.up(:,jj)=Mup*S.adI(:,jj);
% S.Phi(:, jj)=Mp*S.adI(:,jj);
% S.M(:,jj)=Mm*S.adI(:,jj);
S.Q(:,jj)=Ms*S.adI(:,jj);
return