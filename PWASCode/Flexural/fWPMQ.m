function [S, Mw, Mp, Mm, Mq]=fWPMQ(S, x, jj);
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

Mw=[];Mp=[]; Mm=[]; Mq=[];
for ii=1:N
    Mw=[Mw, exp(S.betaI(jj, ii)*x)];
    Mp=[Mp, S.gammaI(jj, ii)*exp(S.betaI(jj, ii)*x)];
    Mm=[Mm, S.muI(jj,ii)*exp(S.betaI(jj,ii)*x)];
    Mq=[Mq, S.qI(jj,ii)*exp(S.betaI(jj,ii)*x)];
end

for ii=1:N
    Mw=[Mw, exp(-S.betaI(jj, ii)*x)];
    Mp=[Mp, -S.gammaI(jj,ii)*exp(-S.betaI(jj,ii)*x)];
    Mm=[Mm, S.muI(jj,ii)*exp(-S.betaI(jj,ii)*x)];
    Mq=[Mq, -S.qI(jj,ii)*exp(-S.betaI(jj,ii)*x)];
end
S.w(:,jj)=Mw*S.adI(:,jj);
S.Phi(:, jj)=Mp*S.adI(:,jj);
S.M(:,jj)=Mm*S.adI(:,jj);
S.Q(:,jj)=Mq*S.adI(:,jj);
return