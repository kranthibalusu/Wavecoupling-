function [S]=fBeam(S, SP)
% this function set up the material properties etc. for the beam section
% without piezo sensor
% Inputs:
%  SP: struct contains info for the section with piezo
%  S: struct holds the info for the beam section
% Outputs:
%  S: modified input

% beam info
S.h=SP.hS;  % height
S.b=SP.bS; % beam width should be the same as piezo width
S.E=SP.Es;  %beam Young's modulus; No loss in the beam. but in the piezo section
S.Rho=SP.RhoS; % beam density
S.ita=SP.itaS;
S.Ec=SP.ESc;  % complex stiffness, including losses
S.Gc=SP.GSc;
S.I=SP.IS;%% Not used in extensional mode!!
S.A=SP.AS;
return