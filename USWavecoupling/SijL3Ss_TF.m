% this program calculates the S-parameter of a bonded section
% longitudinal mode
% First created by Prof. Haiying Huang at University of Texas Arlington
% on Jan. 25th, 2019
%
close all;clear; clc;
tstart = tic;

% Default Material properties
Lbond=2e-3;
Hbond=100e-6;
Ebond=0.3e9; 
xi=[0,1]; % location of S-parameter xi=x/L
jSect=3; % section at which the calculating is done
ToB=1; %for bonding section only, 1=top 0=bottom
%Segment=fSjMesh1SSa(xi, jSect, ToB, Lbond, Hbond, Ebond); %fSMeshSsBondSs(); %fMesh2SsF();
Segment=fSjMeshSBS(xi, jSect, ToB, Lbond,Hbond, Ebond);

% construct the structure
S=fSectionS(Segment);

% Set up frequency sampling points
fStart=50e3;
fEnd=1000e3;
[t, Vt, fk, Vfc, index]=fToneBurst(800e3);   % for getting the frequency points only
index1=find(fk>=fStart &fk<=fEnd);
fkS=fk(index1);  % select freq for S-parameter calculation

% calculate & plot dispersion curve of the bonding section
fBond_dispersion(fkS,Segment, S);

% calculate S-parameter
Segment.Sij=fSijL(fkS, S, Segment);

load Halleluya
sound(y1,Fs);


%% time-frequency analysis
fii=[50e3:10e3:800e3].';  % frequency for tone-burst signal
    
figure(100); xlabel('Freq. (kHz)');ylabel('Transmission |S|'); grid on; hold on
title(['Lbond=', num2str(Lbond*1e3), 'mm, Ebond=', num2str(Ebond/1e9), 'GPa']);
LegendTxt={};
for ii=1:length(xi)
    Sii=squeeze(Segment.Sij(ii,:,jSect));

    % plot the S-parameter
    figure(100); plot(fkS/1e3, abs(Sii), 'linewidth', 2);

    [ttf, SiiTF]=fSij_TF(fkS.', Sii.', 5.5, fii);
    LegendTxt(ii)={['xi=', num2str(xi(ii))]};
    
%     %plot the time-frequency signal
%     if(length(fii)==1)
%         figure(100+ii);plot(ttf*1e6, SiiTF,  'linewidth', 2); xlabel('Time (us)'); ylabel(['u_4(t), f=',num2str(fii/1e3), ' kHz']); title(['jSect=', num2str(jSect), ', xi=', num2str(xi(ii))]);
%         axis([0 200 min(SiiTF)*1.1 max(SiiTF)*1.1]);
%     else
%         figure(100+ii); mesh(fii/1e3, ttf*1e6, SiiTF); xlabel('Freq (kHz)'); ylabel('Time (us)');
%         axis([min(fii/1e3), max(fii/1e3), 0, 500]); view(102,32);
%         title(['jSect=', num2str(jSect), ', xi=', num2str(xi(ii))]);
%     end
end
figure(100); legend(LegendTxt); hold off;
