% this program calculate longitudinal wave propagation 
% propgating in a structure specified by a mesh file fMesh_Sub2Bond...
% First created by Prof. Haiying Huang at University of Texas Arlington
% on Jan. 25th, 2019
%Kranthi edits from May 21st 
close all;clear; clc;
tstart = tic;

% Default Material properties
%Hadhesive=300e-6;
%Segment=fMeshStep(Hadhesive)
Lbond=50e-3;
Hbond=100*1e-6;
Ebond=20e9;
%set up the segments of the whole body, define parameters 
Segment=fMesh_Sub2Bond(Lbond,Hbond, Ebond); %fMesh(Lbond,Hbond, Ebond);

% construct the structure
S=fSectionS(Segment);

% calculate displacement
f=300e3;  % tone-burst signal
[tu, Uxjj]=fTD_long(f,S,Segment);
Segment.tu=tu;
Segment.Uxjj=Uxjj;

% find amplitude at each location
if (length(Segment.xi(:,1))>1)  % multiple sections
    [Segment, Xmax, Umax]=fplotmax(Segment);
    %xlswrite('UxMax.xls', {'Xmax', 'Umax'}, '50um', 'A1'); 
    %xlswrite('UxMax.xls', [Xmax, Umax], '50um', 'A2'); 
else  % single location - for parameter studies
   figure(100); plot(tu*1e6, Uxjj); title('time-domain displacement'); xlabel('Time (us)'); ylabel('Displacement');
   [Segment, Xmax, Umax]=fUmax(Segment);
   % Umaxik(ii,kk)=Umax;
end
%figure(20); legend('structure', 'structure', 'Top plate'); xlabel('x position (m)');ylabel('Amplitude');
%figure(30); plot(tu*1e6,squeeze(Uxjj(1, :,3)), tu*1e6,squeeze(Uxjj(101, :,3)), 'r');
%legend('Left', 'Right'); xlabel('Time (us)');ylabel('Amplitude');

load Halleluya
sound(y1,Fs);
telapsed = toc(tstart)

%% save data
DataSave=1; %=input('Save data? Yes: 1, No: 0 \n');
if (DataSave==1)
    save data.mat Segment -v7.3
end

Ani=1; %input('Plot wave propagation animation? Yes: 1, No: 0 \n');
if (Ani==1)
   % fplot(Segment); % animation
   tsNew=1e-6;
   tend=200e-6;
   fplot_movie(Segment, tsNew, tend, min(Xmax), max(Xmax), max(Umax));
end

