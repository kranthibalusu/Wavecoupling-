% this program produces the electric impedance of a piezowafer bonded
% on a beam with the same dimension as the piezo wafer using RMM
% follwing the derivation of Yan et al, 2007, An EMI approach for
% quantitative damage detection in Timoshenko beams with pz patches, Smart
% Mater. Struct., v16, p1390-1400
% First Created by Prof. Haiying Huang at University of Texas Arlington
% on Feb. 1st, 2013
% The program was modified by Md Mazharul Islam for the case of
% longitudinal mode of the beam and piezo
clear; clc;clear all;close all;
format long e;
% frequency
f=[10:0.1:500]*1e3; % Frequency in Hz
%f=171e3;
% define structures for the three sections
SP=struct;  % structure for the piezo bonded section
SL=struct;  % left section of the beam
SR=struct; % right section of the beam
%% SP.thetaA=0.4; is turned off in the file fPzBondBeam.m
% piezo bonded section
alpha=0.01;  % shear lag coefficient
[SP]=fPzBondBeam(SP, alpha);

%parameters unit for left and right portion
SL.l=9.5*SP.lP;  % length
SR.l=9.5*SP.lP;
SL.Ns=1; %  number of divisions
SR.Ns=1; 
SL.li=[1*SL.l];%, 0.2*SL.l];  % length of division starting at left
SR.li=[1*SR.l];%, 0.8*SL.l];
SL.ai=SL.li/2;
SR.ai=SR.li/2;

% Arbitrarily selected Timoshenho shear correction coeff for each division
SL.k=[0.1*SP.ks]; %, 0.5*SP.ks];  
SR.k=[0.1*SP.ks]; %, 0.1*SP.ks];
[SL]=fBeam(SL, SP);
[SR]=fBeam(SR, SP);
%% Coefficient of boundary condition
l=1.0; % displacement
m=1.0;% Force
n=1.0;% strain


for ii=1:length(f)
    Omega=2*pi*f(ii);  % angular frequency
    
    % set up matrix
    Nd=2*(SL.Ns+SR.Ns)+4; % # of unknows
    S=zeros(Nd, Nd);
    Pu=zeros(Nd, Nd);
    Q=zeros(Nd,1);
   
    % Setting S materix
    %% left side beam section and node(s)
    Srow=0;
    for jj=1:SL.Ns % loop over sections
        [SL]=fBeamBeta(SL, Omega, jj);
        
        
        if jj==1  % joint 0
            A=[m*SL.KI(jj,:)];
            B=[m*SL.KI(jj,:)];
%             A=[SL.muI(jj,:); SL.qI(jj,:)];
%             B=[-SL.muI(jj, :); SL.qI(jj,:)];


%         else % the following equation is not executed when jj=1.
%             %SL.gammaI(jj,:) is the solution of the beam only
%             A=[ones(1,4); SL.gammaI(jj-1,:), -SL.gammaI(jj,:);...
%                 SL.muI(jj-1,:), SL.muI(jj, :); ...
%                 SL.qI(jj-1,:), -SL.qI(jj-1,:)];
%             B=[-ones(1,4); SL.gammaI(jj-1,:), -SL.gammaI(jj,:);...
%                 -SL.muI(jj-1,:), -SL.muI(jj,:);...
%                 SL.qI(jj-1, :), -SL.qI(jj, :)];    
            
            
        end
        [Na, Ma]=size(A);
        S(Srow+[1:Na], Srow+[1:Na])=B\A;  %inv(B)*A
        %Q(Srow+[1:Na])=-B\[0];
        Srow=Srow+Na;
                
        % Phase matrix. The following equation uses one solution of betaI
        Pi=-exp(-SL.betaI(jj,:)*SL.li(jj));
        Pu((jj-1)*2+[1:2], (jj-1)*2+[1:2])=[zeros(1,1),diag(Pi);...
        diag(Pi), zeros(1,1)];

        
%         Pu((jj-1)*4+[1:4], (jj-1)*4+[1:4])=[zeros(2,2),diag(Pi);...
%            diag(Pi), zeros(2,2)];
    end
    

    
    %% left side of the piezo. jj is still 1.
    [SP]=fPzBondBeta(SP, Omega);
    %SL.gammaI is updated and it is the solution of the piezo bonded to
    %beam section only
    A=[ones(1,1),l*ones(1,2);...
       SL.KI(jj,:),-m*SP.KI;...
       zeros(1,1), SP.itaI];
    B=[-ones(1,1),-l*ones(1,2);...
       SL.KI(jj,:),-m*SP.KI;...
       zeros(1,1), SP.itaI];
    
    S(Srow+[1:3], Srow+[1:3])=B\A;
    Q(Srow+[1:3])=-B\[0;0;SP.d31P/SP.hP];
    Srow=Srow+3;
    
    % Phase matrix of the piezo section
    Pi=-exp(-SP.betaI*SP.lP);
    Pu(SL.Ns*2+[1:4], SL.Ns*2+[1:4])=[zeros(2,2),diag(Pi);...
           diag(Pi), zeros(2,2)];

    %% Right side of piezo
    jj=1;
    [SR]=fBeamBeta(SR, Omega, jj); % 1st division of rigth beam
    A=[ones(1,2),l*ones(1,1);...
       SP.KI,-m*SR.KI(jj,:);...
       SP.itaI,zeros(1,1),];
    B=[-ones(1,2),-l*ones(1,1);...
       SP.KI,-m*SR.KI(jj,:);...
       SP.itaI,zeros(1,1)];
    
%     A=[ones(1,5); SP.gammaI, -SR.gammaI(jj,:);...
%                 SP.muI, SR.muI(jj,:); ...
%                 SP.qI,-SR.qI(jj,:);
%                 -SP.kxiI,zeros(1,2)];
%     B=[-ones(1,5); SP.gammaI, -SR.gammaI(jj,:);...
%                 -SP.muI, -SR.muI(jj,:); ...
%                 SP.qI,-SR.qI(jj,:);
%                 SP.kxiI,zeros(1,2)];
    S(Srow+[1:3], Srow+[1:3])=B\A;
    Q(Srow+[1:3])=-B\[0;0;SP.d31P/SP.hP];
    Srow=Srow+3;

   % Phase matrix for right section of the beam
    Pi=-exp(-SR.betaI(jj,:)*SR.li(jj));
    Pu(SL.Ns*2+4+[1:2], SL.Ns*2+4+[1:2])=[zeros(1,1),diag(Pi);...
           diag(Pi), zeros(1,1)]; % the left end Pu is same as right end. It is derived considering as section not node.
    
    % right beam
%     for jj=2:SR.Ns % loop over sections. SR.NS=1. So it never be executed.
%         [SR]=fBeamBeta(SR, Omega, jj);
%         A=[ones(1,4); SR.gammaI(jj-1,:), -SR.gammaI(jj,:);...
%                 SR.muI(jj-1,:), SR.muI(jj, :); ...
%                 SR.qI(jj-1,:), -SR.qI(jj-1,:)];
%         B=[-ones(1,4); SR.gammaI(jj-1,:), -SR.gammaI(jj,:);...
%                 -SR.muI(jj-1,:), -SR.muI(jj,:);...
%                 SR.qI(jj-1, :), -SR.qI(jj, :)];                
%         [Na, Ma]=size(A);
%         S(Srow+[1:Na], Srow+[1:Na])=B\A;
%         Srow=Srow+Na;
%         
%         % Phase matrix
%         Pi=-exp(-SR.betaI(jj,:)*SR.li(jj));
%         Pu(SL.Ns*4+6+(jj-1)*4+[1:4], SL.Ns*4+6+(jj-1)*4+[1:4])=[zeros(2,2),diag(Pi);...
%            diag(Pi), zeros(2,2)];
%     end
%     
    %% last node 
    jj=1;
    [SR]=fBeamBeta(SR, Omega, jj);
    if jj==1  % joint 0
            A=[m*SL.KI(jj,:)];
            B=[m*SL.KI(jj,:)];
    end
    [Na, Ma]=size(A);
%     A=[SR.muI(SR.Ns,:); SR.qI(SR.Ns,:)];
%     B=[-SR.muI(SR.Ns, :); SR.qI(SR.Ns,:)];
    S(Srow+[1:Na], Srow+[1:Na])=B\A;
   
    %% solve for the constant d *a
    Dd=(eye(size(S))-S*Pu)\Q;
    Aa=Pu*Dd;
    
    % set adI={a,d} for left coordinates
    for jj=1:SL.Ns
        SL.adI(:,jj)=[Aa(2*(jj-1)+[1:1]);Dd(4*(jj-1)+[1:1])];
    end  
    SP.adI=[Aa(2*SL.Ns+[1:2]);Dd(2*SL.Ns+[1:2])];
    
    for jj=1:SR.Ns
        SR.adI(:,jj)=[Aa(SL.Ns*2+4+2*(jj-1)+[1:1]);Dd(SL.Ns*2+4+2*(jj-1)+[1:1])];
    end

    % Calculate the impedance see derivation P8
    Eu=[SP.gammaI, SP.gammaI].*SP.adI.';% Page 8
    Upl=Eu*[exp(SP.betaI*0), exp(-SP.betaI*0)].';  % Up(0);
    Upr=Eu*[exp(SP.betaI*SP.lP), exp(-SP.betaI*SP.lP)].';  % Up(lP);

    Y(ii)=1i*Omega*SP.wP*SP.lP*(SP.e33Pc-SP.d31P^2*SP.EPc)/SP.hP+...
        1i*Omega*SP.wP*SP.d31P*SP.EPc*(Upr-Upl);
    %Iuf= 1i*Omega*SP.wP*SP.d31P*SP.EPc*(Upr-Upl);
    Iuf=1i*Omega*SP.wP*SP.lP*(SP.e33Pc-SP.d31P^2*SP.EPc)/SP.hP;
end

figure(11);plot(f,abs(Y));
%return

if (length(f)==1)  %calculate displacement, moment etc
    figure(1);
    % Left section
    x0=0
    for jj=1:SL.Ns
        Li=SL.li(jj)
        x=[0:Li/100:Li].';
        [SL]=fWPMQ(SL,x,jj);
        subplot(2,1,1); plot(x+x0,real(SL.u(:,jj)), 'r');title('displacement'); grid on;hold on;
%         subplot(2,2,2); plot(x+x0,real(SL.Phi(:, jj)), 'r');title('Slope'); grid on;hold on;
%         subplot(2,2,3); plot(x+x0,real(SL.M(:, jj)), 'r');title('Moment'); grid on;hold on;
        subplot(2,1,2); plot(x+x0,real(SL.Q(:, jj)), 'r');title('Normal stress'); grid on;hold on
        x0=x0+Li;
    end
            
    % Piezo section
    Li=SP.lP;
    x=[0:Li/100:Li].';
    [SP]=fWPMQ(SP, x, 1);
    subplot(2,1,1); plot(x+x0,real(SP.u), 'b');title('displacement'); grid on;
%     subplot(2,2,2); plot(x+x0,real(SP.Phi), 'b');title('Slope'); grid on;
%     subplot(2,2,3); plot(x+x0,real(SP.M), 'b');title('Moment'); grid on;
    subplot(2,1,2); plot(x+x0,real(SP.Q), 'b');title('Normal stress'); grid on;
    x0=x0+Li
    
    %adhesive shear stress
    for xii=1:length(x)
        Tau(xii)=SP.GA/SP.hA*(SP.adI(1:2).'.*(SP.gammaI+1)*exp(SP.betaI.'*x(xii))+...
            SP.adI(3:4).'.*(SP.gammaI+1)*exp(-SP.betaI.'*x(xii)));
    end
    figure(2); subplot(2,1,1);plot(x, real(Tau)); title('Shear stress in adhesive layer'); grid on; ylabel('Real(Tau)');
    subplot(2,1,2); plot(x, imag(Tau)); xlabel('x-coord'); ylabel('imag(Tau)');grid on;
    
    figure(1);
    % right section
    for jj=1:SR.Ns
        Li=SR.li(jj);
        x=[0:Li/100:Li].';
        [SR]=fWPMQ(SR,x,jj);
        subplot(2,1,1); plot(x+x0,real(SR.u(:,jj)),'r');title('displacement'); grid on;
%         subplot(2,2,2); plot(x+x0,real(SR.Phi(:,jj)),'r');title('Slope'); grid on;
%         subplot(2,2,3); plot(x+x0,real(SR.M(:,jj)),'r');title('Moment'); grid on;
        subplot(2,1,2); plot(x+x0,real(SR.Q(:,jj)),'r');title('Normal stress'); grid on;
        x0=x0+Li
    end
    
else
     
    
    figure(1); subplot(2,1,1); plot(f/1e3, real(Y),'r'); xlabel('Freq (kHz)'); ylabel('Real (Y)'); grid on;
    subplot(2,1,2); plot(f/1e3, imag(Y),'r'); xlabel('Freq (kHz)'); ylabel('imag(Y)'); grid on;
    [SP.thetaA, max(real(Y)), f(real(Y)==max(real(Y)))];
     
    % calculate time domain signal from Y
    [t, It]=f_ifft_t(Y, f(2)-f(1), max(length(Y), 2^16));
    figure(5); plot(t, real(It)); xlabel('Time (s)'); ylabel('It calculated from Y');
    
    % calculated spectrogram
    WinLen=round(length(It)/100);
    Overlap=floor(WinLen/2);
    figure (10); spectrogram(It, WinLen, Overlap, [10:1:500]*1e3, 1/(t(2)-t(1)));
    %fprintf(1, '%10.4e \t %10.4e \n', [real(Y);imag(Y)]);
end

data = abs (Y)'
fid=fopen('Longitudinal.txt','w+');
fprintf(fid,'% Freq (Hz) \t abs (Y)\n');
fprintf(fid,'%6.4f \t %10.5e \t %10.5e\n',[f;abs(Y);real(Y)]);
fclose(fid);

 SL.C=sqrt(6.895E10/2715);
 1*SL.C./(2*(SL.l*2+SP.lP))
% 
