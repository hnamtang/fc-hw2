% Flight Control - Group 2 - 2nd Homework
% Task 2
clearvars; close all; clc

% Import state-space model of Talon UAV (source: AlphaLink Engineering)
run('vfteStateSpace');
sys=G;

%% Setup filters
% Connect low pass filter
T_LP=.1;
A_lp=-1/T_LP*eye(2);
B_lp=zeros(2,size(sys.C,1));B_lp(1,5)=1/T_LP;B_lp(2,7)=1/T_LP;
C_lp=zeros(size(sys.C,1),2);C_lp(5,1)=1;C_lp(7,2)=1;
D_lp=eye(size(sys.C,1));D_lp(5,5)=0;D_lp(7,7)=0;
lp=ss(A_lp,B_lp,C_lp,D_lp);
lp.InputName=sys.OutputName;
lp.OutputName=sys.OutputName;
lp.OutputName([5 7])={'r_{lp}','p_{lp}'};
lp.StateName={'x_{rlp}','x_{plp}'};

sys1=series(sys,lp); % plant + LP

% Connect high-pass (washout) filter
T_HP=2;
A_wash=-1/T_HP;
B_wash=zeros(1,size(sys.C,1));B_wash(1,5)=-1/T_HP;
C_wash=zeros(size(sys.C,1),1);C_wash(5,1)=1;
D_wash=eye(size(sys.C,1));
wash=ss(A_wash,B_wash,C_wash,D_wash);
wash.InputName=lp.OutputName;
wash.OutputName=lp.OutputName;
wash.OutputName(5)={'r_{wash}'};
wash.StateName={'x_{wash}'};

sys2=series(sys1,wash); % plant + LP + HP

[A,B,C,D]=ssdata(sys2);


%% Design Yaw Damper
% Determine controller gain k_zetar using root locus technique
[A_r,B_r,C_r,D_r]=ssdata(sys2(5,4));
[z_r,p_r,k_r]=ss2zp(A_r,B_r,C_r,D_r);

% k=linspace(0,5,10000);
% r=rlocus(A_r,B_r,-C_r,D_r,k);
% figure; hold on
% plot(r,'LineWidth',1.5);
% plot(real(p_r),imag(p_r),'kx','MarkerSize',9);
% plot(real(z_r),imag(z_r),'ko','MarkerSize',9);
% sgrid
% hold off

%k_zetar=2; % from root locus -> DR damping ratio of 0.387
k_zetar=.07; % from root locus -> DR damping ratio of 0.284

% Close yaw rate loop
Acl_r=A+B(:,4)*k_zetar*C(5,:); % positive feedback


%% Design Roll Damper
% Determine controller gain k_xip using root locus technique
[z_p,p_p,k_p]=ss2zp(Acl_r,B(:,3),C(7,:),D(7,3));

% k=linspace(0,.06,10000);
% r=rlocus(Acl_r,B(:,3),C(7,:),D(7,3),k);
% figure; hold on
% plot(r,'LineWidth',1.5);
% plot(real(p_p),imag(p_p),'kx','MarkerSize',9);
% plot(real(z_p),imag(z_p),'ko','MarkerSize',9);
% sgrid;
% hold off

k_xip=.024;
 
% Close the roll rate loop (with yaw rate loop) -> SAS
Kcl_SAS=zeros(size(B,2),size(C,1));
Kcl_SAS(4,5)=k_zetar;Kcl_SAS(3,7)=k_xip; % positive feedback
Acl_SAS=A+B*Kcl_SAS*C;
sys_SAS=ss(Acl_SAS,B,C,D);
sys_SAS.InputName=sys2.InputName;
sys_SAS.OutputName=sys2.OutputName;
sys_SAS.StateName=sys2.StateName;

[A,B,C,D]=ssdata(sys_SAS);

%plotResponse(sys([5 6 7 8],3:4),sys_SAS([5 6 7 8],3:4),10,'doublet');


%% Design Roll Angle Control (Phi feedback loop)
% Only k_phiP
H=zeros(1,size(A,2));H(8)=1;
F=0;G=0;D=0;J=1;

% Augmented dynamics with compensator
Aa_phi=[A zeros(size(A,1),1);-G*H F];
Ba_phi=[B(:,3);0];
Ga_phi=[zeros(size(B,1),1);G];
%Ca_phi=[C(8,:) zeros(8,1);-J*H D];
Ca_phi=[-J*H D];
Fa_phi=J;
Ha_phi=[H 0];

% figure;rlocus(Aa_phi,Ba_phi,Ca_phi,Fa_phi);set(findall(gca,'Type','Line'),'LineWidth',1.5,'MarkerSize',9);

%k_phiP=.07;
k_phiP=.1;

% Close the phi loop
Aac_phi=Aa_phi-Ba_phi*k_phiP*Ca_phi;
Bac_phi=Ga_phi-Ba_phi*k_phiP*Fa_phi;
sysac_phi=ss(Aac_phi,Bac_phi,Ha_phi,0);
sysac_phi.InputName={'e_{Phi}'};
sysac_phi.OutputName={'Phi'};

% figure; margin(sysac_phi);

% Step response
t=0:.01:20;
u=(t>0)*deg2rad(5);
phi=lsim(sysac_phi,u,t);
%[phi,t]=step(deg2rad(5)*sysac_phi,50);
phi=rad2deg(phi);
phi_ss=rad2deg(dcgain(deg2rad(5)*sysac_phi));

% figure; hold on
% plot(t,phi,'LineWidth',1.5);
% yline(phi_ss,'LineWidth',1,'LineStyle','--','Color','k');
% yline(phi_ss+1,'LineWidth',.5,'LineStyle','--','Color','k');
% yline(phi_ss-1,'LineWidth',.5,'LineStyle','--','Color','k');
% yline(5,'LineWidth',.5,'LineStyle','--','Color','r');
% xlabel('$t$ in s','Interpreter','latex','FontSize',12);
% ylabel('$\Phi$ in deg','Interpreter','latex','FontSize',12);
% xlim([0 10])
% grid on
% hold off

% PI
F=0;G=1;D=1;J=k_phiP;

% Augmented dynamics with compensator
Aa_phi2=[A zeros(size(A,1),1);-G*H F];
Ba_phi2=[B(:,3);0];
Ga_phi2=[zeros(size(B,1),1);G];
Ca_phi2=[-J*H D];
Fa_phi2=J;
Ha_phi2=[H 0];

figure;rlocus(Aa_phi2,Ba_phi2,Ca_phi2,Fa_phi2);set(findall(gca,'Type','Line'),'LineWidth',1.5,'MarkerSize',9);

k_phiI=1.1e-4;
%k_phiI=2.6e-5;

% Close the phi loop
Aac_phi2=Aa_phi2-Ba_phi2*k_phiI*Ca_phi2;
Bac_phi2=Ga_phi2-Ba_phi2*k_phiI*Fa_phi2;
sysac_phi2=ss(Aac_phi2,Bac_phi2,Ha_phi2,0);
sysac_phi2.InputName={'e_{Phi}'};
sysac_phi2.OutputName={'Phi'};