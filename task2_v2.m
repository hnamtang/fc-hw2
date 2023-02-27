% Flight Control - Group 2 - 2nd Homework
% Task 2
clearvars; close all; clc

% Import state-space model of Talon UAV (source: AlphaLink Engineering)
run('vfteStateSpace');
sys=G;

%% Setup filters
% Low-pass filter
T_LP=.1;
%LP=tf([T_LP],[T_LP 1]);
A_lp=zeros(3);A_lp(1,1)=-1/T_LP;A_lp(2,2)=-1/T_LP;
B_lp=zeros(3);B_lp(1,1)=1;B_lp(2,2)=1;
C_lp=zeros(3);C_lp(1,1)=1;C_lp(2,2)=1;
D_lp=eye(3);D_lp(1,1)=0;D_lp(2,2)=0;
lp=ss(A_lp,B_lp,C_lp,D_lp);
lp.InputName=sys.OutputName([5 7 8]);
lp.OutputName(3)=sys.OutputName(8);
lp.OutputName([1 2])={'r_{lp}','p_{lp}'};

sys1=series(sys,lp,[5 7 8]); % plant + LP

% High-pass (washout) filter
T_HP=2;
%HP=tf([T_HP 0],[T_HP 1]);
A_wash=zeros(3);A_wash(1,1)=-1/T_HP;
B_wash=zeros(3);B_wash(1,1)=-1/T_HP;
C_wash=zeros(3);C_wash(1,1)=1;
D_wash=eye(3);
wash=ss(A_wash,B_wash,C_wash,D_wash);
wash.InputName=lp.OutputName;
wash.OutputName=lp.OutputName;
wash.OutputName(1)={'r_{wash}'};

sys2=series(sys1,wash); % plant + LP + HP

% % Connect low pass filter for attenuating noise
% A_lp=zeros(size(sys.C,1),size(sys.C,1));A_lp(5,5)=-1/T_LP;A_lp(7,7)=-1/T_LP;
% %A_lp=[-1/T_LP 0;0 -1/T_LP];
% B_lp=zeros(size(sys.C,1),size(sys.C,1));B_lp(5,5)=1;B_lp(7,7)=1;
% %B_lp=[0 0 0 0 1 0 0 0 0 0 0;
% %      0 0 0 0 0 0 1 0 0 0 0];
% C_lp=zeros(size(sys.C,1),size(sys.C,1));C_lp(5,5)=1;C_lp(7,7)=1;
% %C_lp=[1 0;0 1];
% D_lp=ones(size(C_lp,1),size(B_lp,2));D_lp(1,5)=0;D_lp(2,7)=0;
% lp=ss(A_lp,B_lp,C_lp,D_lp);
% lp.InputName=sys.OutputName;
% lp.OutputName=sys.OutputName;
% lp.OutputName([5 7])={'r_{lp}','p_{lp}'};
%lp.StateName={'x_{rlp}','x_{plp}'};

%sys1=series(sys,lp);

% Connect washout filter
% A_wash=-1/T_HP;
% B_wash=[-1/T_HP 0];
% C_wash=[1;0];
% D_wash=[1 0;0 1];
% wash=ss(A_wash,B_wash,C_wash,D_wash);
% wash.InputName=lp.OutputName;
% wash.OutputName={'r_{wash}','\delta p'};
% wash.StateName={'x_{wash}'};
% A_wash=zeros(size(sys.C,1),size(sys.C,1));A_wash(5,5)=-1/T_HP;
% B_wash=zeros(size(sys.C,1),size(sys.C,1));B_wash(5,5)=-1/T_HP;
% C_wash=zeros(size(sys.C,1),size(sys.C,1));C_lp(5,5)=1;
% D_wash=eye(size(C_wash,1));
% wash=ss(A_wash,B_wash,C_wash,D_wash);
% wash.InputName=lp.OutputName;
% wash.OutputName=lp.OutputName;
% wash.OutputName(5)={'r_{wash}'};


%sys2=series(sys1,wash); % plant + LP + HP

[A,B,C,D]=ssdata(sys2);


%% Design Yaw Damper
% Determine controller gain k_zetar using root locus technique
[A_r,B_r,C_r,D_r]=ssdata(minreal(sys2(1,4)));
[z_r,p_r,k_r]=ss2zp(A_r,B_r,C_r,D_r);

% k=linspace(0,5,10000);
% r=rlocus(A_r,B_r,-C_r,D_r,k);
% figure; hold on
% plot(r,'LineWidth',1.5);
% plot(real(p_r),imag(p_r),'kx','MarkerSize',9);
% plot(real(z_r),imag(z_r),'ko','MarkerSize',9);
% sgrid
% hold off

k_zetar=1.95; % from root locus

% Close yaw rate loop
Acl_r=A+B(:,4)*k_zetar*C(1,:); % positive feedback


%% Design Roll Damper (on top of yaw damper)
% Determine controller gain k_xip using root locus technique
%[A_p,B_p,C_p,D_p]=ssdata(minreal(sys2(2,3)));
[z_p,p_p,k_p]=ss2zp(Acl_r,B(:,3),C(2,:),D(2,3));

% k=linspace(0,5,10000);
% r=rlocus(Acl_r,B(:,3),-C(2,:),D(2,3),k); % k>2.96 -> unstable
% figure; hold on
% plot(r,'LineWidth',1.5);
% plot(real(p_p),imag(p_p),'kx','MarkerSize',9);
% plot(real(z_p),imag(z_p),'ko','MarkerSize',9);
% sgrid;
% hold off

k_xip=.00615; % from root locus -> slow down roll response

% Close the roll rate loop (with yaw rate loop)
Kcl_p=zeros(size(B,2),size(C,1));
Kcl_p(4,1)=k_zetar;Kcl_p(3,2)=k_xip; % positive feedback
Acl_p=A+B*Kcl_p*C;

sys_SAS=ss(Acl_p,B,C,D);
sys_SAS.InputName=sys2.InputName;
sys_SAS.OutputName=sys2.OutputName;

plotDoubletResponse(sys([5 7 8],3:4),sys_SAS(:,3:4),5);


%% Design Roll Angle Control
