% Flight Control - Group 2 - 2nd Homework
% Task 1
clearvars; close all; clc

% Import state-space model of Talon UAV (source: AlphaLink Engineering)
run('vfteStateSpace');

% controllability
Qc=ctrb(G.A,G.B(:,1:4));
if rank(Qc) == size(G.A,1)
    disp('system controllable.')
else
    disp('system uncontrollable.')
end

%disp(G.StateName)

%% Analysis eigenvalues, eigenvectors
format long

[V,D]=eig(G.A);
disp(G.StateName)

% Short period
V_sp=V(:,1)
D_sp=D(1,1)

% Phugoid
V_ph=V(:,3)
D_ph=D(3,3)

% Roll
V_r=V(:,5)
D_r=D(5,5)

% Dutch roll
V_dr=V(:,6)
D_dr=D(6,6)

% Spiral
V_spir=V(:,8)
D_spir=D(8,8)

%% Eigenstructure assignment

% Only lateral-directional motion affects the longitudinal motion
D=zeros(4,8,5);
D(:,:,1)=[1 0 0 0 0 0 0 0;0 0 1 0 0 0 0 0;
          0 0 0 0 0 1 0 0;0 0 0 0 0 0 0 1]; %D boolean matrix for short period eigenvector
D(:,:,2)=[0 1 0 0 0 0 0 0;0 0 1 0 0 0 0 0;
          0 0 0 0 0 1 0 0;0 0 0 0 0 0 0 1]; %D boolean matrix for phugoid eigenvector
D(:,:,3)=[1 0 0 0 0 0 0 0;0 1 0 0 0 0 0 0;
          0 0 0 0 0 1 0 0;0 0 0 0 0 0 1 0]; % D boolean matrix for roll subsidence eigenvector
D(:,:,4)=[0 0 1 0 0 0 0 0;0 0 0 1 0 0 0 0;
          0 0 0 0 0 1 0 0;0 0 0 0 0 0 0 1]; % D boolean matrix for Dutch roll eigenvector
D(:,:,5)=[0 1 0 0 0 0 0 0;0 0 0 1 0 0 0 0;
          0 0 0 0 0 1 0 0;0 0 0 0 0 0 0 1]; % D boolean matrix for spiral eigenvector
vd=[1 0 0 0;0 1 0 0;0 0 0 1;0 0 1 0;0 0 0 1]'; % desired lateral-directional eigenvectors (shortened)

% MIL-F-8785C (roll and spiral mode already satisfactory)
[wn_sp,~]=damp(D_sp);
zeta_sp=.707; w_sp=wn_sp/4; % short-period
eval_sp=-zeta_sp*w_sp+1j*w_sp*sqrt(1-zeta_sp^2);

[wn_ph,~]=damp(D_ph);
zeta_ph=.707; w_ph=wn_ph; % phugoid
eval_ph=-zeta_ph*w_ph+1j*w_ph*sqrt(1-zeta_ph^2);

[wn_dr,~]=damp(D_dr);
zeta_dr=.9; w_dr=wn_dr/2; % Dutch roll
eval_dr=-zeta_dr*w_dr+1j*w_dr*sqrt(1-zeta_dr^2);

desired_evals=[eval_sp eval_ph D_r eval_dr D_spir];

v=zeros(size(G.A,1),5); v_conj=zeros(size(G.A,1),5); u=zeros(4,5); u_conj=zeros(4,5);
for i=1:5
    A=[desired_evals(i)*eye(size(G.A))-G.A G.B(:,1:4);D(:,:,i) zeros(size(D,1),4)];
    b=[zeros(size(G.A,1),1);vd(:,i)];
    res=A\b;
    v(:,i)=res(1:8);
    v_conj(:,i)=conj(v(:,i));
    u(:,i)=res(9:end);
    u_conj(:,i)=conj(u(:,i));
end
K=real([u u_conj]/[v v_conj])

% verify
A_cl=G.A-G.B(:,1:4)*K;
[V_cl,D_cl]=eig(A_cl)

sys_cl=ss(A_cl,G.B(:,1:4),G.C,G.D(:,1:4));
sys_cl.InputName=G.InputName(1:4);
sys_cl.OutputName=G.OutputName;
sys_cl.StateName=G.StateName;

t=0:.01:10;
doublet=(t>1)-2*(t>2)+(t>3);
% figure('Name','elevator input'); lsim(sys_cl(:,1),doublet,t); grid on
% figure('Name','thrust input'); lsim(sys_cl(:,2),doublet,t); grid on
% figure('Name','aileron input'); lsim(sys_cl(:,3),doublet,t); grid on
figure('Name','rudder input'); lsim(sys_cl(:,4),doublet,t); grid on