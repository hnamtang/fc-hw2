% Flight Control - Group 2 - 2nd Homework
% Task 1
clearvars; close all; clc

% Import state-space model of Talon UAV (source: AlphaLink Engineering)
run('vfteStateSpace');

%% Setup the state-space
% longitudinal
Along=G.A(1:4,1:4);
Blong=G.B(1:4,[1 2 5 7 9]);
Clong=G.C([1:4 9 10],1:4);
Dlong=G.D([1:4 9 10],[1 2 5 7 9]);
Glong=ss(Along,Blong,Clong,Dlong);
Glong.InputName=G.InputName([1 2 5 7 9]);
Glong.OutputName=G.OutputName([1 2 3 4 9 10]);
Glong.StateName=G.StateName(1:4);

% lateral
Alat=G.A(5:8,5:8);
Blat=G.B(5:8,[3 4 6 8]);
Clat=G.C([5:8 11],[5:8]);
Dlat=G.D([5:8 11],[3 4 6 8]);
Glat=ss(Alat,Blat,Clat,Dlat);
Glat.InputName=G.InputName([3 4 6 8]);
Glat.OutputName=G.OutputName([5:8 11]);
Glat.StateName=G.StateName(5:8);

% controllability
Qc_long=ctrb(Along,Blong(:,1:2));
if rank(Qc_long) == size(Along,1)
    disp('system controllable (longitudinal).')
else
    disp('system uncontrollable (longitudinal).')
end

Qc_lat=ctrb(Alat,Blat(:,1:2));
if rank(Qc_lat) == size(Alat,1)
    disp('system controllable (lateral-directional).')
else
    disp('system uncontrollable (lateral-directional).')
end

%% Eigenstructure assignment - longitudinal
% eigenvalues and eigenvectors
[V_long,D_long]=eig(Along);


D_long=zeros(2,4,2);
D_long(:,:,1)=[1 0 0 0;0 0 1 0]; % D boolean matrix for short-period eigenvector
D_long(:,:,2)=[0 1 0 0;0 0 0 1]; % D boolean matrix for phugoid eigenvector
vd=[1 0;0 1]'; % desired longitudinal eigenvectors (shortened)

% MIL-F-8785C
[wn_long,~]=damp(Glong);
zeta_sp=.707; w_sp=max(wn_long)/4; % short-period
eval_sp=-zeta_sp*w_sp+1j*w_sp*sqrt(1-zeta_sp^2);

zeta_ph=.707; w_ph=min(wn_long); % phugoid
eval_ph=-zeta_ph*w_ph+1j*w_ph*sqrt(1-zeta_ph^2);

desired_evals=[eval_sp eval_ph];

v=zeros(4,2); v_conj=zeros(4,2); u=zeros(2,2); u_conj=zeros(2,2);
for i=1:2
    A=[desired_evals(i)*eye(size(Along))-Along Blong(:,1:2);D_long(:,:,i) zeros(size(D_long,1),2)];
    b=[zeros(size(Along,1),1);vd(:,i)];
    res=A\b;
    v(:,i)=res(1:4);
    v_conj(:,i)=conj(v(:,i));
    u(:,i)=res(5:end);
    u_conj(:,i)=conj(u(:,i));
end
K_long=real([u u_conj]/[v v_conj])

% verify
Along_cl=Along-Blong(:,1:2)*K_long;
[V_long,D_long]=eig(Along_cl)

sys_long=ss(Along_cl,Blong(:,1:2),Clong(:,1:4),Dlong(:,1:2));
sys_long.InputName=Glong.InputName(1:2);
sys_long.OutputName=Glong.OutputName;
sys_long.StateName=Glong.StateName(1:4);

t=0:.01:10;
doublet=(t>1)-2*(t>2)+(t>3);
figure; lsim(sys_long(:,2),doublet,t); grid on

%% Eigenstructure assignment - lateral
% eigenvalues and eigenvectors
[V_lat,DD_lat]=eig(Alat);

D_lat=zeros(2,4,3);
D_lat(:,:,1)=[0 1 0 0;0 0 1 0]; % D boolean matrix for roll subsidence eigenvector
D_lat(:,:,2)=[0 1 0 0;0 0 1 0]; % D boolean matrix for Dutch roll eigenvector
D_lat(:,:,3)=[0 1 0 0;0 0 0 1]; % D boolean matrix for spiral eigenvector
vd=[0 1;1 0;0 1]'; % desired lateral-directional eigenvectors (shortened)

% MIL-F-8785C (only Dutch roll, other modes already satisfactory)
[wn_lat,~]=damp(Glat);
zeta_dr=.9; w_dr=wn_lat(2)/2;
eval_dr=-zeta_dr*w_dr+1j*w_dr*sqrt(1-zeta_dr^2);

desired_evals=[max(diag(DD_lat)) eval_dr min(diag(DD_lat))];

v=zeros(4,3); %v_conj=zeros(4,3);
u=zeros(2,3); %u_conj=zeros(2,3);
for i=1:3
    A=[desired_evals(i)*eye(size(Alat))-Alat Blat(:,1:2);D_lat(:,:,i) zeros(size(D_lat,1),2)];
    b=[zeros(size(Alat,1),1);vd(:,i)];
    res=A\b;
    v(:,i)=res(1:4);
    u(:,i)=res(5:end);
end
K_lat=real([u(:,1) u(:,2) conj(u(:,2)) u(:,3)]/[v(:,1) v(:,2) conj(v(:,2)) v(:,3)])
%K_lat(2,2)=-1;

% verify
Alat_cl=Alat-Blat(:,1:2)*K_lat;
[V_lat,D_lat]=eig(Alat_cl)

sys_lat=ss(Alat_cl,Blat(:,1:2),Clat(:,1:4),Dlat(:,1:2));
sys_lat.InputName=Glat.InputName(1:2);
sys_lat.OutputName=Glat.OutputName;
sys_lat.StateName=Glat.StateName(1:4);

t=0:.01:6;
doublet=(t>1)-2*(t>2)+(t>3);
figure; lsim(sys_lat(:,1),doublet,t); grid on