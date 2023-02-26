% Flight Control - Group 2 - 2nd Homework
% Task 2
clearvars; close all; clc
format default

% Import state-space model of Talon UAV (source: AlphaLink Engineering)
run('vfteStateSpace');

% Time constants filter
T_LP=.1;
LP=tf([T_LP],[T_LP 1]);
T_HP=2;
HP=tf([T_HP 0],[T_HP 1]);

%% 1. Design yaw damper (with washout and LP filter)
% Transfer function
TF_dr=minreal(tf(G(5,4))); % trans. fn. from rudder to yaw rate

% Determine controller gain k_zetar using root locus technique
%figure; rlocus(-TF_dr*LP*HP,.1:.01:2); set(findall(gca,'Type','Line'),'LineWidth',1.5,'MarkerSize',9);
k_zetar=.5;

TF_dr_cl=feedback(TF_dr,LP*HP*k_zetar,1);


%% 2. Design roll damper (with LP filter)
% Transfer function
TF_r=minreal(tf(G(7,3))); % trans. fn. from aileron to roll rate

% Determine controller gain k_xip using root locus technique
%figure; rlocus(-TF_r*LP); set(findall(gca,'Type','Line'),'LineWidth',1.5,'MarkerSize',9);
k_xip=.005;

TF_r_cl=feedback(TF_r,LP*k_xip,1);