%% PROJECT - 2
%% TASK 5
% In this task we are calculating the forces exerted by the fluid on the
% walls of the Pipes and plotting the x and y component of that force vs
% Time.

clc; clear; close all;
format long; format compact;

% Defining initial variables and mesh components for Tank0 and Tank1

% Tank0 (Coarse Grid)
mesh0=readgri('tank0.gri');
tri_nodes0=mesh0.Elem;
[IE0,BE0]=edgehash(tri_nodes0,mesh0.B.nodes);
coordinates0=mesh0.Node;

% Tank1 (Fine Grid)
mesh1=readgri('tank1.gri');
tri_nodes1=mesh1.Elem;
[IE1,BE1]=edgehash(tri_nodes1,mesh1.B.nodes);
coordinates1=mesh1.Node;

Cfl=0.9;
T=0.5;
p=2;

%Calculating the forces on the pipe
% Tank0
[~,~,k0,F1_x0,F1_y0,F2_x0,F2_y0,F3_x0,F3_y0]=statecalc(tri_nodes0,coordinates0,IE0,BE0,T,Cfl,p);

% Tank1
[~,~,k1,F1_x1,F1_y1,F2_x1,F2_y1,F3_x1,F3_y1]=statecalc(tri_nodes1,coordinates1,IE1,BE1,T,Cfl,p);

%% Plotting Forces on pipe 
k0=linspace(0,0.5,k0-1);
k1=linspace(0,0.5,k1-1);
% X-component of the Forces on the Pipes
figure()
force_plot(k0,F1_x0,F2_x0,F3_x0,1)
hold on
force_plot(k1,F1_x1,F2_x1,F3_x1,1)
grid on
legend('Coarse Pipe 1','Coarse Pipe 2','Coarse Pipe 3','Fine Pipe 1','Fine Pipe 2','Fine Pipe 3','Location','eastoutside')
title('Figure 19: X-component of the Forces on the Pipes in Coarse and Fine Grid vs Time ')

% Y-component of the Forces on the Pipes
figure()
force_plot(k0,F1_y0,F2_y0,F3_y0,2)
hold on
force_plot(k1,F1_y1,F2_y1,F3_y1,2)
grid on
legend('Coarse Pipe 1','Coarse Pipe 2','Coarse Pipe 3','Fine Pipe 1','Fine Pipe 2','Fine Pipe 3','Location','eastoutside')
title('Figure 20: Y-component of the Forces on the Pipes in Coarse and Fine Grid vs Time ')



