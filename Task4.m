%% PROJECT - 2
%% TASK 4
% In this task we are calculating Conservative state vectors at a fixed
% time using an initial condition provided to us for the h vector. Contours
% are plotted for all three variables (h, u, v) in the state vector to show 
% the state at half time and at the final time.

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
T=0.5;          % Final time
p=2;            % Condition variable

%% Calculating the Conservative state vectors
% Tank0 (Coarse Grid)
[U_final0,U_half0]=statecalc(tri_nodes0,coordinates0,IE0,BE0,T,Cfl,p);

% Tank1 (Fine Grid)
[U_final1,U_half1]=statecalc(tri_nodes1,coordinates1,IE1,BE1,T,Cfl,p);

%% Plotting the States
% Tank0 (Coarse Grid)
% Plotting for final state at T=0.5 seconds
figure(1)
Patching(U_final0,tri_nodes0,coordinates0,1)
title('Figure 7: Coarse Grid Contours of h at T=0.5s')
figure(2)
Patching(U_final0,tri_nodes0,coordinates0,2)
title('Figure 8: Coarse Grid Contours of u (velocity in x-direction) at T=0.5s')
figure(3)
Patching(U_final0,tri_nodes0,coordinates0,3)
title('Figure 9: Coarse Grid Contours of v (velocity in y-direction) at T=0.5s')

% Plotting State for half time T~0.25 seconds
figure(4)
Patching(U_half0,tri_nodes0,coordinates0,1)
title('Figure 10: Coarse Grid Contours of h at T=0.25s')
figure(5)
Patching(U_half0,tri_nodes0,coordinates0,2)
title('Figure 11: Coarse Grid Contours of u (velocity in x-direction) at T=0.25s')
figure(6)
Patching(U_half0,tri_nodes0,coordinates0,3)
title('Figure 12: Coarse Grid Contours of v (velocity in y-direction) at T=0.25s')

% Tank1 (Fine Grid)
% Plotting for final state at T=0.5 seconds
figure(7)
Patching(U_final1,tri_nodes1,coordinates1,1)
title('Figure 13: Fine Grid Contours of h at T=0.5s')
figure(8)
Patching(U_final1,tri_nodes1,coordinates1,2)
title('Figure 14: Fine Grid Contours of u (velocity in x-direction) at T=0.5s')
figure(9)
Patching(U_final1,tri_nodes1,coordinates1,3)
title('Figure 15: Fine Grid Contours of v (velocity in y-direction) at T=0.5s')

% Plotting State for half time T~0.25 seconds
figure(10)
Patching(U_half1,tri_nodes1,coordinates1,1)
title('Figure 16: Fine Grid Contours of h at T=0.25s')
figure(11)
Patching(U_half1,tri_nodes1,coordinates1,2)
title('Figure 17: Fine Grid Contours of u (velocity in x-direction) at T=0.25s')
figure(12)
Patching(U_half1,tri_nodes1,coordinates1,3)
title('Figure 18: Fine Grid Contours of v (velocity in y-direction) at T=0.25s')



