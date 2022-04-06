%% PROJECT - 2
%% Task 3
% In this task we are performing a free stream test where height of the
% fluid in the tank is set to be 1 and not changing. In this task we look
% at the flux residuals to verify our codes.

clc; clear; close all;
format long; format compact;

% Defining initial variables and mesh components
mesh=readgri('tank0.gri');
tri_nodes=mesh.Elem;
[IE,BE]=edgehash(tri_nodes,mesh.B.nodes);
coordinates=mesh.Node;
Cfl=0.9;
T=4;       % Longer time to get more iterations for the Free stream preservation test 
p=1;       % condition variable

% Calculating the State vectors for Free Stream Tests
[~,~,k,~,~,~,~,~,~,Residual_first,Residual_last]=statecalc(tri_nodes,coordinates,IE,BE,T,Cfl,p);

% Calculating the L_2 Residual norm
% (a) Residual norm on the first iteration
A=cell2mat(Residual_first');
Resid_L2_norm=[sqrt(sum(A(1,:).^2)/(length(A(1,:)))),sqrt(sum(A(2,:).^2)/(length(A(1,:))))...
    sqrt(sum(A(3,:).^2)/(length(A(1,:))))]

% (b) Residual norm after ~2000 iterations
Number_of_iterations=k-1;
B=cell2mat(Residual_last');
Resid_L2_norm_preservation=[sqrt(sum(B(1,:).^2)/(length(B(1,:)))),sqrt(sum(B(2,:).^2)/(length(B(1,:))))...
    sqrt(sum(B(3,:).^2)/(length(B(1,:))))]