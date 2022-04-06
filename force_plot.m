function force_plot(k,F1,F2,F3,p)
% PURPOSE: Plots the forces on the pipes vs the total time of the
%          simulation.
%
% INPUTS: 
%    k: Total time divided into total number of iterations
%    F1: Forces in the x or y direction on Pipe 1.
%    F2: Forces in the x or y direction on Pipe 2.
%    F3: Forces in the x or y direction on Pipe 3.
%    p: condition to determine x or y direction plots.
%
% OUTPUTS:
%    Two plots one for the x-direction forces and another for the
%    y-direction forces.
if p==1
% Plotting for the x-direction forces on all the pipes.
    hold on
    plot(k,F1,'LineWidth',1.3)
    plot(k,F2,'LineWidth',1.3)
    plot(k,F3,'LineWidth',1.3)
    xlabel('Time (Seconds)')
    ylabel('Force in x-direction (N)')
elseif p==2
% Plotting for the y-direction forces on all the pipes.
    hold on
    plot(k,F1,'LineWidth',1.3)
    plot(k,F2,'LineWidth',1.3)
    plot(k,F3,'LineWidth',1.3)
    xlabel('Time (Seconds)')
    ylabel('Force in y-direction (N)')
end
end