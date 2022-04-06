%% PROJECT 2
%% TASK 1

clc; clear; close all;
% Defining initial variables and mesh components
mesh=readgri('test.gri');
tri_nodes=mesh.Elem;
[IE,BE]=edgehash(tri_nodes,mesh.B.nodes);
coordinates=mesh.Node;
 
%% Calculating Boundary edge normals
[M,~]=size(BE);
dl_BE=zeros();
for i=1:M
    a1=coordinates(BE(i,1),:);
    a2=coordinates(BE(i,2),:);
    dl_BE(i)=norm(a1-a2);
    normal_BE{i,1}=[(a2(2)-a1(2)) (a1(1)-a2(1))]/dl_BE(i);
end

%% initializing the mesh
for i=1:(mesh.nElem)
    h=1;
    u=0; v=0;
    U{i,1}=[h,h*u,h*v];
end

%% Calculating the boundary flux uing wallflux function
for ib=1:length(BE(:,3))
    u=U{BE(ib,3),1};
    [boundary_flux{ib,1},smagBE]=wallflux(u,normal_BE{ib,1});
end
disp(cell2mat(boundary_flux'))