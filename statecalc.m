function[U_final,U_half,k,F1_x,F1_y,F2_x,F2_y,F3_x,F3_y,Residual_FS,Residual_FSP]=statecalc(tri_nodes,coordinates,IE,BE,T,Cfl,p)
% PURPOSE: This function calculates the State vector of the shallow water
% equations using Forward Euler time marching method. This also calculates
% the Flux Residuals and the Forces on the three Oil Pipes
%
% INPUTS:
%    tri_nodes: The nodes matrix of the triangular cells.
%    coordinates: The x and y coordinates of thr triangular cell's nodes.
%    IE: Internal Edge matrix found using edgehash.
%    BE: Boundary Edge matrix found using edgehash. 
%    T: Final time.
%    Cfl: 
%    p: condition to specify wether to run free stream test or the normal
%       tests.
%
% OUTPUTS:
%    U_final: Conservative state vector at the final time step.
%    U_half: Conservative state vector at the halfway point of the
%            simulation.
%    k: Number of iteration/time steps performed.
%    F1_x,F1_y: Forces in the x and y direction on Pipe 1.
%    F2_x,F2_y: Forces in the x and y direction on Pipe 2.
%    F3_x,F3_y: Forces in the x and y direction on Pipe 3.
%    Residual_FS: Residual of the first iteration in the free stream test.
%    Residual_FSP: Residual of the last iteration in the free stream
%                  preservation test.

    %% Calculating the Areas and centroids
    [M,~]=size(tri_nodes);
    centroid{:}=[];
    Area=zeros();
    for i=1:M
        a1=coordinates(tri_nodes(i,1),:);
        a2=coordinates(tri_nodes(i,2),:);
        a3=coordinates(tri_nodes(i,3),:);
        centroid{i,1}=[((a1(1)+a2(1)+a3(1))/3) ((a1(2)+a2(2)+a3(2))/3)];
        len=[norm(a1-a2) norm(a2-a3) norm(a3-a1)];
        s=(len(1)+len(2)+len(3))/2;
        Area(i,1)=sqrt(s*(s-len(1))*(s-len(2))*(s-len(3)));
    end

    %% Calculating Boundary edge normals
    [M,~]=size(BE);
    dl_BE=zeros();
    normal_BE{:}=[];
    for i=1:M
        a1=coordinates(BE(i,1),:);
        a2=coordinates(BE(i,2),:);
        dl_BE(i)=norm(a1-a2);
        normal_BE{i,1}=[(a2(2)-a1(2)) (a1(1)-a2(1))]/dl_BE(i);
    end

    %% Calculating Internal edge normals
    [M,~]=size(IE);
    dl_IE=zeros();
    normal_IE{:}=[];
    for i=1:M
        a1=coordinates(IE(i,1),:);
        a2=coordinates(IE(i,2),:);
        dl_IE(i)=norm(a1-a2);
        normal_IE{i,1}=[(a2(2)-a1(2)) a1(1)-a2(1)]/dl_IE(i);
    end

    %% Intializing the mesh
    U{:}=[];
    for i=1:length(centroid)
        if p==1
        h=1;
        u=0; v=0;
        U{i,1}=[h,h*u,h*v];
        elseif p==2
        h=1+0.3*exp(-50*(centroid{i,1}(1)-1.3)^2-50*(centroid{i,1}(2)-0.9)^2);
        u=0; v=0;
        U{i,1}=[h,h*u,h*v];
        end
    end
    
    k=1;
    Residual{:}=[];
    S_sum{:}=[];
    dt=zeros();
    t=0;
    while t<=T    
        % Intialising the Residual and Speed to zeros 
        for i=1:length(centroid)
            Residual{i,1}=zeros(3,1);
            S_sum{i,1}=0;
        end

        % Calculating the Boundary Edges Flux and Speeds
        q=1; w=1; e=1;
        for ib=1:length(BE(:,3))
            u=U{BE(ib,3),1};
            [boundary_flux,smagBE]=wallflux(u,normal_BE{ib,1});
            Residual{BE(ib,3),1}=(boundary_flux).*dl_BE(ib)+Residual{BE(ib,3),1};
            S_sum{BE(ib,3),1}=smagBE*dl_BE(ib)+S_sum{BE(ib,3),1};

            % Force on pipes calculation
            if BE(ib,4)==2  % Pipe 1
                F1{q,1}=(boundary_flux)'.*dl_BE(ib);
                q=q+1;
            elseif BE(ib,4)==3  % Pipe 2
                F2{w,1}=(boundary_flux)'.*dl_BE(ib);
                w=w+1;
            elseif BE(ib,4)==4  % Pipe 3
                F3{e,1}=(boundary_flux)'.*dl_BE(ib);
                e=e+1;
            end
        end
            Q=cell2mat(F1);
            F1_x(k)=sum(Q(:,2)); F1_y(k)=sum(Q(:,3));
            W=cell2mat(F2);
            F2_x(k)=sum(W(:,2)); F2_y(k)=sum(W(:,3));
            E=cell2mat(F3);
            F3_x(k)=sum(E(:,2)); F3_y(k)=sum(E(:,3));

        % Calculating the Interior Edges Flux and Speeds 
        for ie=1:length(IE(:,3))
            uL=U{IE(ie,3),1};
            uR=U{IE(ie,4),1};
            [internal_flux,smagIE]=flux(uL,uR,normal_IE{ie,1}); 

            % Adding the residual to the left cell based on the normal
            % direction and adding the Speed 
            Residual{IE(ie,3),1}=Residual{IE(ie,3),1}+(internal_flux)'.*dl_IE(ie);
            S_sum{IE(ie,3),1}=smagIE*dl_IE(ie)+S_sum{IE(ie,3),1};

            % Subtracting the residual to the right cell based on the normal
            % direction and adding the Speed 
            Residual{IE(ie,4),1}= Residual{IE(ie,4),1}-(internal_flux)'.*dl_IE(ie);
            S_sum{IE(ie,4),1}=smagIE*dl_IE(ie)+S_sum{IE(ie,4),1};
        end

        % Calculating the Local time steps 
        for i=1:length(centroid)
            dt(i)=(2*Area(i)*Cfl)/S_sum{i,1};
        end

        % Global Time Stepping
        DT=min(dt);

        % Calculating the state at the current time step
        for i=1:length(centroid)
            U{i,1}=U{i,1}-(DT/Area(i))*Residual{i,1}';
        end

        % Increasing the time step 
        t=t+DT;
        
        % Finding the state at half time
        if (0.2480<=t)&&(t<=0.250)
            U_half=U;
        end
        
        % Storing the Residual for first iteration
        if k==1
            Resid_first_iter=Residual;
        end
        
        k=k+1;
    end
 U_final=U;
 Residual_FS=Resid_first_iter;
 Residual_FSP=Residual;
end