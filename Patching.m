function Patching(U,tri_nodes,coordinates,j)
% PURPOSE: Plots the countour graph of the Conservative state vector's
%          values (h, u, v) on the grid provided.
%
% INPUTS:
%    U: The overall Conservative state vector 
%    tri_nodes: The nodes matrix of the triangular cells.
%    coordinates: The x and y coordinates of thr triangular cell's nodes.

    if j==1
    v = coordinates;
    f = tri_nodes;
    a=cell2mat(U);
    col=a(:,1);
    patch('Faces',f,'Vertices',v,'FaceVertexCData',col,'FaceColor','flat','EdgeColor','none');
    hold on
    h=colorbar('southoutside');
    xlabel(h,'Contours of h = fluid height (m)')
    
    elseif j==2
    v = coordinates;
    f = tri_nodes;
    a=cell2mat(U);
    col=a(:,2)./a(:,1);
    patch('Faces',f,'Vertices',v,'FaceVertexCData',col,'FaceColor','flat','EdgeColor','none');
    hold on
    h=colorbar('southoutside');
    xlabel(h,'Contours of u = fluid velocity in x-direction (m/s)')
    
    elseif j==3
    v = coordinates;
    f = tri_nodes;
    a=cell2mat(U);
    col=a(:,3)./a(:,1);
    patch('Faces',f,'Vertices',v,'FaceVertexCData',col,'FaceColor','flat','EdgeColor','none');
    hold on
    h=colorbar('southoutside');
    xlabel(h,'Contours of v = fluid velocity in y-direction (m/s)');
    end
end
