function [IE,BE] = edgehash(E2N, Bnodes)
% This function identifies interior and boundary edges, and their
% connectivities, in a triangular mesh given an element-to-node array.
%
% INPUT : E2N = [nelem x 3] array mapping triangles to nodes
%         Bnodes = cell array of boundary edge nodes, as identified
%                  by readgri.m and stored in Mesh.B.nodes (optional)
% OUTPUT: IE = [niedge x 4] array giving (n1, n2, elem1, elem2)
%              information for each interior edge
%         BE = [nbedge x 4] array giving (n1, n2, elem, bgroup)
%              information for each boundary edge

nelem = size(E2N,1);            % number of elements
nnode = max(max(E2N));          % number of nodes
H = sparse(nnode,nnode);        % Create a hash list to identify edges
IE = zeros(ceil(nelem*3/2), 4); % (over) allocate interior edge array
niedge = 0;                     % number of interior edges (running total)

% Loop over elements and identify all edges
for elem = 1:nelem
  nv = E2N(elem,1:3);
  for edge = 1:3
    n1 = nv(edge);
    n2 = nv(mod(edge,3)+1);
    if (H(n2,n1) == 0)
      H(n1,n2) = elem;
    else
      elemR = H(n2,n1);
      niedge = niedge+1;
      IE(niedge,:) = [n1,n2, elem, elemR];
      H(n2,n1) = 0;
    end
  end
end

IE = IE(1:niedge,:);  % clip IE

% find boundary edges
if isempty(Bnodes)
  [I,J] = find(triu(H)>0);  
  BE = [I, J, zeros(size(I)), zeros(size(I))];
  for b = 1:size(I,1), BE(b, 3) = H(I(b),J(b)); end
else
  nb0 = 0; nb = 0;
  for g = 1:length(Bnodes), nb0 = nb0 + size(Bnodes{g}, 1); end
  BE = zeros(nb0, 4);
  for g = 1:length(Bnodes)
    Bi = Bnodes{g};
    for b = 1:size(Bi,1)
      n1 = Bi(b,1); n2 = Bi(b,2);
      if (H(n1,n2) == 0), nt = n1; n1 = n2; n2 = nt; end
      nb = nb + 1;
      BE(nb,:) = [n1, n2, H(n1,n2), g];
    end
  end
end


