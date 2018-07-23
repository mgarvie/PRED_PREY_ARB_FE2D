function cpp = subsetconnectivity (t, p, nds)

% Discussion
%
%    Determines the connectivity of a subset of boundary nodes in a
%    triangulation
%
% Modified: 
%       
%    June 20, 2009
%
% Author:
%
%    Marcus Garvie
%
% Parameters
%
%    Input, real P(NP,2), the coordinates of a set of nodes.
%
%    Input, integer T(NT,3), a list of the nodes which make up each 
%    triangle of a triangulation of the nodes in P.
%
%    Input, integer NDS(NN,1), a list of nodes comprising a subset of the
%    boundary nodes
%
%    Output, CPP(NE,2), a list of edges, each edge defined by 2 nodes from 
%    the list NDS

cpp = [];

% Work out connectivity for the whole boundary
edges = boundedges (p,t);
% Work out some dimensions
[no_edges,junk] = size(edges);  % Number of edges in whole boundary
NN = length(nds); % Number of nodes in the subset of boundary nodes
% Find the edges in the subset of the boundary NDS 
count = 0;
for i = 1:no_edges
    node1 = edges(i,1);
    node2 = edges(i,2);
    bin1 = ismember(node1,nds);
    bin2 = ismember(node2,nds);
    if (bin1==1)&(bin2==1)
        count = count + 1;
        cpp(count,:) = edges(i,:);
    end
end