function G = parseTrackmateEdges(this,T)

sources = T.SPOT_SOURCE_ID;
targets = T.SPOT_TARGET_ID;
costs = T.LINK_COST;

nanmask = isnan(sources);
sources = sources(~nanmask);
targets = targets(~nanmask);
costs = costs(~nanmask);

%initialize sparse graph with one node for each cell in each frame
%find a robust way to keep up with indexing
ntotal = sum(this.ncells);
ntime = this.nTime;
cidxs = zeros(ntotal,2);
idx = 1;
for ti = 1:ntime
    n = this.ncells(ti);
    cidxs(idx:idx+n-1, :) = [ti*ones(n,1), (1:n)'];
    idx = idx + n;
end
%initialize a directed graph to store tracking info (this could potentially
%be done after the frame-frame linking step instead)
G = digraph(sparse(ntotal, ntotal));
%give each node in the graph an index for the timepoint and for the index
%of the cell at that timepoint
G.Nodes.frame = cidxs(:,1);
G.Nodes.cellidx = cidxs(:,2);

G = addedge(G, sources, targets, ones(size(sources)));
G.Edges.Cost = costs;

% H = addedge(G,s,t,w);

end