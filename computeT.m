% ===============================================================
% Computes a cell array where each element T(u,v) contains all
% paths of shortest length from vertex u to vertex v.
% Uses a the included pathtrack_rec.m and takes the shortest 
% distance matrix (D) as input.
% ================================================================

function T = computeT(D) 

T = {};
[m,n] = size(D);
for u=1:m 
    for v=1:n
        T(u,v) = {pathtrack_rec(u,v,D)};
    end
end
end

 
% ====================================================================
% Mohammad Motamed, 2/17/2022
% Finds all paths of shortest path length from u to v. Based on 
% looptrack algorithm in Oliva (2004), but uses recursion and trees. 
% Takes u, v, and shortest distance matrix D.
% ====================================================================

function xcell=pathtrack_rec(u,v,D)
global xcell
xcell={}; %generate an empty cell object
xchilds(u,v,D);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function []=xchilds(x,v,D)
global xcell
i=length(x);
u=x(1);
if i < D(u,v)
    rh=find(D(:,v)==D(u,v)-i); %{w: D(w,v)==D(u,v)-i}
    lh=find(D(x(end),:)==1); %{immediate successors of x(end)}
    rlh = intersect(rh,lh);
    q = length(rlh);
    for k=1:q
        xk=[x rlh(k)];
        xchilds(xk,v,D); %recursive calls
    end
else
    xcell{end+1}=[x v];
end
end


