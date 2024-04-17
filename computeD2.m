% ===============================================================
% Serafina Middleton, 2/11/2022
% Implementation of SSP algorithm for use in CSILS algorithm.
% Takes the shortest distance matrix (D) and the cell array T as inputs.
% ================================================================

function D2 = computeD2(D,T)

n = length(D);
D2 = inf(n);  % creates nxn D2 matrix of inf values

% ------------------------
% sets diag of D2 to zeros
for i=1:n  
    for j=1:n
        if i==j
            D2(i,j)=0;
        end
    end
end

% ------------------------
% find all edges a-->b
[row,col]=find(D==1);
edges = [row,col];

% ------------------------
% look for SSP >= 3

flag = 0;
for e = 1:length(edges)     % loop over all edges
    ab = edges(e,:);
    a = ab(1); b = ab(2);
    for u = 1:n             % loop over all (u,v) s.t. (u,v)~=(a,b) and u~=v
        for v = 1:n
            if u~=v
                if length([u,v,a,b]) == length(unique([u,v,a,b]))
                    D_abv = D(u,a)+D(a,b)+D(b,v);
                    if D2(u,v) > D_abv && D_abv > D(u,v) % checks first condition
                        T_ua = T{u,a};   % cell array of paths from u-->a
                        T_bv = T{b,v};   % cell array of paths from b-->v
                        for i = 1:size(T_ua,1) 
                            for j = 1:size(T_bv,1) % looping through all path combos
                                    combo = [T_ua{i}, T_bv{j}];
                                    if length(combo) == length(unique(combo))  % checks for reuse of same vertex
                                        D2(u,v) = D_abv;   % replaces D2(u,v) with SSP just found
                                        flag = 1;
                                        break % only need one path
                                    end    
                            end
                            if flag==1
                                break
                            end
                        end
                    end
                end
            end
        end
    end
end

% ================================
% look for all SSP = 2
% ================================
for i = 1:length(edges)
    uv = edges(i,:);
    u = uv(1); v=uv(2);
    row = D(u,:); row(u)=nan; row(v)=nan;
    col = D(:,v); col(u)=nan; col(v)=nan;
    loc1 = find(row==1);
    loc2 = find(col==1);
    exists = find(loc1 == loc2);
    if exists>0
        D2(u,v)=2;
    end

end
