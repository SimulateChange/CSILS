%=======================================================
% Serafina Lawles, Fall 2021
% This code creates the adjacency matrix for n agents.
% The n-agent model is from the paper we titled
% "Observations 2".
%=======================================================

function A = n_agent_adj(n)

    A = zeros(3*n,3*n);
    for i = 1:n
        A(3*i-2,3*i-1) = 1;
        A(3*i-1,3*i-2) = 1;
        A(3*i,3*i-1) = 1;
        
        ii = 1+3*(i-1);
        for j = 1:n
            jj = 3*j;
            A(ii,jj) = 1;
        end
    end
end
        
        
        