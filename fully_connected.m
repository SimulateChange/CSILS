% Creates a fully connected adjacency matrix of size nxn
% with no self loops. i.e. a matrix of ones with zeros on the diag

function F = fully_connected(n)
    F = ones(n);
    F(logical(eye(n))) = 0;
end

