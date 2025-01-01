% M. Motamed & S. Middleton, 2/22/22
% 
% Modified version of SILS code which takes a cycle partition adjacency 
% matrix and calculates the "Shortest Independent Loop Set" as presented 
% in Oliva (2004). Resolves incomplete ILS example for an n-agent model 
% by searching second shortest paths.

function [LPS, D1, D2] = csils(C)

% ====================================================================
% Calculates distance and length matrices. 
% ====================================================================
    
    n=length(C);
    B=C+eye(n); 
    D=C;
    for i=2:n-1
        D=D+i*(logical(B^i)-logical(B^(i-1)));
    end
    L=tril(D)'+triu(D);

% ====================================================================
% For each length >=2, finds all geodetic loops of that length
% Outputs array of SILS.
% ====================================================================

    LPS = [];
    c = 1; %loop counter
    for i = 2:max(L(:))
        [row, col] = find(L==i);
        pairs = [row,col];
        for j = 1:size(pairs,1) % gives number of pairs
            if pairs(j,:) ~= 0
                u = pairs(j,1);
                v = pairs(j,2);
                loop = looptrack(u,v,D);
                for k = 1:size(pairs,1) %sets pairs in loop to zero
                    if sum(ismember(pairs(k,:),loop)) == 2
                       pairs(k,:) = 0;
                    end
                end
                if length(loop) == length(unique(loop))
                    LPS{c} = loop;
                    c = c+1;
                end
            end

        end
    end
D1 = D;
D2 = [];
% ====================================================================
% Test for complete ILS. If incomplete, search again using second
% shortest loops. (D2 gives second shortest paths.)
% ====================================================================
    
   if length(LPS) ~= sum(sum(C)) - length(C) + 1
        disp('ILS is incomplete. Searching L12.')
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Create T (formerly P) matrix
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        T = computeT(D1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Create D2 matrix
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        D2 = computeD2(D1,T);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Search for second-shortest loops in D12.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % D21=tril(D1)+triu(D2);
        % D12=tril(D2)+triu(D1);
        
        L12=tril(D2)'+triu(D1);
        L21=tril(D1)'+triu(D2);
        L2=tril(D2)' + triu(D2); %these are no longer second-shortest loops, should we remove?
                
        % max(data(~isinf(data))
        for i = 2:max(L12(~isinf(L12(:)))) %is there a better way to get max over non-inf elements?
            [row, col] = find(L12==i);
            pairs = [row,col];
            for j = 1:size(pairs,1)
                if pairs(j,:) ~= 0
                    u = pairs(j,1);
                    v = pairs(j,2);
                    loop = looptrack2(u,v,D1,D1,D2); 
                    for k = 1:size(pairs,1) %sets pairs in loop to zero
                        if sum(ismember(pairs(k,:),loop)) == 2
                            pairs(k,:) = 0;
                        end
                    end
                    if length(loop) == length(unique(loop))
                        [~,n] = size(LPS);
                        hasmatch = 0;
                        for i = 1:n
                            z = isequal(sort(loop),sort(LPS{i}));
                            hasmatch = hasmatch + z;
                        end
                        if hasmatch == 0
                            LPS{c} = loop;
                            c = c+1; 
                        end
                    end
                end

            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if length(LPS) < sum(sum(C)) - length(C) + 1
            disp('ILS is incomplete. Searching L21.')
            cont = true;
        else if length(LPS) > sum(sum(C)) - length(C) + 1
            disp('Capturing too many loops.')
            cont = false;
        else if length(LPS) == sum(sum(C)) - length(C) + 1
            disp('Full ILS Captured.')
            cont = false;
        end
        end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Search for second-shortest loops in D21.%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if cont == true
            for i = 2:max(L21(~isinf(L21(:))))
                [row, col] = find(L21==i);
                pairs = [row,col];
                for j = 1:size(pairs,1)
                    if pairs(j,:) ~= 0
                        u = pairs(j,1);
                        v = pairs(j,2);
                        if D2(v,u) ~= 0
                            loop = looptrack2(u,v,D1,D2,D1); 
                        end
                        for k = 1:size(pairs,1) %sets pairs in loop to zero
                            if sum(ismember(pairs(k,:),loop)) == 2
                                pairs(k,:) = 0;
                            end
                        end
                        if length(loop) == length(unique(loop)) && length(loop) == i
                            [~,n] = size(LPS);
                            hasmatch = 0;
                            for i = 1:n
                                z = isequal(sort(loop),sort(LPS{i}));
                                hasmatch = hasmatch + z;
                            end
                            if hasmatch == 0
                                LPS{c} = loop;
                                c = c+1;
%                                 if length(LPS) == sum(sum(C)) - length(C) + 1
%                                     break
%                                 end
                            end
                        end
                    end
                end
%                 if length(LPS) == sum(sum(C)) - length(C) + 1
%                     break
%                 end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Checking for complete ILS. %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            if length(LPS) < sum(sum(C)) - length(C) + 1
                disp('ILS is incomplete. Searching L2.')
                cont = true;
            else if length(LPS) > sum(sum(C)) - length(C) + 1
                disp('Capturing too many loops.')
                cont = false;
            else if length(LPS) == sum(sum(C)) - length(C) + 1
                disp('Full ILS Captured.')
                cont = false;
            end
            end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Search for second-shortest loops in D2.%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if cont == true
            for i = 2:max(L2(:))
                [row, col] = find(L2==i);
                pairs = [row,col];
                for j = 1:size(pairs,1)
                    if pairs(j,:) ~= 0
                        u = pairs(j,1);
                        v = pairs(j,2);
                        if D2(v,u) ~= 0 %% <<< if statement may not be necessary, CHECK
                            loop = looptrack2(u,v,D1,D2,D2); 
                        end
                        for k = 1:size(pairs,1) %sets pairs in loop to zero
                            if sum(ismember(pairs(k,:),loop)) == 2
                                pairs(k,:) = 0;
                            end
                        end
                        if length(loop) == length(unique(loop)) && length(loop) == i
                            [~,n] = size(LPS);
                            hasmatch = 0;
                            for i = 1:n
                                z = isequal(sort(loop),sort(LPS{i}));
                                hasmatch = hasmatch + z;
                            end
                            if hasmatch == 0
                                LPS{c} = loop;
                                c = c+1;
%                                 if length(LPS) == sum(sum(C)) - length(C) + 1
%                                     break
%                                 end
                            end
                        end
                    end
                end
%                 if length(LPS) == sum(sum(C)) - length(C) + 1
%                     break
%                 end
            end
            
        end
   end
  

% ====================================================================
% Function for looptrack algorithm.
% ====================================================================

function x = looptrack(u,v,D)

k = 1;
x(1) = u;
for i = 1:(D(u,v)-1)
    rh = find(D(:,v)==D(u,v)-i);
    lh = find(D(x(k),:)==1);
    k = k+1;
    for w = 1:length(rh)
        if sum(find(lh == rh(w))) ~= 0
            x(k) = rh(w);
            break
        end
    end
end

k = k+1;
x(k) = v;
for i = 1:(D(v,u)-1)
    rh = find(D(:,u)==D(v,u)-i);
    lh = find(D(x(k),:)==1);
    k = k+1;
    for w = 1:length(lh)
        if sum(find(rh == lh(w))) ~= 0
            x(k) = lh(w);
            break
        end
    end
end

% ====================================================================
% Function for looptrack2 algorithm searching SSLs.
% ====================================================================

function x = looptrack2(u,v,D1,Duv,Dvu)

k = 1;
x(1) = u;
for i = 1:(Duv(u,v)-1)
    rh1 = find(Duv(:,v)==Duv(u,v)-i);
    rh2 = find(D1(:,v)==Duv(u,v)-i);
    rh = [rh1, rh2];
    lh = find(D1(x(k),:)==1);
    %k = k+1;
    for w = 1:length(rh)
        if sum(find(lh == rh(w))) ~= 0
            k = k+1;
            x(k) = rh(w);
            break
        end
    end
end

k = k+1;
x(k) = v;
for i = 1:(Dvu(v,u)-1)
    rh1 = find(Dvu(:,u)==Dvu(v,u)-i);
    rh2 = find(D1(:,u)==Dvu(v,u)-i);
    rh = [rh1, rh2];
    lh = find(D1(x(k),:)==1);
    %k = k+1;
    for w = 1:length(lh)
        if sum(find(rh == lh(w))) ~= 0
            k = k+1;
            x(k) = lh(w);
            break
        end
    end
end

     
