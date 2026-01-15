%% state-aware mutation
function V = State_aware_mutation(P,F,xl,xu,n_dec,n_obj,nsize,K)
    dec_dist = pdist2(P(:,1:n_dec),P(:,1:n_dec));      
    NP = size(P,1);
    V = [];
    for i = 1:NP
        [~, sortedIndices] = sort(dec_dist(i, :));
        nearestIndices = sortedIndices(2:nsize+1); 
        Neighbor = P([i,nearestIndices], :); 
        D_dec = KN_dec(Neighbor,n_dec,K);  
        aver_D_dec = mean(D_dec);    
        Neighbor = NDSort(Neighbor,n_dec,n_obj);
        PS = Neighbor(Neighbor(:, end) == 1, 1:end-1);  
        isPS = ismember(P(i,:), PS, 'rows');  
        % four different states
        if isPS  
            if D_dec(1) < aver_D_dec   % good convergence and good diversity
                selected = randperm(nsize,4);  
                v = P(i,1:n_dec) + F * (Neighbor(selected(1)+1,1:n_dec)-Neighbor(selected(2)+1,1:n_dec)) +...
                    F * (Neighbor(selected(3)+1,1:n_dec)-Neighbor(selected(4)+1,1:n_dec)); 
            else    % good convergence but bad diversity
                validRange = setdiff(1:NP, i); 
                selected = datasample(validRange, 5, 'Replace', false);
                v = P(selected(1),1:n_dec) + F * (P(selected(2),1:n_dec)-P(selected(3),1:n_dec)) +...
                    F * (P(selected(4),1:n_dec)-P(selected(5),1:n_dec));
            end
        else     
            if D_dec(1) < aver_D_dec   % bad convergence but good diversity
                diff = PS(:,1:n_dec) - P(i,1:n_dec);
                dist = sum(diff .^ 2, 2);  
                [~, nearestIndex] = min(dist);  
                selected = randperm(nsize,4); 
                v = P(i,1:n_dec) + F * (PS(nearestIndex,1:n_dec)-P(i,1:n_dec)) + ...
                    F * (Neighbor(selected(1),1:n_dec)-Neighbor(selected(2),1:n_dec)) + ...
                    F * (Neighbor(selected(3),1:n_dec)-Neighbor(selected(4),1:n_dec));

            else    % bad convergence and bad diversity
                corresponding_KN_var = D_dec(Neighbor(:, end) == 1); 
                [~, min_index_in_corresponding] = min(corresponding_KN_var);
                min_index_in_HAD_var = find(Neighbor(:, end) == 1);
                minIndex = min_index_in_HAD_var(min_index_in_corresponding); 
                validRange = setdiff(1:NP, i); 
                selected = datasample(validRange, 4, 'Replace', false);
                v = P(i,1:n_dec) + F * (PS(minIndex,1:n_dec)-P(i,1:n_dec)) + ...
                    F * (P(selected(1),1:n_dec)-P(selected(2),1:n_dec)) + ...
                    F * (P(selected(3),1:n_dec)-P(selected(4),1:n_dec));
            end
        end
        for j = 1:length(v)
            if v(j) > xu(j) || v(j) < xl(j)
                v(j) = xl(j) + (xu(j) - xl(j)) * rand();
            end
        end
        V = [V;v];
    end
end