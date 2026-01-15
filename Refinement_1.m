%% refinement scheme under the first condition
function [PS,remain] = Refinement_1(PS,remain,n_dec,n_obj,func_name,xl,xu,K)
    D_dec = KN_dec(PS,n_dec,K);
    [~,maxIndex_var] = max(D_dec); 
    bias = zeros(1,n_dec); 
    [neg_min_dist, pos_min_dist] = find_distances(PS, n_dec, maxIndex_var);
    for i = 1:n_dec
        r = rand();
        if r >=0 && r < neg_min_dist(i)/(neg_min_dist(i) + pos_min_dist(i))   
            bias(i) = -abs(randn() * neg_min_dist(i));
        else      
            bias(i) = abs(randn() * pos_min_dist(i));
        end
    end
    
    sparse_position = PS(maxIndex_var,1:n_dec) + bias; 
    for j = 1:length(sparse_position)
        if sparse_position(j) > xu(j) || sparse_position(j) < xl(j)
            sparse_position(j) = xl(j) + (xu(j) - xl(j)) * rand();
        end
    end
    fitness = zeros(1,n_obj);
    fitness(:) = feval(func_name,sparse_position);
    new = [sparse_position,fitness];
    PS = [PS;new];  
    C = Convergence_metric(remain,n_dec,n_obj);
    [~,maxIndex_convergence] = max(C);
    remain(maxIndex_convergence,:) = []; 
end

%% find σpos，σneg
function [neg_min_dist, pos_min_dist] = find_distances(PS, n_dec, i)
    neg_min_dist = zeros(1, n_dec);
    pos_min_dist = zeros(1, n_dec);
    for dim = 1:n_dec
        dimensions_value_var = PS(:, dim);
        target_value = dimensions_value_var(i);
        neg_values = dimensions_value_var(dimensions_value_var < target_value);
        if ~isempty(neg_values)
            neg_min_dist(dim) = min(abs(neg_values - target_value));
        else
            neg_min_dist(dim) = 0;
        end
        pos_values = dimensions_value_var(dimensions_value_var > target_value);
        if ~isempty(pos_values)
            pos_min_dist(dim) = min(abs(pos_values - target_value));
        else
            pos_min_dist(dim) = 0;
        end
    end
end