%% refinement scheme under the second condition
function PS = Refinement_2(PS,n_dec,n_obj,func_name,xl,xu,K)
    D_dec = KN_dec(PS,n_dec,K);
    D_obj = KN_obj(PS,n_dec,n_obj,K);

    [~,maxIndex_obj] = max(D_obj);  
    [~, sortedIndices] = sort(D_dec, 'ascend');
    minIndex_var = sortedIndices(1:2);
    maxIndex_var = sortedIndices(end); 
    
    bias1 = zeros(1,n_dec); 
    [neg_min_dist, pos_min_dist] = find_distances(PS, n_dec, minIndex_var(1));
    for i = 1:n_dec
        r = rand();
        if r >=0 && r < neg_min_dist(i)/(neg_min_dist(i) + pos_min_dist(i))   
            bias1(i) = -abs(randn() * neg_min_dist(i));
        else       
            bias1(i) = abs(randn() * pos_min_dist(i));
        end
    end
    
    bias2 = zeros(1,n_dec);  
    [neg_min_dist, pos_min_dist] = find_distances(PS, n_dec, minIndex_var(2));
    for i = 1:n_dec
        r = rand();
        if r >=0 && r < neg_min_dist(i)/(neg_min_dist(i) + pos_min_dist(i))   
            bias2(i) = -abs(randn() * neg_min_dist(i));
        else       
            bias2(i) = abs(randn() * pos_min_dist(i));
        end
    end
    
    sparse_position1 = PS(minIndex_var(1),1:n_dec) + bias1; 
    sparse_position2 = PS(minIndex_var(2),1:n_dec) + bias2;  
    for j = 1:length(sparse_position1)
        if sparse_position1(j) > xu(j) || sparse_position1(j) < xl(j)
            sparse_position1(j) = xl(j) + (xu(j) - xl(j)) * rand();
        end
        if sparse_position2(j) > xu(j) || sparse_position2(j) < xl(j)
            sparse_position2(j) = xl(j) + (xu(j) - xl(j)) * rand();
        end
    end
    fitness1 = zeros(1,n_obj);
    fitness2 = zeros(1,n_obj);
    fitness1(:) = feval(func_name,sparse_position1);
    fitness2(:) = feval(func_name,sparse_position2);
    new1 = [sparse_position1,fitness1];
    new2 = [sparse_position2,fitness2];
    
    PS([maxIndex_var,maxIndex_obj],:) = [];  
    PS = [PS;new1;new2];
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