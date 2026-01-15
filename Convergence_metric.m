%% convergence metric
function C = Convergence_metric(P,n_dec,n_obj)
    obj_values = P(:, n_dec + 1 : n_dec + n_obj); 
    min_vals = min(obj_values); 
    max_vals = max(obj_values);
    normalized_obj = (obj_values - min_vals) ./ (max_vals - min_vals);
    square_sum = sqrt(sum(normalized_obj.^2, 2));  
    min_val = min(square_sum);
    max_val = max(square_sum);
    C = (square_sum - min_val) / (max_val - min_val);
end