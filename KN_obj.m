%% k-nearest neighbor method in the objective space
function D_obj = KN_obj(P,n_dec,n_obj,K)
    dist_obj = pdist2(P(:,n_dec+1:n_dec+n_obj),P(:,n_dec+1:n_dec+n_obj));
    d_obj = zeros(size(P, 1), 1); 
    for i = 1:size(P, 1)
        dist_obj_i = dist_obj(i, :);
        dist_obj_i(i) = inf;   
        [~, idx] = sort(dist_obj_i);
        nearest_k_obj = dist_obj_i(idx(1:K));
        d_obj(i) = mean(nearest_k_obj);
    end
    d_obj_mean = mean(d_obj);
    d_obj = d_obj/d_obj_mean;
    D_obj = 1./(1 + d_obj);
end