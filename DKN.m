%% double k-nearest neighbor method
function D = DKN(P,n_dec,n_obj,K)
    dist_dec = pdist2(P(:,1:n_dec),P(:,1:n_dec));
    dist_obj = pdist2(P(:,n_dec+1:n_dec+n_obj),P(:,n_dec+1:n_dec+n_obj));
    d_dec = zeros(size(P, 1), 1); 
    d_obj = zeros(size(P, 1), 1);
    for i = 1:size(P, 1)
        dist_dec_i = dist_dec(i, :);
        dist_dec_i(i) = inf;   
        [~, idx] = sort(dist_dec_i);    
        nearest_k_dec = dist_dec_i(idx(1:K));
        d_dec(i) = mean(nearest_k_dec);

        dist_obj_i = dist_obj(i, :);
        dist_obj_i(i) = inf;  
        [~, idx] = sort(dist_obj_i);
        nearest_k_obj = dist_obj_i(idx(1:K));
        d_obj(i) = mean(nearest_k_obj);
    end

    d_dec_mean = mean(d_dec);
    d_obj_mean = mean(d_obj);
    d_dec = d_dec/d_dec_mean;
    d_obj = d_obj/d_obj_mean;
    
    D = 1./(1+d_dec+d_obj);
end