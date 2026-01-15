%% k-nearest neighbor method in the decision space
function D_dec = KN_dec(P,n_dec,K)
    dist_dec = pdist2(P(:,1:n_dec),P(:,1:n_dec));
    d_dec = zeros(size(P, 1), 1);  
    for i = 1:size(P, 1)
        dist_dec_i = dist_dec(i, :);
        dist_dec_i(i) = inf;  
        [~, idx] = sort(dist_dec_i);   
        nearest_k_dec = dist_dec_i(idx(1:K));
        d_dec(i) = mean(nearest_k_dec);
    end
    d_dec_mean = mean(d_dec);  
    d_dec = d_dec/d_dec_mean;
    D_dec = 1./(1 + d_dec);
end