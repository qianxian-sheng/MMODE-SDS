%% diversity-enhanced environmental selection
function P = Diversity_enhanced_environmental_selection(U,NP,n_dec,n_obj,K,func_name,xl,xu)
    P = [];
    P = [P;U(U(:, end) == 1, 1:end-1)]; 
    U(U(:, end) == 1, :) = [];
    U = U(:,1:end-1); 
    if size(P,1) < NP  % when the number of non-dominated solutions is less than NP
        C = Convergence_metric(U,n_dec,n_obj);  % convergence metric
        D = DKN(U,n_dec,n_obj,K);  %  diversity metric
        mean_C = mean(C);        
        mean_D = mean(U(:,end)); 
        I = zeros(size(U,1),1);  % hybrid convergence-diversity metric
        for i = 1:size(U,1)
            if C(i) < mean_C || D(i) < mean_D
                I(i) = min(C(i),D(i));
            else
                I(i) = max(C(i),D(i));
            end
        end
        U = [U,I];
        U = sortrows(U, size(U, 2), 'ascend'); 
        remain = U(1:NP-size(P,1),1:end-1);   
        [P,remain] = Refinement_1(P,remain,n_dec,n_obj,func_name,xl,xu,K);
        P = [P;remain];
        
    else  % when the number of non-dominated solutions exceeds NP
        D = DKN(P,n_dec,n_obj,K);
        P = [P,D];
        P = sortrows(P, size(P, 2),'ascend'); 
        P = P(1:NP,1:end-1);
        P = Refinement_2(P,n_dec,n_obj,func_name,xl,xu,K);
    end
end