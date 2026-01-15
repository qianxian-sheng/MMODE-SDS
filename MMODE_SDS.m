%% MMODE/SDS
function [ps,pf]=MMODE_SDS(func_name,xl,xu,n_obj)
    %% parameters setup
    n_dec = size(xl, 2); % dimention of decision vector
    FEs = 0;   % current fitness evaluations
    MaxFEs = 5000*n_dec;  % maximum of fitness evaluations
    NP = 100*n_dec;  % population size
    F = 0.5;   % scale factor
    CR = 0.9;   % crossover rate
    nsize = 12;  % neighborhood size
    K = 3;   % the K in the double k-nearest neighbor method
    
    %% initialization
    positions = zeros(NP,n_dec);
    for i = 1:NP
        positions(i,:) = xl + (xu - xl) .* rand(1, n_dec);
    end
    fitness = zeros(NP,n_obj);
    for i = 1:NP
        fitness(i,:) = feval(func_name,positions(i,:));
        FEs = FEs+1;
    end
    P = [positions,fitness];  
    
   %% main loop
    while FEs < MaxFEs
        V = State_aware_mutation(P,F,xl,xu,n_dec,n_obj,nsize,K);
        Offspring = Crossover(P,V,CR,n_dec,n_obj,func_name);
        FEs = FEs + NP;
        U = [P;Offspring];
        U = NDSort(U,n_dec,n_obj);
        P = Diversity_enhanced_environmental_selection(U,NP,n_dec,n_obj,K,func_name,xl,xu);
        FEs = FEs + 2;
    end
    
   %% output the results
    S = NDSort(P,n_dec,n_obj);
    non_dominated_S = S(S(:, n_dec+n_obj+1) == 1, :);  
    ps = non_dominated_S(:,1:n_dec);
    pf = non_dominated_S(:,n_dec+1:n_dec+n_obj);
end

%% crossover
function Offspring = Crossover(P,V,CR,n_dec,n_obj,func_name)
    N = size(P,1);
    Offspring = [];
    for i = 1:N
        x = P(i,1:n_dec);
        v = V(i,1:n_dec);
        u = zeros(1,length(v));
        for j = 1:length(v)
            r = rand;
            if r < CR
                u(j) = v(j);
            else
                u(j) = x(j);
            end
        end
        Offspring = [Offspring;u];
    end
    
    fitness = zeros(N,n_obj);
    for i = 1:N
        fitness(i,:) = feval(func_name,Offspring(i,:));
    end
    Offspring = [Offspring,fitness];
end