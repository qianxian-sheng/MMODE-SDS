%% non-dominated sorting
function P = NDSort(population,n_var,n_obj) 
    [NP,~] = size(population);
    front = 1;
    F(front).f = [];
    individual = [];
    for i = 1 : NP
        individual(i).n = 0; 
        individual(i).p = []; 
        fi = population(i,n_var+1:n_var+n_obj); 
        for j = 1 : NP
            fj = population(j,n_var+1:n_var+n_obj); 
            domination = isDominated(fi,fj); 
            if domination == 1 
                individual(i).p = [individual(i).p j];
            elseif domination == -1 
                individual(i).n = individual(i).n + 1;
            end
        end
        if individual(i).n == 0 
            population(i,n_var+n_obj+1) = 1;  
            F(front).f = [F(front).f i];
        end
    end
    
    while ~isempty(F(front).f)
        Q = []; 
        for i = 1:length(F(front).f) 
            if ~isempty(individual(F(front).f(i)).p)
                for j = 1:length(individual(F(front).f(i)).p)
                    individual(individual(F(front).f(i)).p(j)).n = ... 
                        individual(individual(F(front).f(i)).p(j)).n - 1;
                    if individual(individual(F(front).f(i)).p(j)).n == 0  
                        population(individual(F(front).f(i)).p(j),n_var+n_obj+1) = front + 1;
                        Q = [Q individual(F(front).f(i)).p(j)]; 
                    end
                end
            end
        end
        front = front + 1;
        F(front).f = Q;
    end

    [~,index_of_fronts] = sort(population(:,n_obj + n_var + 1));
    for i = 1 : length(index_of_fronts)
        sorted_based_on_front(i,:) = population(index_of_fronts(i),:);
    end
P = sorted_based_on_front;
end

%% determine the dominance relationship between two solutions
function f = isDominated(fa,fb)
    is_a_dominate_b = all(fa <= fb) && any(fa < fb);
    is_b_dominate_a = all(fb <= fa) && any(fb < fa);
    if is_a_dominate_b 
        f = 1;
    elseif is_b_dominate_a
        f = -1;
    else
        f = 0;
    end
end