% Now we have the actual initial condition, the actual UBCs and the actual
% DBCs so we can use them to calculate the density values anywhere in the
% grid that we want
% We calculate the density values at 'num_steps_sim' that will be the
% initial density for the next run of the optimization
% num_steps_sim = num_steps;

num_steps_sim = N_c;
dx = 2;

for link = 1 : num_links
    linkStr = sprintf('link_%d',link);
    L = network_links.(linkStr).L;
    network_links.(linkStr).Ini_DenM = NaN*zeros(length(0:dx:L),1);
    network_links.(linkStr).Ini_Den = NaN*zeros(length(1:dx:L),1);
    
    %to obtain values at time step num_steps_sim, we need to consider UBC
    %at num_steps_sim-1 and just one before a switch if present after a
    %switch
    Sim_T = num_steps_sim*T;
    %the following 3 'for' loops store the M values at the required time
    %step throughout the link space in the Ini_DenM array
    for cond = 0 : num_steps_sim-1 %for upstream conditions
       
        if(cond ~= num_steps_sim-1 && ~any(cond==(sim_links.(linkStr).Sw-1)))
            continue;
        end
        
        M1 = 0 + T*sum(network_links.(linkStr).UBC(1:cond)); %Initial Cumulative Number of vehicles (t = 0, x = 0)

        %Link + Simulation Object initialization
        si = sim_config(network_links.(linkStr).xi, network_links.(linkStr).chi, network_links.(linkStr).L,...
            Sim_T, T, sim_links.(linkStr).Sw, sim_links.(linkStr).v, network_links.(linkStr).w, network_links.(linkStr).km,...
            network_links.(linkStr).UBC, network_links.(linkStr).DBC);

        [X, K, N] = si.initial_projection(cond, M1); 
        currsec = find(si.Sw > cond,1);

        MasterX = []; MasterK = []; MasterN = [];
        if(~isempty(currsec))
            [MasterX, MasterK, MasterN] = si.rec_dens_proj(currsec, X, K, N);
        end
        if(isempty(currsec))
            currsec = length(sim_links.(linkStr).Sw) + 1;
        end
        if(isempty(MasterX))
            MasterX = X; MasterK = K; MasterN = N;
        elseif(size(X,2) >= size(MasterX,2))
            extra = abs(size(X,2)-size(MasterX,2));
            MasterX = [X; MasterX ones(size(MasterX,1),extra)*NaN];
            MasterK = [K; MasterK ones(size(MasterK,1),extra)*NaN];
            MasterN = [N; MasterN ones(size(MasterN,1),extra)*NaN];
        elseif(size(X,2) < size(MasterX,2))
            extra = abs(size(X,2)-size(MasterX,2));
            MasterX = [X ones(1,extra)*NaN; MasterX];
            MasterK = [K ones(1,extra)*NaN; MasterK];
            MasterN = [N ones(1,extra)*NaN; MasterN];
        end

        i = num_steps_sim-1;
        trsec = find(sim_links.(linkStr).Sw > i,1);
        if(isempty(trsec))
            trsec = length(sim_links.(linkStr).Sw) + 1;
        end
        if(trsec > currsec)
            count = 1;
            for linkPos = 1 : dx : L
                intc1 = linkPos - sim_links.(linkStr).v(trsec)*((i+1)*T - sim_links.(linkStr).Sw(trsec-1)*T);
                intc2 = linkPos - si.w*((i+1)*T - sim_links.(linkStr).Sw(trsec-1)*T);
                K = MasterK(trsec-currsec,:); 
                X = MasterX(trsec-currsec,:);
                N = MasterN(trsec-currsec,:);
                if(any(isnan(K)))
                    K = K(1:find(isnan(K))-1);
                    X = X(1:find(isnan(X))-1);
                    N = N(1:find(isnan(N))-1);
                end
                x_c = find(K < si.k_c(trsec),1,'last')+1;
                if(isempty(x_c)) 
                    x_c = 1; 
                end
                if(intc1 >= X(x_c) && intc1 <= X(1))
                    l = find(X >= intc1, 1, 'last');
                    if(l == length(X))
                        l = l-1;
                    end
                    temp = N(l) + K(l)*(X(l)-intc1);
                    network_links.(linkStr).Ini_DenM(count,1) = min(temp,network_links.(linkStr).Ini_DenM(count,1));
                elseif(intc1 < X(x_c) && intc2 >= X(x_c))
                    temp = N(x_c) + si.k_c(trsec)*(sim_links.(linkStr).v(trsec)*((i+1)*T - sim_links.(linkStr).Sw(trsec-1)*T) - linkPos + X(x_c));
                    network_links.(linkStr).Ini_DenM(count,1) = min(temp,network_links.(linkStr).Ini_DenM(count,1));
                elseif(intc2 < X(x_c) && intc2 >= X(end))
                    u = find(X <= intc2, 1, 'first');
                    if(u == 1)
                        u = u+1;
                    end
                    temp = N(u) + K(u-1)*(X(u)-intc2) + si.k_c(trsec)*((i+1)*T - sim_links.(linkStr).Sw(trsec-1)*T)*(sim_links.(linkStr).v(trsec) - w);
                    network_links.(linkStr).Ini_DenM(count,1) = min(temp,network_links.(linkStr).Ini_DenM(count,1));
                end
                count = count + 1;
            end
        elseif(trsec == currsec)
            count = 1;
            for linkPos = 1 : dx : L
                intc = (i+1)*T - linkPos / sim_links.(linkStr).v(trsec);
                if(intc > cond*T && intc <= (cond+1)*T)
                    network_links.(linkStr).Ini_DenM(count,1) = M1 + network_links.(linkStr).UBC(cond+1)*(intc - cond*T);
                elseif(intc > (cond+1)*T)
                    network_links.(linkStr).Ini_DenM(count,1) = M1 + T * network_links.(linkStr).UBC(cond+1) + si.k_c(trsec)*(sim_links.(linkStr).v(trsec)*((i+1)*T-(cond+1)*T) - linkPos + network_links.(linkStr).xi);
                end
                count = count + 1;
            end
        end
    end
    
    for cond = 0 : num_steps_sim-1 %for downstream conditions
        
        M2 = original.(linkStr).ID_n(end) + T*sum(network_links.(linkStr).DBC(1:cond)); %Cumulative Number of vehicles (t = 0, x = L)

        %Link + Simulation Object initialization
        si = sim_config(network_links.(linkStr).xi, network_links.(linkStr).chi, network_links.(linkStr).L,...
            Sim_T, T, sim_links.(linkStr).Sw, sim_links.(linkStr).v, network_links.(linkStr).w, network_links.(linkStr).km,...
            network_links.(linkStr).UBC, network_links.(linkStr).DBC);


        [X, K, N] = si.initial_projection2(cond, M2);
        currsec = find(si.Sw > cond,1);

        MasterX = []; MasterK = []; MasterN = [];
        if(~isempty(currsec))
            [MasterX, MasterK, MasterN] = si.rec_dens_proj(currsec, X, K, N);
        end
        if(isempty(currsec))
            currsec = length(sim_links.(linkStr).Sw) + 1;
        end
        if(isempty(MasterX))
            MasterX = X; MasterK = K; MasterN = N;
        elseif(size(X,2) >= size(MasterX,2))
            extra = abs(size(X,2)-size(MasterX,2));
            MasterX = [X; MasterX ones(size(MasterX,1),extra)*NaN];
            MasterK = [K; MasterK ones(size(MasterK,1),extra)*NaN];
            MasterN = [N; MasterN ones(size(MasterN,1),extra)*NaN];
        elseif(size(X,2) < size(MasterX,2))
            extra = abs(size(X,2)-size(MasterX,2));
            MasterX = [X ones(1,extra)*NaN; MasterX];
            MasterK = [K ones(1,extra)*NaN; MasterK];
            MasterN = [N ones(1,extra)*NaN; MasterN];
        end 

        i = num_steps_sim-1;
        trsec = find(sim_links.(linkStr).Sw > i,1);
        if(isempty(trsec))
            trsec = length(sim_links.(linkStr).Sw) + 1;
        end
        if(trsec > currsec)
            count = 1;
            for linkPos = 1 : dx : L
                intc1 = linkPos - si.w*((i+1)*T - sim_links.(linkStr).Sw(trsec-1)*T);
                intc2 = linkPos - sim_links.(linkStr).v(trsec)*((i+1)*T - sim_links.(linkStr).Sw(trsec-1)*T);
                K = MasterK(trsec-currsec,:); 
                X = MasterX(trsec-currsec,:);
                N = MasterN(trsec-currsec,:);
                if(any(isnan(K)))
                    K = K(1:find(isnan(K))-1);
                    X = X(1:find(isnan(X))-1);
                    N = N(1:find(isnan(N))-1);
                end
                x_c = find(K < si.k_c(trsec),1,'last')+1;
                if(isempty(x_c)) 
                    x_c = 1; 
                end
                if(intc1 <= X(x_c) && intc1 >= X(length(X)))
                    u = find(X <= intc1, 1, 'first');
                    if(u == 1)
                        u = u+1;
                    end
                    temp = N(u) + K(u-1)*(X(u)-intc1) + si.k_c(trsec)*((i+1)*T - sim_links.(linkStr).Sw(trsec-1)*T)*(sim_links.(linkStr).v(trsec) - si.w);   
                    network_links.(linkStr).Ini_DenM(count,1) = min(network_links.(linkStr).Ini_DenM(count,1),temp);
                elseif(intc1 > X(x_c) && intc2 <= X(x_c))
                    temp = N(x_c) + si.k_c(trsec)*(sim_links.(linkStr).v(trsec)*((i+1)*T - sim_links.(linkStr).Sw(trsec-1)*T) - linkPos + X(x_c));
                    network_links.(linkStr).Ini_DenM(count,1) = min(network_links.(linkStr).Ini_DenM(count,1),temp);
                elseif(intc2 >= X(x_c) && intc2 <= X(1))
                    l = find(X >= intc2, 1, 'last');
                    if(l == length(X))
                        l = l-1;
                    end
                    temp = N(l) + K(l)*(X(l)-intc2);
                    network_links.(linkStr).Ini_DenM(count,1) = min(network_links.(linkStr).Ini_DenM(count,1),temp);
                end
                count = count + 1;
            end
        elseif(trsec == currsec)
            count = 1;
            for linkPos = 1 : dx : L
                intc = (i+1)*T + (L-linkPos)/network_links.(linkStr).w;
                if(intc > cond*T && intc <= (cond+1)*T)
                    temp = M2 + network_links.(linkStr).DBC(cond+1)*(intc - cond*T) + si.k_c(trsec)*(sim_links.(linkStr).v(trsec)*((i+1)*T-intc) + (L-linkPos));
                    network_links.(linkStr).Ini_DenM(count,1) = min(network_links.(linkStr).Ini_DenM(count,1),temp);
                elseif(intc > (cond+1)*T)
                    temp = M2 + T * network_links.(linkStr).DBC(cond+1) + si.k_c(trsec)*(sim_links.(linkStr).v(trsec)*((i+1)*T-(cond+1)*T) + (L-linkPos));
                    network_links.(linkStr).Ini_DenM(count,1) = min(network_links.(linkStr).Ini_DenM(count,1),temp);
                end
                count = count + 1;
            end
        end

    end
    
    for k = 1 : length(original.(linkStr).ID_k) %for initial density conditions

        X = original.(linkStr).ID_x(k+1:-1:k);
        K = original.(linkStr).ID_k(k);
        N = original.(linkStr).ID_n(k+1:-1:k);

        currsec = 0;

        [MasterX, MasterK, MasterN] = si.rec_dens_proj(currsec, X, K, N);

        if(isempty(MasterX))
            MasterX = X; MasterK = K; MasterN = N;
        elseif(size(X,2) >= size(MasterX,2))
            extra = abs(size(X,2)-size(MasterX,2));
            MasterX = [X; MasterX ones(size(MasterX,1),extra)*NaN];
            MasterK = [K; MasterK ones(size(MasterK,1),extra)*NaN];
            MasterN = [N; MasterN ones(size(MasterN,1),extra)*NaN];
        elseif(size(X,2) < size(MasterX,2))
            extra = abs(size(X,2)-size(MasterX,2));
            MasterX = [X ones(1,extra)*NaN; MasterX];
            MasterK = [K ones(1,extra)*NaN; MasterK];
            MasterN = [N ones(1,extra)*NaN; MasterN];
        end

        i = num_steps_sim-1;
        trsec = find(sim_links.(linkStr).Sw > i,1);
        if(isempty(trsec))
            trsec = length(sim_links.(linkStr).Sw) + 1;
        end
        if(trsec > currsec)
            count = 1;
            for linkPos = 1 : dx : L
                if(trsec ~= 1)
                    intc1 = linkPos - sim_links.(linkStr).v(trsec)*((i+1)*T - sim_links.(linkStr).Sw(trsec-1)*T);
                    intc2 = linkPos - si.w*((i+1)*T - sim_links.(linkStr).Sw(trsec-1)*T);
                else
                    intc1 = linkPos - sim_links.(linkStr).v(trsec)*(i+1)*T;
                    intc2 = linkPos - si.w*(i+1)*T;
                end
                K = MasterK(trsec-currsec,:); 
                X = MasterX(trsec-currsec,:);
                N = MasterN(trsec-currsec,:);
                if(any(isnan(K)))
                    K = K(1:find(isnan(K))-1);
                    X = X(1:find(isnan(X))-1);
                    N = N(1:find(isnan(N))-1);
                end
                x_c = find(K < si.k_c(trsec),1,'last')+1;
                if(isempty(x_c)) 
                    x_c = 1; 
                end
                if(intc1 >= X(x_c) && intc1 < X(1))
                    l = find(X >= intc, 1, 'last');
                    if(l == length(X))
                        l = l-1;
                    end
                    temp = N(l) + K(l)*(X(l)-intc1);
                    network_links.(linkStr).Ini_DenM(count,1) = min(network_links.(linkStr).Ini_DenM(count,1),temp);
                elseif(intc1 < X(x_c) && intc2 >= X(x_c))
                    if(trsec ~= 1)
                        temp = N(x_c) + si.k_c(trsec)*(sim_links.(linkStr).v(trsec)*((i+1)*T - sim_links.(linkStr).Sw(trsec-1)*T) - linkPos + X(x_c));
                        network_links.(linkStr).Ini_DenM(count,1) = min(network_links.(linkStr).Ini_DenM(count,1),temp);
                    else
                        temp = N(x_c) + si.k_c(trsec)*(sim_links.(linkStr).v(trsec)*(i+1)*T - linkPos + X(x_c));
                        network_links.(linkStr).Ini_DenM(count,1) = min(network_links.(linkStr).Ini_DenM(count,1),temp);
                    end
                elseif(intc2 <= X(x_c) && intc2 >= X(length(X)))
                    u = find(X <= intc2, 1, 'first');
                    if(u == 1)
                        u = u+1;
                    end
                    if(trsec ~= 1)
                        temp = N(u) + K(u-1)*(X(u)-intc2) + si.k_c(trsec)*((i+1)*T - sim_links.(linkStr).Sw(trsec-1)*T)*(sim_links.(linkStr).v(trsec) - si.w);  
                        network_links.(linkStr).Ini_DenM(count,1) = min(network_links.(linkStr).Ini_DenM(count,1),temp);
                    else
                        temp = N(u) + K(u-1)*(X(u)-intc2) + si.k_c(trsec)*(i+1)*T*(sim_links.(linkStr).v(trsec) - si.w);
                        network_links.(linkStr).Ini_DenM(count,1) = min(network_links.(linkStr).Ini_DenM(count,1),temp);
                    end
                end
                count = count + 1;
            end
        end
    end
    
    %Now we can calculate and store the values of the initial density in
    %the Ini_Den array
    for d = 1 : length(network_links.(linkStr).Ini_DenM)-1
        network_links.(linkStr).Ini_Den(d) = -(network_links.(linkStr).Ini_DenM(d+1)-network_links.(linkStr).Ini_DenM(d))/dx;
    end
    
    network_links.(linkStr).Ini_Den = network_links.(linkStr).Ini_Den(1:find(isnan(network_links.(linkStr).Ini_Den),1)-1);
    %Now we want to split the initial density in atmost three pieces
    %depending upon the change in density across the link. In the following
    %code we try to obtain three (or less) spliting positions for the array 
    %such that the density within the split up pieces has least variance.
    
    changepts = findchangepts(round(network_links.(linkStr).Ini_Den,5),'MaxNumChanges',2);
    %the MATLAB function findchangepts() provides the points where the densities
    %change in the Ini_Den array
    %rounding to 5 sd is done to avoid machine caused changes
    
    if(isempty(changepts))%if the function does not give any change points so the density is uniform
        sim_links.(linkStr).newID_x = [0 L];
        sim_links.(linkStr).newID_k = network_links.(linkStr).Ini_Den(1);
        continue;
    elseif(length(changepts)==1 && changepts == length(network_links.(linkStr).Ini_Den))%if the function returns the end of the array
                                                                                        %which also means that there is no significant change in the
                                                                                        %density throughout the link
        sim_links.(linkStr).newID_x = [0 L];
        sim_links.(linkStr).newID_k = network_links.(linkStr).Ini_Den(1);
        continue;
    else %if there are actual change points one or two in number
        if(length(changepts)==2 && changepts(1) == changepts(2)-1)
            changepts = changepts(1);
        end
        sim_links.(linkStr).newID_x = [0 changepts'*dx L];
        if(length(changepts)==1)
            sim_links.(linkStr).newID_k = [network_links.(linkStr).Ini_Den(max(1,changepts-2)) network_links.(linkStr).Ini_Den(min(changepts+2,length(network_links.(linkStr).Ini_Den)))];
        else
%             sim_links.(linkStr).newID_k = [network_links.(linkStr).Ini_Den(max(1,changepts(1)-2)) network_links.(linkStr).Ini_Den(changepts(1)+2) network_links.(linkStr).Ini_Den(min(changepts(2)+2,length(network_links.(linkStr).Ini_Den)))];
%             sim_links.(linkStr).newID_k = [trimmean(network_links.(linkStr).Ini_Den(1:changepts(1)),10) trimmean(network_links.(linkStr).Ini_Den(changepts(1):changepts(2)),10) trimmean(network_links.(linkStr).Ini_Den(changepts(2):end),10)];
            sim_links.(linkStr).newID_k = [trimmean(network_links.(linkStr).Ini_Den(1:changepts(1)),1) trimmean(network_links.(linkStr).Ini_Den(changepts(1):changepts(2)),1) trimmean(network_links.(linkStr).Ini_Den(changepts(2):end),1)];
        end
    end
end

%% Objective Calculation 
% not considering penalty terms in this calculation
fval_all = 0;
for link = 1 : num_links
    linkStr = sprintf('link_%d',link);
    M1 = 0 + T*sum(network_links.(linkStr).UBC(1:cond)); %Initial Cumulative Number of vehicles (t = 0, x = 0)
    M2 = original.(linkStr).ID_n(end) + T*sum(network_links.(linkStr).DBC(1:cond));
    if(any(link==VSL_link))
        fval_link = network_links.(linkStr).L*max((M1-M2)/network_links.(linkStr).L,sim_links.(linkStr).kcrit(1));
    else
        fval_link = network_links.(linkStr).L*max((M1-M2)/network_links.(linkStr).L,sim_links.(linkStr).kcrit(1));
    end
    fval_all = fval_all + fval_link;
end
oTAC(Simulation_tracker,1) = fval_all;

fval_all = 0;
for input = 1 : length(InputID)
    link = InputID(input);
    linkStr = sprintf('link_%d',link);
    
    carry_demand = NSupply_real(input, num_steps_sim+1)/T - sum(network_links.(linkStr).UBC(1:num_steps_sim));
    Qin_all(input,Simulation_tracker+1) = Qin_all(input,Simulation_tracker+1) + carry_demand;
    for n = 1 : num_steps_sim
        fval_all = fval_all + we(input,n)*((NSupply_real(input,n+1)-NSupply_real(input,n))-T*network_links.(linkStr).UBC(n));
%         fval_all = fval_all + T*((NSupply_real(input,n+1)-NSupply_real(input,n))-T*network_links.(linkStr).UBC(n));    
    end
end
oInput(Simulation_tracker,1) = fval_all;


%Total Travel Time Objective
fval_all = 0;
for link = 1 : num_links
    linkStr = sprintf('link_%d',link);
    fval_link = 0;
    for n = 1 : num_steps_sim    
        M1 = 0 + T*sum(network_links.(linkStr).UBC(1:n)); %Initial Cumulative Number of vehicles (t = 0, x = 0)
        M2 = original.(linkStr).ID_n(end) + T*sum(network_links.(linkStr).DBC(1:n));

        fval_link = fval_link + T*(M1-M2);

    end
    fval_all = fval_all + fval_link;
end
oTTT(Simulation_tracker,1) = fval_all;
