%M values due to DBC on UB

M2 = original.(linkStr).ID_n(end) + T*sum(network_links.(linkStr).DBC(1:n)); %Cumulative Number of vehicles (t = 0, x = L)

%Link + Simulation Object initialization
si = sim_config(network_links.(linkStr).xi, network_links.(linkStr).chi, network_links.(linkStr).L,...
    Sim_T, T, sim_links.(linkStr).Sw, sim_links.(linkStr).v, network_links.(linkStr).w, network_links.(linkStr).km,...
    network_links.(linkStr).UBC, network_links.(linkStr).DBC);


[X, K, N] = si.initial_projection2(n, M2); 
currsec = find(si.Sw > n,1);

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

for i = 0: num_steps_sim-1
    trsec = find(sim_links.(linkStr).Sw > i,1);
    if(isempty(trsec))
        trsec = length(sim_links.(linkStr).Sw) + 1;
    end
    if(trsec > currsec)
        intc = network_links.(linkStr).xi - si.w*((i+1)*T - sim_links.(linkStr).Sw(trsec-1)*T);
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
        if(intc <= X(x_c) && intc > X(length(X)))
            u = find(X <= intc, 1, 'first');
            if(u == 1)
                u = u+1;
            end
            network_links.(linkStr).gridN4(n+1,i+1) = N(u) + K(u-1)*(X(u)-intc) + si.k_c(trsec)*((i+1)*T - sim_links.(linkStr).Sw(trsec-1)*T)*(sim_links.(linkStr).v(trsec) - w);   
        elseif(intc > X(x_c))
            network_links.(linkStr).gridN4(n+1,i+1) = N(x_c) + si.k_c(trsec)*(sim_links.(linkStr).v(trsec)*((i+1)*T - sim_links.(linkStr).Sw(trsec-1)*T) - network_links.(linkStr).xi + X(x_c));
        else
            network_links.(linkStr).gridN4(n+1,i+1) = NaN;
        end
    elseif(trsec == currsec)
        intc = (i+1)*T + network_links.(linkStr).L/network_links.(linkStr).w;
        if(intc > n*T && intc <= (n+1)*T)
            network_links.(linkStr).gridN4(n+1,i+1) = M2 + network_links.(linkStr).DBC(n+1)*(intc - n*T) + si.k_c(trsec)*(sim_links.(linkStr).v(trsec)*((i+1)*T-intc) + network_links.(linkStr).L);
        elseif(intc > (n+1)*T)
            network_links.(linkStr).gridN4(n+1,i+1) = M2 + T * network_links.(linkStr).DBC(n+1) + si.k_c(trsec)*(sim_links.(linkStr).v(trsec)*((i+1)*T-(n+1)*T) + network_links.(linkStr).L);
        else
            network_links.(linkStr).gridN4(n+1,i+1) = NaN;
        end
    else
        network_links.(linkStr).gridN4(n+1,i+1) = NaN;
    end
end

for i = 0: num_steps_sim-1
    trsec = find(sim_links.(linkStr).Sw > i,1);
    if(isempty(trsec))
        trsec = length(sim_links.(linkStr).Sw) + 1;
    end
    if(trsec > currsec)
        intc = network_links.(linkStr).chi - sim_links.(linkStr).v(trsec)*((i+1)*T - sim_links.(linkStr).Sw(trsec-1)*T);
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
        if(intc >= X(x_c) && intc < X(1))
            l = find(X >= intc, 1, 'last');
            if(l == length(X))
                l = l-1;
            end
            network_links.(linkStr).gridN3(n+1,i+1) = N(l) + K(l)*(X(l)-intc);
        elseif(intc < X(x_c))
            network_links.(linkStr).gridN3(n+1,i+1) = N(x_c) + si.k_c(trsec)*(sim_links.(linkStr).v(trsec)*((i+1)*T - sim_links.(linkStr).Sw(trsec-1)*T) - network_links.(linkStr).chi + X(x_c));
        else
            network_links.(linkStr).gridN3(n+1,i+1) = NaN;
        end
    elseif(trsec == currsec)
        intc = (i+1)*T;
        if(intc > n*T && intc <= (n+1)*T)
            network_links.(linkStr).gridN3(n+1,i+1) = M2 + network_links.(linkStr).DBC(n+1)*(intc - n*T);
        elseif(intc > (n+1)*T)
            network_links.(linkStr).gridN3(n+1,i+1) = M2 + T * network_links.(linkStr).DBC(n+1) + si.k_c(trsec)*(sim_links.(linkStr).v(trsec)*((i+1)*T-(n+1)*T) - 0);
        else
            network_links.(linkStr).gridN3(n+1,i+1) = NaN;
        end
    else
        network_links.(linkStr).gridN3(n+1,i+1) = NaN;
    end
end