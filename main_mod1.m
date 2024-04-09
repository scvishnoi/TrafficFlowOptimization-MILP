num_steps_sim = 1;
Sim_T = num_steps_sim*T;

%defining a new struct for simulation purpose called 'sim_links'

for link = 1 : num_links
    linkStr = sprintf('link_%d',link);

    sim_links.(linkStr).v = [];
    sim_links.(linkStr).Sw = [];

    if(any(link==VSL_link))
        sim_links.(linkStr).Sw = 5;
        sim_links.(linkStr).v = vsel(link,:);
        
    else
        sim_links.(linkStr).Sw = num_steps;
        sim_links.(linkStr).v = vsel(link,1);
    end
    sim_links.(linkStr).kcrit = (network_links.(linkStr).km.*network_links.(linkStr).w)./(network_links.(linkStr).w-sim_links.(linkStr).v); %critical density
    sim_links.(linkStr).q_max = sim_links.(linkStr).kcrit.*sim_links.(linkStr).v; %manetwork_links.(linkStr).ximum/critical flow
end

for link = 1: num_links
    linkStr = sprintf('link_%d',link);
    network_links.(linkStr).gridN5 = NaN*ones(length(original.(linkStr).ID_k), num_steps_sim); %downstream points from IDC
    network_links.(linkStr).gridN6 = NaN*ones(length(original.(linkStr).ID_k), num_steps_sim); %upstream points from IDC

    %Boundary Flow Conditions
    network_links.(linkStr).UBC = zeros(1, num_steps_sim);
    network_links.(linkStr).DBC = zeros(1, num_steps_sim);

    for n = 0 : num_steps_sim-1
        currsec = find(sim_links.(linkStr).Sw > n,1);
        if(isempty(currsec))
            currsec = length(sim_links.(linkStr).Sw) + 1;
        end
        network_links.(linkStr).UBC(n+1) = sim_links.(linkStr).q_max(currsec);
        network_links.(linkStr).DBC(n+1) = sim_links.(linkStr).q_max(currsec);
    end

    %Link + Simulation Object initialization
    si = sim_config(network_links.(linkStr).xi, network_links.(linkStr).chi, network_links.(linkStr).L,...
        Sim_T, T, sim_links.(linkStr).Sw, sim_links.(linkStr).v, network_links.(linkStr).w, network_links.(linkStr).km,...
        network_links.(linkStr).UBC, network_links.(linkStr).DBC);

    %initial Density conditions
    for k = 1 : length(original.(linkStr).ID_k) %original density values not the updated ones

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

        for i = 0: num_steps_sim-1
            trsec = find(sim_links.(linkStr).Sw > i,1);
            if(isempty(trsec))
                trsec = length(sim_links.(linkStr).Sw) + 1;
            end
            if(trsec > currsec)
                if(trsec ~= 1)
                    intc = network_links.(linkStr).chi - sim_links.(linkStr).v(trsec)*((i+1)*T - sim_links.(linkStr).Sw(trsec-1)*T);
                else
                    intc = network_links.(linkStr).chi - sim_links.(linkStr).v(trsec)*(i+1)*T;
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
                if(intc >= X(x_c) && intc < X(1))
                    l = find(X >= intc, 1, 'last');
                    if(l == length(X))
                        l = l-1;
                    end
                    network_links.(linkStr).gridN5(k,i+1) = N(l) + K(l)*(X(l)-intc);
                elseif(intc < X(x_c))
                    if(trsec ~= 1)
                        network_links.(linkStr).gridN5(k,i+1) = N(x_c) + si.k_c(trsec)*(sim_links.(linkStr).v(trsec)*((i+1)*T - sim_links.(linkStr).Sw(trsec-1)*T) - network_links.(linkStr).chi + X(x_c));
                    else
                        network_links.(linkStr).gridN5(k,i+1) = N(x_c) + si.k_c(trsec)*(sim_links.(linkStr).v(trsec)*(i+1)*T - network_links.(linkStr).chi + X(x_c));
                    end
                else
                    network_links.(linkStr).gridN5(k,i+1) = NaN;
                end
            end
        end

        for i = 0: num_steps_sim-1
            trsec = find(sim_links.(linkStr).Sw > i,1);
            if(isempty(trsec))
                trsec = length(sim_links.(linkStr).Sw) + 1;
            end
            if(trsec > currsec)
                if(trsec ~= 1)
                    intc = network_links.(linkStr).xi - si.w*((i+1)*T - sim_links.(linkStr).Sw(trsec-1)*T);
                else
                    intc = network_links.(linkStr).xi - si.w*(i+1)*T;
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
                if(intc <= X(x_c) && intc > X(length(X)))
                    u = find(X <= intc, 1, 'first');
                    if(u == 1)
                        u = u+1;
                    end
                    if(trsec ~= 1)
                        network_links.(linkStr).gridN6(k,i+1) = N(u) + K(u-1)*(X(u)-intc) + si.k_c(trsec)*((i+1)*T - sim_links.(linkStr).Sw(trsec-1)*T)*(sim_links.(linkStr).v(trsec) - si.w);   
                    else
                        network_links.(linkStr).gridN6(k,i+1) = N(u) + K(u-1)*(X(u)-intc) + si.k_c(trsec)*(i+1)*T*(sim_links.(linkStr).v(trsec) - si.w);
                    end
                elseif(intc > X(x_c))
                    if(trsec ~= 1)
                        network_links.(linkStr).gridN6(k,i+1) = N(x_c) + si.k_c(trsec)*(sim_links.(linkStr).v(trsec)*((i+1)*T - sim_links.(linkStr).Sw(trsec-1)*T) - network_links.(linkStr).xi + X(x_c));
                    else
                        network_links.(linkStr).gridN6(k,i+1) = N(x_c) + si.k_c(trsec)*(sim_links.(linkStr).v(trsec)*(i+1)*T - network_links.(linkStr).xi + X(x_c));
                    end
                else
                    network_links.(linkStr).gridN6(k,i+1) = NaN;
                end
            end
        end
    end

    network_links.(linkStr).gridN1 = NaN*ones(num_steps_sim, num_steps_sim); %downstream points from upstream conditions
    network_links.(linkStr).gridN2 = NaN*ones(num_steps_sim, num_steps_sim); %upstream from upstream
    network_links.(linkStr).gridN3 = NaN*ones(num_steps_sim, num_steps_sim); %downstream points from downstream conditions
    network_links.(linkStr).gridN4 = NaN*ones(num_steps_sim, num_steps_sim); %upstream from downstream
end

for n = 0 : num_steps_sim-1
    
    for i = 1 : length(InputID)           
        linkStr = sprintf('link_%d',InputID(i));
        NSupply_real(i,:) = [0 cumsum(Qin_all(i,Simulation_tracker:Simulation_tracker+num_steps-1)*T)];
        
        UBCt = NSupply_real(i,n+2)/T - sum(network_links.(linkStr).UBC(1:n)); 
        M1 = 0 + T*sum(network_links.(linkStr).UBC(1:n)); %Initial Cumulative Number of vehicles (t = 0, x = 0)
        network_links.(linkStr).UBC(n+1) = min((min(min([network_links.(linkStr).gridN2(:,n+1) network_links.(linkStr).gridN4(:,n+1)],[],'all'),min(network_links.(linkStr).gridN6(:,n+1))) - M1)/T, UBCt);  
    end
    
    for junc = 1 : num_juncs
        juncStr = sprintf('junc_%d',junc);
        if(network_junc.(juncStr).type==1)%one-one
            linkin = network_junc.(juncStr).inlabel;
            linkout = network_junc.(juncStr).outlabel;
            
            %Calculate downstream M values of input link and then get
            %possible inflow into the junction
            linkStr = sprintf('link_%d',linkin);
            M_UBConDB;
            M2 = original.(linkStr).ID_n(end) + T*sum(network_links.(linkStr).DBC(1:n));
            Qin = (min(min([network_links.(linkStr).gridN1(:,n+1) network_links.(linkStr).gridN3(:,n+1)],[],'all'),min(network_links.(linkStr).gridN5(:,n+1))) - M2)/T;
            
            %Calculate possible outflow from the junction
            linkStr = sprintf('link_%d',linkout);   
            M1 = 0 + T*sum(network_links.(linkStr).UBC(1:n)); %Initial Cumulative Number of vehicles (t = 0, x = 0)
            network_links.(linkStr).UBC(n+1) = (min(min([network_links.(linkStr).gridN2(:,n+1) network_links.(linkStr).gridN4(:,n+1)],[],'all'),min(network_links.(linkStr).gridN6(:,n+1))) - M1)/T;
            Qout = network_links.(linkStr).UBC(n+1);
            
            %Set the link flows according to final flow through the
            %junction
            linkStr = sprintf('link_%d',linkin);
            network_links.(linkStr).DBC(n+1) = min(Qin,Qout);
            linkStr = sprintf('link_%d',linkout);
            network_links.(linkStr).UBC(n+1) = min(Qin,Qout);
            
        end

        if(network_junc.(juncStr).type==2)%2-1 merge
            linkin = network_junc.(juncStr).inlabel;
            linkout = network_junc.(juncStr).outlabel;
            
            linkStr = sprintf('link_%d',linkin(1)); 
            M_UBConDB;
            M2 = original.(linkStr).ID_n(end) + T*sum(network_links.(linkStr).DBC(1:n));
            Qin1 = (min(min([network_links.(linkStr).gridN1(:,n+1) network_links.(linkStr).gridN3(:,n+1)],[],'all'),min(network_links.(linkStr).gridN5(:,n+1))) - M2)/T;
       
            linkStr = sprintf('link_%d',linkin(2)); 
            M_UBConDB;
            M2 = original.(linkStr).ID_n(end) + T*sum(network_links.(linkStr).DBC(1:n));
            Qin2 = (min(min([network_links.(linkStr).gridN1(:,n+1) network_links.(linkStr).gridN3(:,n+1)],[],'all'),min(network_links.(linkStr).gridN5(:,n+1))) - M2)/T;
            
            if(any(linkin(2)==RM_link))
               Qin2 = min(Qin2,...
                   x(network_links.(linkStr).offset+...
                   sum(network_links.(linkStr).sizes(1:iRM-1))+n+1)*sim_links.(linkStr).q_max); 
            end
            
            linkStr = sprintf('link_%d',linkout);
            M1 = 0 + T*sum(network_links.(linkStr).UBC(1:n)); %Initial Cumulative Number of vehicles (t = 0, x = 0)
            network_links.(linkStr).UBC(n+1) = (min(min([network_links.(linkStr).gridN2(:,n+1) network_links.(linkStr).gridN4(:,n+1)],[],'all'),min(network_links.(linkStr).gridN6(:,n+1))) - M1)/T;
            Qout = network_links.(linkStr).UBC(n+1);
            
            f = [-1 -1];
            Aineq = [1 0; 0 1; 1 0; 0 1];
            bineq = [];
            
            bineq(1,1) = network_junc.(juncStr).Alpha_supply(1,1)*Qout;
            bineq(2,1) = network_junc.(juncStr).Alpha_supply(1,2)*Qout;
            bineq(3,1) = Qin1;
            bineq(4,1) = Qin2;
             
            Aeq = [];
            beq = [];
            
            trsec = find(sim_links.(linkStr).Sw > n-1,1);
            if(isempty(trsec))
                trsec = length(sim_links.(linkStr).Sw) + 1;
            end
            lb = [0 0];
            linkStr = sprintf('link_%d',linkin(1));
            trsec = find(sim_links.(linkStr).Sw > n-1,1);
            if(isempty(trsec))
                trsec = length(sim_links.(linkStr).Sw) + 1;
            end
            ub = sim_links.(linkStr).q_max(trsec);
            linkStr = sprintf('link_%d',linkin(2));
            trsec = find(sim_links.(linkStr).Sw > n-1,1);
            if(isempty(trsec))
                trsec = length(sim_links.(linkStr).Sw) + 1;
            end
            ub(end+1) = sim_links.(linkStr).q_max(trsec);
            
            [x1, fval, exitflag, output] = cplexlp(f, Aineq, bineq, Aeq, beq, lb, ub);  
            
            linkStr = sprintf('link_%d',linkin(1));
            network_links.(linkStr).DBC(n+1) = x1(1);
            linkStr = sprintf('link_%d',linkin(2));
            network_links.(linkStr).DBC(n+1) = x1(2);
            linkStr = sprintf('link_%d',linkout);
            network_links.(linkStr).UBC(n+1) = x1(1)+x1(2);
        end

        if(network_junc.(juncStr).type==3)%1-2 diverge
            linkin = network_junc.(juncStr).inlabel;
            linkout = network_junc.(juncStr).outlabel;
            
            linkStr = sprintf('link_%d',linkin); 
            M_UBConDB;
            M2 = original.(linkStr).ID_n(end) + T*sum(network_links.(linkStr).DBC(1:n));
            Qin = (min(min([network_links.(linkStr).gridN1(:,n+1) network_links.(linkStr).gridN3(:,n+1)],[],'all'),min(network_links.(linkStr).gridN5(:,n+1))) - M2)/T;
            
            linkStr = sprintf('link_%d',linkout(1));
            M1 = 0 + T*sum(network_links.(linkStr).UBC(1:n)); %Initial Cumulative Number of vehicles (t = 0, x = 0)
            network_links.(linkStr).UBC(n+1) = (min(min([network_links.(linkStr).gridN2(:,n+1) network_links.(linkStr).gridN4(:,n+1)],[],'all'),min(network_links.(linkStr).gridN6(:,n+1))) - M1)/T;
            Qout1 = network_links.(linkStr).UBC(n+1);
            
            linkStr = sprintf('link_%d',linkout(2));
            M1 = 0 + T*sum(network_links.(linkStr).UBC(1:n)); %Initial Cumulative Number of vehicles (t = 0, x = 0)
            network_links.(linkStr).UBC(n+1) = (min(min([network_links.(linkStr).gridN2(:,n+1) network_links.(linkStr).gridN4(:,n+1)],[],'all'),min(network_links.(linkStr).gridN6(:,n+1))) - M1)/T;
            Qout2 = network_links.(linkStr).UBC(n+1);
            
            f = [-1 -1];
            Aineq = [1 1; 1 0; 0 1];
            bineq = [];
            
            bineq(1,1) = Qin;
            bineq(2,1) = Qout1;
            bineq(3,1) = Qout2;
            
            Aeq = [network_junc.(juncStr).Alpha(2,1) -network_junc.(juncStr).Alpha(1,1)];
            beq = 0;         
                
            trsec = find(sim_links.(linkStr).Sw > n-1,1);
            
            lb = [0 0];
            linkStr = sprintf('link_%d',linkout(1));
            trsec = find(sim_links.(linkStr).Sw > n-1,1);
            if(isempty(trsec))
                trsec = length(sim_links.(linkStr).Sw) + 1;
            end
            ub = sim_links.(linkStr).q_max(trsec);
            linkStr = sprintf('link_%d',linkout(2));
            trsec = find(sim_links.(linkStr).Sw > n-1,1);
            if(isempty(trsec))
                trsec = length(sim_links.(linkStr).Sw) + 1;
            end
            ub(end+1) = sim_links.(linkStr).q_max(trsec);
            
            [x1, fval, exitflag, output] = cplexlp(f, Aineq, bineq, Aeq, beq, lb, ub);  
            
            linkStr = sprintf('link_%d',linkin);
            network_links.(linkStr).DBC(n+1) = x1(1)+x1(2);
            linkStr = sprintf('link_%d',linkout(1));
            network_links.(linkStr).UBC(n+1) = x1(1);
            linkStr = sprintf('link_%d',linkout(2));
            network_links.(linkStr).UBC(n+1) = x1(2);
        end

    end
    
    for ithis = 1 : length(OutputID)        
        linkStr = sprintf('link_%d',OutputID(ithis));
        NDemand_real(ithis, :) = network_links.(linkStr).ID_n(end) + [0 cumsum(Qout_all(ithis,Simulation_tracker:Simulation_tracker+num_steps-1)*T)];
        
        M_UBConDB;
        DBCt = (NDemand_real(ithis, n+2) - original.(linkStr).ID_n(end))/T - sum(network_links.(linkStr).DBC(1:n));
        M2 = original.(linkStr).ID_n(end) + T*sum(network_links.(linkStr).DBC(1:n)); %Cumulative Number of vehicles (t = 0, x = L)
        network_links.(linkStr).DBC(n+1) = min((min(min([network_links.(linkStr).gridN1(:,n+1) network_links.(linkStr).gridN3(:,n+1)],[],'all'),min(network_links.(linkStr).gridN5(:,n+1))) - M2)/T,DBCt);
    end
    
    for link = 1 : num_links
        linkStr = sprintf('link_%d',link);
        M_DBConUB; % this will populate the upstream and downstream boundary grids due to the downstream boundary consition blocks
    end

end

% After running the above piece of code we get all the boundary condition
% block values. After this we can run the 'ini_den_split' code to obtain
% the split initial density for each link at the desired time step that is
% 'num_steps_sim' which will serve as the starting point for the next
% optimization.