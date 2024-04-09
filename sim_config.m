%Units: Length: m, Time: s, Speed: m/s, Density: veh/m, Flow: veh/s.
classdef sim_config

    properties
        xi
        chi
        L
        Sim_T
        T
        Sw
        num_steps
        v
        w
        k_m
        k_c
        q_max
        UBC
        DBC
    end      

    methods
        
        function si = sim_config(xi, chi, L, Sim_T, T, Sw, v, w, k_m, UBC, DBC)
            %Link + Simulation Configuration
            si.xi = xi;
            si.chi = chi;
            si.L = L;
            si.Sim_T = Sim_T;
            si.T = T;
            si.num_steps = Sim_T / T;
            si.Sw = Sw;

            %Traffic Flow Parameters
            si.v = v;
            si.w = w;
            si.k_m = k_m;
            si.k_c = (k_m*w)./(w-v);
            si.q_max = si.k_c.*v;

            %Boundary Flow Conditions
            si.UBC = UBC;
            si.DBC = DBC;
        end

        function [X, K, N] = initial_projection(si, n, M)
            currsec = find(si.Sw>n,1);
            if(isempty(currsec))
                X=[]; K=[]; N=[];
                return
            end

            k = si.UBC(n+1) / si.v(currsec);
            
            xup = min(si.xi+(si.Sw(currsec)-n)*si.T*si.v(currsec),si.chi);
            xdown = min(si.xi+(si.Sw(currsec)-(n+1))*si.T*si.v(currsec),si.chi);

            if(xdown ~= si.chi)
                nup = M + si.UBC(n+1)*(si.T - ((xup - si.xi)/si.v(currsec) - (si.Sw(currsec) - (n+1))*si.T));
                ndown = M + si.UBC(n+1)*si.T + si.k_c(currsec)*(si.v(currsec)*(si.Sw(currsec)-(n+1))*si.T - xdown +si.xi); 
            else
                ndown = M + si.UBC(n+1)*si.T + si.k_c(currsec)*(si.v(currsec)*(si.Sw(currsec)-(n+1))*si.T - si.chi +si.xi);
                nup = ndown;
            end
            
            ni = M + si.UBC(n+1)*si.T + si.k_c(currsec)*si.v(currsec)*(si.Sw(currsec)-(n+1))*si.T;
            
            X = [xup xdown si.xi];
            K = [k si.k_c(currsec)];
            N = [nup ndown ni];
            
            if(X(1)==X(2))
                X = X(2:3);
                K = K(2);
                N = N(2:3);
            elseif(X(2)==X(3))
                X = X(1:2);
                K = K(1);
                N = N(1:2);
            end
        end

        function [X, K, N] = initial_projection2(si, n, M)
            currsec = find(si.Sw > n,1);
            if(isempty(currsec))
                X=[]; K=[]; N=[];
                return
            end

            k = si.k_m + si.DBC(n+1) / si.w;
            
            xup = max(si.chi+(si.Sw(currsec)-(n+1))*si.T*si.w,si.xi);
            xdown = max(si.chi+(si.Sw(currsec)-n)*si.T*si.w,si.xi);

            if(xup ~= si.xi)
                nup = M + si.DBC(n+1)*si.T + si.k_c(currsec)*(si.v(currsec)*(si.Sw(currsec)-(n+1))*si.T - xup + si.chi);
                ndown = M + si.DBC(n+1)*(si.Sw(currsec)*si.T -n*si.T +(si.chi-xdown)/si.w) + si.k_c(currsec)*(si.v(currsec)*(+xdown-si.chi)/si.w + si.chi -xdown); 
            else
                nup = M + si.DBC(n+1)*si.T + si.k_c(currsec)*(si.v(currsec)*(si.Sw(currsec)-(n+1))*si.T - si.xi + si.chi);
                ndown = nup;
            end
            
            ni = M + si.DBC(n+1)*si.T + si.k_c(currsec)*si.v(currsec)*(si.Sw(currsec)-(n+1))*si.T;
            
            X = [si.chi xup xdown];
            K = [si.k_c(currsec) k];
            N = [ni nup ndown];
            
            if(X(1)==X(2))
                X = X(2:3);
                K = K(2);
                N = N(2:3);
            elseif(X(2)==X(3))
                X = X(1:2);
                K = K(1);
                N = N(1:2);
            end
        end
        
        function [MasterX, MasterK, MasterN] = rec_dens_proj(si, currsec, X, K, N)
            if(currsec == length(si.Sw))
                MasterX = [];
                MasterK = [];
                MasterN = [];
            else
                currsec = currsec + 1;
                
                x_c = find(K < si.k_c(currsec),1,'last')+1;
                if(isempty(x_c)) 
                    x_c = 1; 
                end
                if(currsec ~= 1)
                    X1 = [X(1:x_c)+si.v(currsec)*(si.Sw(currsec)-si.Sw(currsec-1))*si.T X(x_c:length(X))+si.w*(si.Sw(currsec)-si.Sw(currsec-1))*si.T];
                    K1 = [K(1:x_c-1) si.k_c(currsec) K(x_c:length(K))];
                    N1 = [N(1:x_c) N(x_c:length(N))+si.k_c(currsec)*(si.Sw(currsec)-si.Sw(currsec-1))*si.T*(si.v(currsec)-si.w)];
                else
                    X1 = [X(1:x_c)+si.v(currsec)*(si.Sw(currsec))*si.T X(x_c:length(X))+si.w*(si.Sw(currsec))*si.T];
                    K1 = [K(1:x_c-1) si.k_c(currsec) K(x_c:length(K))];
                    N1 = [N(1:x_c) N(x_c:length(N))+si.k_c(currsec)*(si.Sw(currsec))*si.T*(si.v(currsec)-si.w)];
                end
                if(any(X1 >= si.chi))
                    K1(1:find(X1>=si.chi,1,'last')-1) = []; 
                    N1(1:find(X1>=si.chi,1,'last')) = [];
                    X1(1:find(X1>=si.chi,1,'last')) = [];
                    X1 = [si.chi X1];  
                    N1 = [N1(1)-K1(1)*(si.chi-X1(2)) N1];                    
                end
                if(any(X1 <= si.xi))
                    K1(find(X1<=si.xi,1):length(K1)) = [];
                    N1(find(X1<=si.xi,1):length(X1)) = [];
                    X1(find(X1<=si.xi,1):length(X1)) = [];
                    X1 = [X1 si.xi];
                    N1 = [N1 N1(length(N1))+K1(length(K1))*(X1(length(X1)-1)-si.xi)];                   
                end
                [tempX, tempK, tempN] = si.rec_dens_proj(currsec, X1, K1, N1);
                if(isempty(tempX))
                    MasterX = X1; MasterK = K1; MasterN = N1;                    
                elseif(size(X1,2) > size(tempX,2))
                    extra = abs(size(X1,2)-size(tempX,2));
                    MasterX = [X1; tempX ones(size(tempX,1),extra)*NaN];
                    MasterK = [K1; tempK ones(size(tempK,1),extra)*NaN];
                    MasterN = [N1; tempN ones(size(tempN,1),extra)*NaN];
                else
                    extra = abs(size(X1,2)-size(tempX,2));
                    MasterX = [X1 ones(1,extra)*NaN; tempX];
                    MasterK = [K1 ones(1,extra)*NaN; tempK];
                    MasterN = [N1 ones(1,extra)*NaN; tempN];
                end
            end
        end
        
        function[gridN] = gridcomp(si, n, MasterX, MasterK, MasterN, M, cond)
            if(cond == 0)
                currsec = 0; % initial density conditions
            else
                currsec = find(si.Sw > n, 1); % boundary conditions
            end
            dt = 1;
            dx = 1;
            gridN = zeros(si.L/dx+1, si.Sim_T/dt+1); 
            for t = 0 : dt : si.Sim_T
                t=round(t,1);
                for x = si.xi : dx : si.chi
                    x=round(x,1);
                    trsec = find(si.Sw*si.T > t,1);
                    if(isempty(trsec))
                        trsec = length(si.Sw) + 1;
                    end
                    if(trsec > currsec)
                        if(trsec ~= 1)
                            intd = x - si.v(trsec)*(t - si.Sw(trsec-1)*si.T);
                            intu = x - si.w*(t - si.Sw(trsec-1)*si.T);
                        else
                            intd = x - si.v(trsec)*t;
                            intu = x - si.w*t;
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
                        if(intd > X(1) || intu < X(length(X)))
                            gridN(x*(1/dx)+1,t*(1/dt)+1) = NaN;
                        elseif(intd >= X(x_c))
                            l = find(X >= intd, 1, 'last');
                            if(l == length(X))
                                l = l-1;
                            end
                            gridN(x*(1/dx)+1,t*(1/dt)+1) = N(l) + K(l)*(X(l)-intd);
                        elseif(intu <= X(x_c))
                            u = find(X <= intu, 1, 'first');
                            if(u == 1)
                                u = u+1;
                            end
                            if(trsec ~= 1)
                                gridN(x*(1/dx)+1,t*(1/dt)+1) = N(u) + K(u-1)*(X(u)-intu) + si.k_c(trsec)*(t - si.Sw(trsec-1)*si.T)*(si.v(trsec) - si.w);
                            else
                                gridN(x*(1/dx)+1,t*(1/dt)+1) = N(u) + K(u-1)*(X(u)-intu) + si.k_c(trsec)*(t)*(si.v(trsec) - si.w);
                            end
                        else
                            if(trsec ~= 1)
                                gridN(x*(1/dx)+1,t*(1/dt)+1) = N(x_c) + si.k_c(trsec)*(si.v(trsec)*(t - si.Sw(trsec-1)*si.T) - x + X(x_c));
                            else
                                gridN(x*(1/dx)+1,t*(1/dt)+1) = N(x_c) + si.k_c(trsec)*(si.v(trsec)*(t) - x + X(x_c));
                            end          
                        end
                    elseif(trsec == currsec)
                        if(cond == 1)
                            intd = t - x / si.v(trsec);
                            if(intd >= n*si.T && intd <= (n+1)*si.T)
                                gridN(x*(1/dx)+1,t*(1/dt)+1) = M + si.UBC(n+1)*(intd - n*si.T);
                            elseif(intd > (n+1)*si.T)
                                gridN(x*(1/dx)+1,t*(1/dt)+1) = M + si.T * si.UBC(n+1) + si.k_c(trsec)*(si.v(trsec)*(t-(n+1)*si.T) - x + si.xi);
                            else
                                gridN(x*(1/dx)+1,t*(1/dt)+1) = NaN;
                            end
                        elseif(cond == 2)
                            intu = t + (si.chi-x) / si.w;
                            if(intu >= n*si.T && intu <= (n+1)*si.T)
                                gridN(x*(1/dx)+1,t*(1/dt)+1) = M + si.DBC(n+1)*(intu - n*si.T)+si.k_c(trsec)*(-(si.chi-x) / si.w)*(si.v(trsec)-si.w);
                            elseif(intu > (n+1)*si.T)
                                gridN(x*(1/dx)+1,t*(1/dt)+1) = M + si.T * si.DBC(n+1) + si.k_c(trsec)*(si.v(trsec)*(t-(n+1)*si.T) +si.chi - x);
                            else
                                gridN(x*(1/dx)+1,t*(1/dt)+1) = NaN;
                            end
                        end
                    else
                        gridN(x*(1/dx)+1, t*(1/dt)+1) = NaN;
                    end
                end
            end
        end
           
        function[gridK] = Kcomp(si, gridN)
            dx= 1; dt = 1;
            gridK = zeros(si.L/dx, si.Sim_T/dt);
           
            for x = si.xi+dx : dx : si.chi
               x = round(x,1);
               for t = dt : dt : si.Sim_T
                   t = round(t,1);
                   gridK(x*(1/dx), t*(1/dt)) = -(gridN(x*(1/dx)+1, t*(1/dt)) - gridN(x*(1/dx), t*(1/dt)))/dx;
               end
           end
           
        end
end
end




