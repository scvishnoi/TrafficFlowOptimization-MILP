classdef setConstraints_noVSL
    properties
        L;
        chi;
        xi;

        v;
        w;
        km;
        kc;
        
        sizes;
        size_row;
        
        nb_upstream;
        nb_downstream;
        
        junc_indicator
        network_junc;
        
        ID_x;
        ID_k;
        ID_n;
        
        RM; %whether the link is ramp metered
    end
    
    methods
        
        function si = setConstraints_noVSL(L, chi, xi, v, w, km, kc, sizes, size_row, network_junc, num_steps, T, junc_indicator, ID_x,ID_k,ID_n,RM)
            si.L = L;
            si.chi = chi;
            si.xi = xi;
            si.v = v;
            si.w = w;
            si.km = km;
            si.kc = kc; 
            
            si.junc_indicator = junc_indicator;
            si.network_junc = network_junc;
            
            si.ID_x = ID_x;
            si.ID_k = ID_k;
            si.ID_n = ID_n;
            
            si.RM = RM; %0 or 1
            
            si.sizes = sizes;
            si.size_row = size_row;
            [si.nb_upstream, si.nb_downstream] = si.num_bin(num_steps, T);

            
            new_Obj_or_iRM = 3;
%             si.sizes = [si.sizes si.junc_indicator(1)*(num_steps+1) si.junc_indicator(2)*(num_steps+1) sum(si.nb_downstream) sum(si.nb_upstream)]; %add number of Dem,Sup,nb_upstr,nb_downstr
            si.sizes = [si.sizes(1:new_Obj_or_iRM-1) (num_steps+1) (num_steps+1) sum(si.nb_downstream) sum(si.nb_upstream) si.sizes(new_Obj_or_iRM:end)]; %add number of Dem,Sup,nb_upstr,nb_downstr
%             si.size_row = sum(si.sizes)+1;
                                 
            si.size_row = sum(si.sizes)+1;

        end
        
        function [arrayineq] = mgamma0(si, p, n, x, T)
            
            qUS = 1; %upstream flows
            qDS = 2; %downstream flows

            arrayineq = zeros(1, si.size_row);
            i=1;
            t = p*T;
            if(x==si.xi) %upstream points

                if(t < (n+1)*T)
                    arrayineq = [];
                    return
                elseif(t >= (n+1)*T)
                    % M >= c
                    T_t = t - (n+1)*T;
                    arrayineq(i , sum(si.sizes(1:qUS-1))+n+2:sum(si.sizes(1:qUS-1))+p) = T;
                    arrayineq(i , si.size_row) = si.kc*si.v*T_t;
                    return
                end
                
            elseif(x == si.chi) %downstream points

                if(t < n*T + si.L/(si.v))
                    arrayineq = [];
                    return
                elseif(t >= n*T + si.L/(si.v) && t <= (n+1)*T + si.L/(si.v))
                    T_t = t - (n*T + si.L/(si.v));
                    arrayineq(i , sum(si.sizes(1:qDS-1))+1:sum(si.sizes(1:qDS-1))+p) = T;
                    arrayineq(i , si.size_row) = -si.ID_n(end);
                    arrayineq(i , sum(si.sizes(1:qUS-1))+1:sum(si.sizes(1:qUS-1))+n) = -T;
                    arrayineq(i , sum(si.sizes(1:qUS-1))+n+1) = -T_t;
                    return
                elseif(t > (n+1)*T + si.L/(si.v))
                    T_t = t - (n+1)*T;
                    arrayineq(i , sum(si.sizes(1:qDS-1))+1:sum(si.sizes(1:qDS-1))+p) = T;
                    arrayineq(i , si.size_row) = -si.ID_n(end) + si.kc*si.v*T_t - si.kc*(x-si.xi);
                    arrayineq(i , sum(si.sizes(1:qUS-1))+1:sum(si.sizes(1:qUS-1))+n+1) = -T;
                    return
                end
                
            end
                
        end
        
        function [arrayineq] = mbeta0(si, p, n, x, T)
            qUS = 1; %upstream flows
            qDS = 2; %downstream flows
            
            arrayineq = zeros(1 , si.size_row);
            i=1; %replace this with 1 directly
            t=p*T;
            if(x == si.xi) %upstream points

                if(t < n*T - si.L/si.w)
                    arrayineq = [];
                    return
                elseif(t >= n*T - si.L/si.w && t <= (n+1)*T - si.L/si.w)
                    T_t = t - (n*T - si.L/si.w);
                    arrayineq(i , sum(si.sizes(1:qUS-1))+1:sum(si.sizes(1:qUS-1))+p) = T; %c
                    arrayineq(i , sum(si.sizes(1:qDS-1))+1:sum(si.sizes(1:qDS-1))+n) = -T;
                    arrayineq(i , si.size_row) = si.ID_n(end) + si.kc*si.v*(si.chi-x)/(-si.w) + si.kc*(si.chi-x);
                    arrayineq(i , sum(si.sizes(1:qDS-1))+n+1) = -T_t;
                    return
                elseif(t > (n+1)*T - si.L/si.w)
                    T_t = t - (n+1)*T;
                    arrayineq(i , sum(si.sizes(1:qUS-1))+1:sum(si.sizes(1:qUS-1))+p) = T; %c
                    arrayineq(i , sum(si.sizes(1:qDS-1))+1:sum(si.sizes(1:qDS-1))+n+1) = -T;
                    arrayineq(i , si.size_row) = si.ID_n(end) + si.kc*si.v*T_t + si.kc*(si.chi-x);
                    return
                end
                
            elseif(x == si.chi) %downstream points

                if (t < (n+1)*T)
                    arrayineq = []; 
                    return
                elseif(t >= (n+1)*T)
                    T_t = t - (n+1)*T;
                    arrayineq(i , sum(si.sizes(1:qDS-1))+n+2 : sum(si.sizes(1:qDS-1))+p) = T;
                    arrayineq(i , si.size_row) = si.kc*si.v*T_t;
                    return
                end
                
            end
            
        end
        
        function [arrayineq] = mtau0(si, p, k, x, T)
            qUS = 1; %upstream flows
            qDS = 2; %downstream flows
            
            arrayineq = zeros(1, si.size_row);
            i=0;

            t = p*T;

            int1 = si.ID_x(k+1); %top
            int2 = si.ID_x(k); %bottom

            if(x == si.xi)% Upstream Points

                if(t == 0 && si.ID_x(k) > si.xi)
                    arrayineq = [];
                    return
                elseif(t ==0 && si.ID_x(k) == si.xi)
                    return;
                end
%                     arrayineq =[];return
                %possible point of projection for p
                pru = si.xi - si.w*t;

                if(pru < int2)
                    arrayineq = [];
                    return;

                elseif(pru >= int2 && pru <= int1)

                    % M1 = N_int3 + k_c*t*(v - int3/t)
                    i=i+1;
                    arrayineq(i , sum(si.sizes(1:qUS-1))+1:sum(si.sizes(1:qUS-1))+p) = T;
                    arrayineq(i , si.size_row) = si.ID_n(k) + si.kc*t*(si.v + (int2-si.xi)/t); %N_int3

                    % M2 = N_int3 - k*(pru-int3) + k_c*t*(v - w)
                    i=i+1;
                    arrayineq(i , sum(si.sizes(1:qUS-1))+1:sum(si.sizes(1:qUS-1))+p) = T;
                    arrayineq(i , si.size_row) = si.ID_n(k) - si.ID_k(k)*(pru - int2) + si.kc*t*(si.v - si.w); %N_int3
                    return

                elseif(pru > int1)

                    % M1 = N_int1 + k_c*t*(v - int1/t)
                    i=i+1;
                    arrayineq(i , sum(si.sizes(1:qUS-1))+1:sum(si.sizes(1:qUS-1))+p) = T; %c
                    arrayineq(i , si.size_row) = si.ID_n(k+1) + si.kc*t*(si.v + (int1-si.xi)/t);  %N_int1
                    
                    i=i+1;
                    arrayineq(i , sum(si.sizes(1:qUS-1))+1:sum(si.sizes(1:qUS-1))+p) = T; %c
                    arrayineq(i , si.size_row) = si.ID_n(k) + si.kc*t*(si.v + (int2-si.xi)/t);  %N_int2
                    return

                end

            elseif(x == si.chi)%Downstream Points

                if(t == 0 && si.ID_x(k+1) < si.chi)
                    arrayineq = [];
                    return
                elseif(t == 0 && si.ID_x(k+1) == si.chi)
                    return
                end
%                     arrayineq =[];return
                %prd < xi; %because v > L/T therefore irrlevant

                prd = si.chi - si.v*t;
                
                if(prd > int1)
                    arrayineq = [];
                    return
                elseif(prd <= int1 && prd >= int2)
%                     i=i+1;
                    % M1 = N_int3 + k_c*t*(v - int3/t)
                    i=i+1;
                    arrayineq(i , sum(si.sizes(1:qDS-1))+1:sum(si.sizes(1:qDS-1))+p) = T;
                    arrayineq(i , si.size_row) = si.ID_n(k+1) - si.ID_n(end) + si.kc*t*(si.v - (si.chi-int1)/t); %N_int3

                    % M2 = N_int3 - k*(pru-int3) + k_c*t*(v - v)
                    i=i+1;
                    arrayineq(i , sum(si.sizes(1:qDS-1))+1:sum(si.sizes(1:qDS-1))+p) = T;
                    arrayineq(i , si.size_row) = si.ID_n(k+1) - si.ID_n(end) + si.ID_k(k)*(int1 - prd); %N_int3
                    return
                    
                elseif(prd < int2)
                    % M1 = N_int1 + k_c*t*(v - (chi-int1)/t)
                    i=i+1;
                    arrayineq(i , sum(si.sizes(1:qDS-1))+1:sum(si.sizes(1:qDS-1))+p) = T;
                    arrayineq(i , si.size_row) = si.ID_n(k+1) - si.ID_n(end) + si.kc*t*(si.v - (si.chi-int1)/t); %N_int1

                    % M2 = N_int2 + k_c*t*(v - (si.chi-int2)/t)
                    i=i+1;
                    arrayineq(i , sum(si.sizes(1:qDS-1))+1:sum(si.sizes(1:qDS-1))+p) = T;
                    arrayineq(i , si.size_row) = si.ID_n(k) - si.ID_n(end) + si.kc*t*(si.v - (si.chi-int2)/t); %N_int3
                    return
                end
            end            
        
            
        end
        
        function [arrayineq] = setAllConstraints(si, num_steps, T)
            
            %Defining index variables
            qUS = 1; %upstream flows
            qDS = 2; %downstream flows
            iRM = 7; %rampmetering
            RMpenalty = 8; %Ramp metering penalty
            if(si.RM == 1)
                new_Obj = 9;
            else
                new_Obj = 7;
            end
            
            i = 0;
            arrayineq = zeros(0,0);
            %% NCalc_noswitch

            for p = 1 : num_steps %upstream points

                % from Upstream Condition
                if(p >= 2)
                    n = p-2;
                    array = si.mgamma0(p, n, si.xi, T);
                    arrayineq = [arrayineq; array];
                    i = i + size(array,1);
                end
                
                %from Downstream Condition
                if(p >= ceil(si.L/(-si.w*T)))
                    n = p - ceil(si.L/(-si.w*T));
                    array = si.mbeta0(p, n, si.xi, T);
                    arrayineq = [arrayineq; array];
                    i = i + size(array,1);
                end
                
            end

            for p = 1 : num_steps %downstream points

                % from Upstream Condition
                if(p >= ceil(si.L/(si.v*T)))
                    n = p-ceil(si.L/(si.v*T));
                    array = si.mgamma0(p, n, si.chi, T);
                    arrayineq = [arrayineq; array];
                    i = i + size(array,1);
                end
                
                % from Downstream Condition
                if(p >= 2)
                    n = p-2;
                    array = si.mbeta0(p, n, si.chi, T);
                    arrayineq = [arrayineq; array];
                    i = i + size(array,1);
                end
            end


            %% NIniDen
            
            for k = 1 : length(si.ID_k)
                %upstream points
                for p = 1 : ceil(si.L/(-si.w*T))
                    array = si.mtau0(p, k, si.xi, T);
                    arrayineq = [arrayineq; array];
                    i = i + size(array,1);
                end
                
                %downstream points
                for p = 1 : ceil(si.L/(si.v*T))
                    array = si.mtau0(p, k, si.chi, T);
                    arrayineq = [arrayineq; array];
                    i = i + size(array,1);
                end
            end                        

            
            %% RM
            if(si.RM == 1)
               
               for n = 1 : num_steps
                  i=i+1;
                  arrayineq(i , sum(si.sizes(1:qDS-1))+n) = 1;
                  arrayineq(i , sum(si.sizes(1:iRM-1))+n) = -si.kc*si.v;
               end
                
            end            
            
            %% RM Penalty (Safety)
            if(si.RM==1)
                
                for sec = 1 : num_steps-1
                   i=i+1;
                   arrayineq(i , sum(si.sizes(1:RMpenalty-1))+sec) = -1;
                   arrayineq(i , sum(si.sizes(1:iRM-1))+sec) = -1;
                   arrayineq(i , sum(si.sizes(1:iRM-1))+(sec+1)) = 1;
                   
                   i=i+1;                   
                   arrayineq(i, sum(si.sizes(1:RMpenalty-1))+sec) = -1;
                   arrayineq(i , sum(si.sizes(1:iRM-1))+sec) = 1;
                   arrayineq(i , sum(si.sizes(1:iRM-1))+(sec+1)) = -1;
                   
                end
            end
            
            
            %% new_Obj
            for n = 1 : num_steps
                i = i + 1;
                arrayineq(i , sum(si.sizes(1:new_Obj-1))+n) = -1;
                arrayineq(i , sum(si.sizes(1:qUS-1))+1:sum(si.sizes(1:qUS-1))+n) = T/si.L;
                arrayineq(i , sum(si.sizes(1:qDS-1))+1:sum(si.sizes(1:qDS-1))+n) = -T/si.L;
                arrayineq(i , si.size_row) = si.ID_n(end)/si.L;
                i = i + 1;
                arrayineq(i , sum(si.sizes(1:new_Obj-1))+n) = -1;
                arrayineq(i , si.size_row) = -si.kc;
            end
            
        end
                
        function [list3] = setBoundaryDS(si, NData, num_steps, T, flag)
            
            qUS = 1;
            qDS = 2;
            Dem = 3;
            Sup = 4;
            bDem = 5;
            bSup = 6;
            
            list3 = zeros(0,0);
            
            %Binary variables useful counter as a reference
            countB=0;
            Cmax = 20000; 
            
            % nb: number of binary variables per point
            nb = ceil(log2(2*3*(num_steps)));
            
            % comb_mat: combination matrix for possible binary combinations
            comb_mat = zeros(2^nb,nb);
%             epsilon = 0.00001;
            epsilon = 0;
            
            for i = 1:2^nb
                comb_mat(i,:) = double(dec2bin(i-1,nb)=='1');
            end
            
            nb_upstream = 0;
            nb_downstream = 0;
            
            % when link is a network input - here we are constraining
            % supply
            if (si.junc_indicator(2) == 0 && flag == 1)
                
                for k = 0 : num_steps
                        
                    tempM = zeros(0,0);
                    %Upper bound constraints
                    %L1<=M1
                    %L1<=M2

                    temp_array = zeros(1,si.size_row);
                    temp_array(1,sum(si.sizes(1:Sup-1))+k+1) = 1;
                    temp_array(1, si.size_row) = epsilon;
                    countM = 0;
                    
                   if(k <= ceil(si.L/(-si.w*T)))
                        for id = 1 : length(si.ID_k)

                           arraym = si.mtau0(k, id, si.xi, T);
                           if(~isempty(arraym))
                                arrayq = zeros(1, si.size_row);
                                arrayq(1 , sum(si.sizes(1:qUS-1))+1:sum(si.sizes(1:qUS-1))+k) = T;
%                                     temp_array(1, si.size_row) = 0;
                                array2 = arraym - arrayq + temp_array;
                                %Add a constraint
                                rows = size(list3,1);
                                list3(rows+1:rows+(size(array2,1)),:) = array2;
                                %Assign the solution to temporary Matrix
                                rows = size(tempM,1);
                                tempM(rows+1:rows+(size(array2,1)),:) = arraym-arrayq;
                                countM = countM + (size(array2,1));
                            end                         
                        end
                   end

                    if (k-2 >= 0)    % For upstream in same switch
                        n = k-2;

                        arraym = si.mgamma0(k, n, si.xi, T);

                        if(~isempty(arraym))
                            arrayq = zeros(1, si.size_row);
                            arrayq(1 , sum(si.sizes(1:qUS-1))+1:sum(si.sizes(1:qUS-1))+k) = T;
%                                     temp_array(1, si.size_row) = 0;
                            array2 = arraym - arrayq + temp_array;
                            %Add a constraint
                            rows = size(list3,1);
                            list3(rows+1:rows+(size(array2,1)),:) = array2;
                            %Assign the solution to temporary Matrix
                            rows = size(tempM,1);
                            tempM(rows+1:rows+(size(array2,1)),:) = arraym - arrayq;
                            countM = countM + (size(array2,1));
                        end
                    end

                    %for downstream conditions in same section

                    n = k - ceil(si.L/(-si.w*T));
                    if(n >= 0)
                        arraym = si.mbeta0(k, n, si.xi, T);
                        if(~isempty(arraym))
                            arrayq = zeros(1, si.size_row);
                            arrayq(1 , sum(si.sizes(1:qUS-1))+1:sum(si.sizes(1:qUS-1))+k) = T;
    %                                 temp_array(1, si.size_row) = 0;
                            array2 = arraym - arrayq + temp_array;
                            %Add a constraint
                            rows = size(list3,1);
                            list3(rows+1:rows+(size(array2,1)),:) = array2;
                            %Assign the solution to temporary Matrix
                            rows = size(tempM,1);
                            tempM(rows+1:rows+(size(array2,1)),:) = arraym-arrayq;
                            countM =countM + (size(array2,1));
                        end
                    end
                    %for Supply constraints from data
                    rows = size(list3,1);
                    arraym = zeros(1, si.size_row);
                    arraym(1, si.size_row) = NData(k+1);
                    list3(rows+1, :) = arraym + temp_array;
                    rows = size(tempM,1);
                    tempM(rows+1, :) = arraym; 
                    countM = countM + 1;
                    
                    
                    rowsM = size(tempM,1);

                    nb_upstream = ceil(log2(countM));
                    %Define all the possible combinations according to the solutions that apply
                    %at the specified point

                    %Lower bound constraints
                    %C*b1 +L1>=M1
                    %C*b1 +L1>=M2
                    for i = 1 : rowsM

                        array = temp_array;

                        %Decode the binary combinations

                        for counter = 1 : si.nb_upstream(k+1)
                            if comb_mat(i,nb-si.nb_upstream(k+1)+counter) == 1
                                array(1,sum(si.sizes(1:bSup-1))+ countB+ counter) = -Cmax;

                            elseif comb_mat(i,nb-si.nb_upstream(k+1)+counter) == 0
                                array(1,sum(si.sizes(1:bSup-1))+ countB+ counter) = +Cmax;

                            end
                        end
                        array(1,si.size_row) = - Cmax*(sum(comb_mat(i,:)))-epsilon; 

                        % Add constraint to the MILP matrix

                        array2 = - array - tempM(i,:) ;
                        rows = size(list3,1);
                        list3(rows+1,:) = array2;
                    end

                    %Define the last constraint to bound the binary combination

                    array = zeros(1,si.size_row);

                    for counter = 1 : si.nb_upstream(k+1)

                        array(1,sum(si.sizes(1:bSup-1)) + countB + counter) = 2^(si.nb_upstream(k+1) - counter);

                    end

                    array(1,si.size_row) = (rowsM-1); %RHS (maximum possible value of the binary comb)
                    rows = size(list3,1);
                    list3(rows+1,:) = array;

                    countB = countB + si.nb_upstream(k+1);

                end
            end   
            
            
            % when link is a network output - here we are constraining
            % demand
            if(si.junc_indicator(1) == 0 && flag == 2)
                
                countB = 0;
                
                for k = 0 : num_steps

                    tempM = zeros(0,0);

                    %Upper bound constraints
                    %L1<=M1
                    %L1<=M2

                    temp_array = zeros(1,si.size_row);
                    temp_array(1,sum(si.sizes(1:Dem-1))+k+1)=1;
                    temp_array(1, si.size_row) = epsilon;
                    countM = 0;

                    if(k <= ceil(si.L/(si.v*T)))
                        for id = 1 : length(si.ID_k)

                           arraym = si.mtau0(k, id, si.chi, T);
                           if(~isempty(arraym))
                                arrayq = zeros(1, si.size_row);
                                arrayq(1 , sum(si.sizes(1:qDS-1))+1:sum(si.sizes(1:qDS-1))+k) = T;
                                arrayq(1 , si.size_row) = -si.ID_n(end);
%                                 temp_array(1, si.size_row) = 0;
                                array2 = arraym - arrayq + temp_array;
                                %Add a constraint
                                rows = size(list3,1);
                                list3(rows+1:rows+(size(array2,1)),:) = array2;
                                %Assign the solution to temporary Matrix
                                rows = size(tempM,1);
                                tempM(rows+1:rows+(size(array2,1)),:) = arraym-arrayq;
                                countM = countM + (size(array2,1));
                            end                         
                        end
                    end

                    if (k-2 >= 0)
                        n=k-2;
                        arraym = si.mbeta0(k, n, si.chi, T);
                        if(~isempty(arraym))
                            arrayq = zeros(1, si.size_row);
                            arrayq(1 , sum(si.sizes(1:qDS-1))+1:sum(si.sizes(1:qDS-1))+k) = T;
                            arrayq(1 , si.size_row) = -si.ID_n(end);
%                                  temp_array(1, si.size_row) = 0;
                            array2 = arraym - arrayq + temp_array;
                            %Add a constraint
                            rows = size(list3,1);
                            list3(rows+1:rows+(size(array2,1)),:) = array2;
                            %Assign the solution to temporary Matrix
                            rows = size(tempM,1);
                            tempM(rows+1:rows+(size(array2,1)),:) = arraym-arrayq;
                            countM = countM + (size(array2,1));
                        end
                    end
                    
                    %for upstream conditions in same section

                    n = k - ceil(si.L/(si.v*T));
                    if(n >= 0)
                        arraym = si.mgamma0(k, n, si.chi, T);
                        if(~isempty(arraym))
                            arrayq = zeros(1, si.size_row);
                            arrayq(1 , sum(si.sizes(1:qDS-1))+1:sum(si.sizes(1:qDS-1))+k) = T;
                            arrayq(1 , si.size_row) = -si.ID_n(end);
    %                             temp_array(1, si.size_row) = 0;
                            array2 = arraym - arrayq + temp_array;
                            %Add a constraint
                            rows = size(list3,1);
                            list3(rows+1:rows+(size(array2,1)),:) = array2;
                            %Assign the solution to temporary Matrix
                            rows = size(tempM,1);
                            tempM(rows+1:rows+(size(array2,1)),:) = arraym-arrayq;
                            countM = countM + (size(array2,1));
                        end
                    end                    
                    %for Demand constraints from data
                    rows = size(list3,1);
                    arraym = zeros(1, si.size_row);
                    arraym(1, si.size_row) = NData(k+1);
                    list3(rows+1, :) = arraym + temp_array;
                    rows = size(tempM,1);
                    tempM(rows+1, :) = arraym; 
                    countM = countM + 1;

                    rowsM = size(tempM,1);

                    nb_downstream = ceil(log2(countM));
                    %Define all the possible combinations according to the solutions that apply
                    %at the specified point

                    %Lower bound constraints
                    %C*b1 +L1>=M1
                    %C*b1 +L1>=M2

                    for i = 1:rowsM

                        array = temp_array;
                        %Decode the binary combinations
                        for counter=1:si.nb_downstream(k+1)
                            if comb_mat(i,nb-si.nb_downstream(k+1)+counter) == 1
                                array(1,sum(si.sizes(1:bDem-1))+ countB + counter) = -Cmax;

                            elseif comb_mat(i,nb-si.nb_downstream(k+1)+counter) == 0
                                
                                array(1,sum(si.sizes(1:bDem-1))+  countB + counter) = +Cmax;

                            end
                        end
                        array(1,si.size_row) = - Cmax*(sum(comb_mat(i,:)))-epsilon; %RHS (set negative to be on same side)

                        % Add constraint to the matrix constraints

                        array2 = - array - tempM(i,:);
                        rows = size(list3,1);
                        list3(rows+1,:) = array2;
                    end

                    %Define the last constraint to bound the binary combination

                    array = zeros(1,si.size_row);

                    % Define constraint matrix elements for binary terms
                    for counter=1:si.nb_downstream(k+1)

                        array(1,sum(si.sizes(1:bDem-1))+ countB + counter) = 2^(si.nb_downstream(k+1)-counter);

                    end

                    array(1,si.size_row) = (rowsM-1); %RHS
                    rows = size(list3,1);
                    list3(rows+1,:) = array;

                    countB = countB + si.nb_downstream(k+1);

                end  

                
                
            end
        end
        
        function [list] = setMatrixDS(si, num_steps, T, flag)    
            qUS = 1; %upstream flows
            qDS = 2; %downstream flows
            Dem = 3; %demand variables
            Sup = 4; %supply variables
            
            list = zeros(0,si.size_row);
            
            %Binary variables useful counter as a reference
            countB=0;
            Cmax = 20000; 
            
            % nb: number of binary variables per point THIS CAN BE REDUCED
            nb = ceil(log2(10));
            
            % comb_mat: combination matrix for possible binary combinations
            comb_mat = zeros(2^nb,nb);
%             epsilon = 0.0000001;
            epsilon = 0;
            
            for i = 1:2^nb
                comb_mat(i,:) = double(dec2bin(i-1,nb)=='1');
            end
            
            
            if(si.junc_indicator(1)==0 && flag == 2) % Demand (Downstream Flows)
                  
                for ti = 1 : num_steps

                    tempM = zeros(0,0);
                    tempArray = zeros(1,si.size_row);

%                             linkStr = sprintf('link_%d',si.network_junc.(juncStr).inlabel(1,in));
%                             tmp_offset = si.network_links.(linkStr).offset+...
%                                 sum( si.network_links.(linkStr).sizes(1:qDS-1) );

                    tmp_offset = sum( si.sizes(1:qDS-1) );

                    tempArray(1, tmp_offset + ti)...
                        = 1;

                    tempArray(1, si.size_row) = epsilon;

                    tmp_offset = sum( si.sizes(1:Dem-1) );

                    arraym = zeros(1,si.size_row);

                    arraym(1, tmp_offset + 1)= -1/T;
                    %arraym(1, tmp_offset + ti)= -1/T;
                    arraym(1, tmp_offset + ti+1)= 1/T; 

                    if(ti > 1)
                        tmp_offset = sum( si.sizes(1:qDS-1) );
                        arraym(1, tmp_offset + 1: tmp_offset + ti-1)= -1;
                    end

                    array2 = tempArray - arraym;
                    if(~isempty(array2))
                        rows = size(list,1);
                        list(rows+1,:) = array2;
                        rows = size(tempM,1);
                        tempM(rows+1,:) = arraym;

                    end


                end

            end  
            
            
            if(si.junc_indicator(2)==0 && flag == 1) % Supply (Upstream Flows)
                  
                for ti = 1 : num_steps

                    tempM = zeros(0,0);
                    tempArray = zeros(1,si.size_row);

%                             linkStr = sprintf('link_%d',si.network_junc.(juncStr).inlabel(1,in));
%                             tmp_offset = si.network_links.(linkStr).offset+...
%                                 sum( si.network_links.(linkStr).sizes(1:qDS-1) );

                    tmp_offset = sum( si.sizes(1:qUS-1) );

                    tempArray(1, tmp_offset + ti)...
                        = 1;

                    tempArray(1, si.size_row) = epsilon;

                    tmp_offset = sum( si.sizes(1:Sup-1) );

                    arraym = zeros(1,si.size_row);

                    arraym(1, tmp_offset + 1)= -1/T;
                    %arraym(1, tmp_offset + ti)= -1/T;
                    arraym(1, tmp_offset + ti+1)= 1/T; 

                    if(ti > 1)
                        tmp_offset = sum( si.sizes(1:qUS-1) );
                        arraym(1, tmp_offset + 1: tmp_offset + ti-1)= -1;
                    end

                    array2 = tempArray - arraym;
                    if(~isempty(array2))
                        rows = size(list,1);
                        list(rows+1,:) = array2;
                        rows = size(tempM,1);
                        tempM(rows+1,:) = arraym;

                    end

                end

            end          
            
        end                                                                                
        
        function [nb_upstream, nb_downstream] = num_bin(si,num_steps,T)
            
            nb_upstream = zeros(1,num_steps+1);
            nb_downstream = zeros(1,num_steps+1);
           
            for k = 0 : num_steps  
    
                   %Supply variables ( Upstream, since the road is an output of a junction)------------------
%                    if (si.junc_indicator(2)==1)     

                       countM = 0;

                        %Upper bound constraints
                        %L1<=M1
                        %L1<=M2
                        
                       
                        if(k <= ceil(si.L/(-si.w*T)))
                            for id = 1 : length(si.ID_k)

                               arraym = si.mtau0(k, id, si.xi, T);
                               if(~isempty(arraym))
                                   countM = countM + size(arraym,1);
                               end                            
                            end
                        end
                        
                        if (k-2 >= 0)    
                            n = k-2;
                            
                            arraym = si.mgamma0(k, n, si.xi, T);
                            if(~isempty(arraym))
                                countM = countM + size(arraym,1);
                            end

                        end
                        
                        %for downstream conditions in same section

                        n = k - ceil(si.L/(-si.w*T));
                        if(n >= 0)
                            arraym = si.mbeta0(k, n, si.xi, T);
                            if(~isempty(arraym))
                                countM = countM + size(arraym,1);
                            end
                        end                        
                        if(si.junc_indicator(2) == 0)
                            countM = countM + 1;
                        end
                        
                        nb_upstream(1,k+1) = ceil(log2(countM));
%                    end



             %Demand (Downstream, since is in an input to the junction)-------------------------
%                if(si.junc_indicator(1)==1)
                    countM = 0;
                    
                    %Upper bound constraints
                    %L1<=M1
                    %L1<=M2
                    if(k <= ceil(si.L/(si.v*T)))
                        for id = 1 : length(si.ID_k)
                           arraym = si.mtau0(k, id, si.chi, T);
                           if(~isempty(arraym))
                               countM = countM + size(arraym,1);
                           end                            
                        end
                    end
                        
                    if (k-2>=0)
                        n=k-2;
                        
                        arraym = si.mbeta0(k, n, si.chi, T);

                        if(~isempty(arraym))
                            countM = countM + size(arraym,1);
                        end
                        
                    end
                    %for upstream conditions in same section

                    n = k - ceil(si.L/(si.v*T));
                    if(n >= 0)
                        arraym = si.mgamma0(k, n, si.chi, T);
                        if(~isempty(arraym))
                            countM = countM + size(arraym,1);
                        end
                    end
                    
                    if(si.junc_indicator(1) == 0)
                        countM = countM + 1;
                    end
                    
                    if(si.RM == 1 && k > 0)
                        countM = countM + 1;    
                    end                                        

                    nb_downstream(1,k+1) = ceil(log2(countM));
               
%                end
           end
        end
                        
        function [list3] = setDemandSupply(si,num_steps,T)
            qUS = 1;
            qDS = 2;
            Dem = 3;
            Sup = 4;
            bDem = 5;
            bSup = 6;
            iRM = 7;
            
            list3 = zeros(0,0);
            
            %Binary variables useful counter as a reference
            countB=0;
            Cmax = 20000; 
            
            % nb: number of binary variables per point
            nb = ceil(log2(2*3*(num_steps)));
            
            % comb_mat: combination matrix for possible binary combinations
            comb_mat = zeros(2^nb,nb);
            epsilon = 0;
            for i = 1:2^nb
                comb_mat(i,:) = double(dec2bin(i-1,nb)=='1');
            end
            
            nb_upstream = 0;
            nb_downstream = 0;
        
           %Suppply variables ( Upstream, since the road is an output of a junction)------------------
             if (si.junc_indicator(2)==1)     

                    for k = 0 : num_steps
                        
                        tempM = zeros(0,0);
                        %Upper bound constraints
                        %L1<=M1
                        %L1<=M2

                        temp_array = zeros(1,si.size_row);
                        temp_array(1,sum(si.sizes(1:Sup-1))+k+1) = 1;
                        temp_array(1, si.size_row) = epsilon;
                        countM = 0;
                       if(k <= ceil(si.L/(-si.w*T)))
                            for id = 1 : length(si.ID_k)

                               arraym = si.mtau0(k, id, si.xi, T);
                               if(~isempty(arraym))
                                    arrayq = zeros(1, si.size_row);
                                    arrayq(1 , sum(si.sizes(1:qUS-1))+1:sum(si.sizes(1:qUS-1))+k) = T;
%                                     temp_array(1, si.size_row) = 0;
                                    array2 = arraym - arrayq + temp_array;
                                    %Add a constraint
                                    rows = size(list3,1);
                                    list3(rows+1:rows+(size(array2,1)),:) = array2;
                                    %Assign the solution to temporary Matrix
                                    rows = size(tempM,1);
                                    tempM(rows+1:rows+(size(array2,1)),:) = arraym-arrayq;
                                    countM = countM + (size(array2,1));
                                end                         
                            end
                       end

                        if (k-2 >= 0)    % For upstream in same switch
                            n = k-2;
                            
                            arraym = si.mgamma0(k, n, si.xi, T);

                            if(~isempty(arraym))
                                arrayq = zeros(1, si.size_row);
                                arrayq(1 , sum(si.sizes(1:qUS-1))+1:sum(si.sizes(1:qUS-1))+k) = T;
%                                     temp_array(1, si.size_row) = 0;
                                array2 = arraym - arrayq + temp_array;
                                %Add a constraint
                                rows = size(list3,1);
                                list3(rows+1:rows+(size(array2,1)),:) = array2;
                                %Assign the solution to temporary Matrix
                                rows = size(tempM,1);
                                tempM(rows+1:rows+(size(array2,1)),:) = arraym - arrayq;
                                countM = countM + (size(array2,1));
                            end

                        end
                        
                        %for downstream conditions in same section

                        n = k - ceil(si.L/(-si.w*T));
                        if(n >= 0)
                            arraym = si.mbeta0(k, n, si.xi, T);
                            if(~isempty(arraym))
                                arrayq = zeros(1, si.size_row);
                                arrayq(1 , sum(si.sizes(1:qUS-1))+1:sum(si.sizes(1:qUS-1))+k) = T;
    %                                 temp_array(1, si.size_row) = 0;
                                array2 = arraym - arrayq + temp_array;
                                %Add a constraint
                                rows = size(list3,1);
                                list3(rows+1:rows+(size(array2,1)),:) = array2;
                                %Assign the solution to temporary Matrix
                                rows = size(tempM,1);
                                tempM(rows+1:rows+(size(array2,1)),:) = arraym-arrayq;
                                countM =countM + (size(array2,1));
                            end
                        end                        
                        rowsM = size(tempM,1);

                        nb_upstream = ceil(log2(countM));
                        %Define all the possible combinations according to the solutions that apply
                        %at the specified point

                        %Lower bound constraints
                        %C*b1 +L1>=M1
                        %C*b1 +L1>=M2
                        for i = 1 : rowsM

                            array = temp_array;

                            %Decode the binary combinations

                            for counter = 1 : si.nb_upstream(k+1)
                                if comb_mat(i,nb-si.nb_upstream(k+1)+counter) == 1
                                    array(1,sum(si.sizes(1:bSup-1))+ countB+ counter) = -Cmax;

                                elseif comb_mat(i,nb-si.nb_upstream(k+1)+counter) == 0
                                    array(1,sum(si.sizes(1:bSup-1))+ countB+ counter) = +Cmax;

                                end
                            end
                            array(1,si.size_row) = - Cmax*(sum(comb_mat(i,:)))-epsilon; 

                            % Add constraint to the MILP matrix

                            array2 = - array - tempM(i,:) ;
                            rows = size(list3,1);
                            list3(rows+1,:) = array2;
                        end

                        %Define the last constraint to bound the binary combination

                        array = zeros(1,si.size_row);

                        for counter = 1 : si.nb_upstream(k+1)

                            array(1,sum(si.sizes(1:bSup-1)) + countB + counter) = 2^(si.nb_upstream(k+1) - counter);

                        end

                        array(1,si.size_row) = (rowsM-1); %RHS (maximum possible value of the binary comb)
                        rows = size(list3,1);
                        list3(rows+1,:) = array;

                        countB = countB + si.nb_upstream(k+1);

                    end

             end


             %Demand (Downstream, since is in an input to the junction)-------------------------
           if(si.junc_indicator(1)==1)

                countB = 0;
                for k = 0 : num_steps
                    
                    tempM = zeros(0,0);

                    %Upper bound constraints
                    %L1<=M1
                    %L1<=M2

                    temp_array = zeros(1,si.size_row);
                    temp_array(1,sum(si.sizes(1:Dem-1))+k+1)=1;
                    temp_array(1, si.size_row) = epsilon;
                    countM = 0;

                    if(k <= ceil(si.L/(si.v*T)))
                        for id = 1 : length(si.ID_k)

                           arraym = si.mtau0(k, id, si.chi, T);
                           if(~isempty(arraym))
                                arrayq = zeros(1, si.size_row);
                                arrayq(1 , sum(si.sizes(1:qDS-1))+1:sum(si.sizes(1:qDS-1))+k) = T;
                                arrayq(1 , si.size_row) = -si.ID_n(end);
%                                 temp_array(1, si.size_row) = 0;
                                array2 = arraym - arrayq + temp_array;
                                %Add a constraint
                                rows = size(list3,1);
                                list3(rows+1:rows+(size(array2,1)),:) = array2;
                                %Assign the solution to temporary Matrix
                                rows = size(tempM,1);
                                tempM(rows+1:rows+(size(array2,1)),:) = arraym-arrayq;
                                countM = countM + (size(array2,1));
                            end                         
                        end
                    end

                    if (k-2 >= 0)
                        n=k-2;

                        arraym = si.mbeta0(k, n, si.chi, T);

                        if(~isempty(arraym))
                            arrayq = zeros(1, si.size_row);
                            arrayq(1 , sum(si.sizes(1:qDS-1))+1:sum(si.sizes(1:qDS-1))+k) = T;
                            arrayq(1 , si.size_row) = -si.ID_n(end);
%                                  temp_array(1, si.size_row) = 0;
                            array2 = arraym - arrayq + temp_array;
                            %Add a constraint
                            rows = size(list3,1);
                            list3(rows+1:rows+(size(array2,1)),:) = array2;
                            %Assign the solution to temporary Matrix
                            rows = size(tempM,1);
                            tempM(rows+1:rows+(size(array2,1)),:) = arraym-arrayq;
                            countM = countM + (size(array2,1));
                        end
                    end
                    
                    %for upstream conditions in same section

                    n = k - ceil(si.L/(si.v*T));
                    if(n >=0)
                        arraym = si.mgamma0(k, n, si.chi, T);
                        if(~isempty(arraym))
                            arrayq = zeros(1, si.size_row);
                            arrayq(1 , sum(si.sizes(1:qDS-1))+1:sum(si.sizes(1:qDS-1))+k) = T;
                            arrayq(1 , si.size_row) = -si.ID_n(end);
    %                             temp_array(1, si.size_row) = 0;
                            array2 = arraym - arrayq + temp_array;
                            %Add a constraint
                            rows = size(list3,1);
                            list3(rows+1:rows+(size(array2,1)),:) = array2;
                            %Assign the solution to temporary Matrix
                            rows = size(tempM,1);
                            tempM(rows+1:rows+(size(array2,1)),:) = arraym-arrayq;
                            countM = countM + (size(array2,1));
                        end
                    end
                    
                    %Ramp Metering constraints
                    if(si.RM == 1 && k > 0)
                        arraym = zeros(1,si.size_row);
                        arraym(1, sum(si.sizes(1:qDS-1))+1:sum(si.sizes(1:qDS-1))+k-1) = -T;
                        arraym(1, sum(si.sizes(1:iRM-1))+k) = -si.kc*si.v*T;
                        arraym(1,si.size_row) = si.ID_n(end); 
                        if(~isempty(arraym))
                            array2 = arraym + temp_array;
                            %Add a constraint
                            rows = size(list3,1);
                            list3(rows+1:rows+(size(array2,1)),:) = array2;
                            %Assign the solution to temporary Matrix
                            rows = size(tempM,1);
                            tempM(rows+1:rows+(size(array2,1)),:) = arraym;
                            countM = countM + (size(array2,1));
                        end
                    end
                    
                    rowsM = size(tempM,1);

                    nb_downstream = ceil(log2(countM));
                    %Define all the possible combinations according to the solutions that apply
                    %at the specified point

                    %Lower bound constraints
                    %C*b1 +L1>=M1
                    %C*b1 +L1>=M2

                    for i = 1:rowsM

                        array = temp_array;
                        %Decode the binary combinations
                        for counter=1:si.nb_downstream(k+1)
                            if comb_mat(i,nb-si.nb_downstream(k+1)+counter) == 1
                                array(1,sum(si.sizes(1:bDem-1))+ countB + counter) = -Cmax;

                            elseif comb_mat(i,nb-si.nb_downstream(k+1)+counter) == 0
                                
                                array(1,sum(si.sizes(1:bDem-1))+  countB + counter) = +Cmax;

                            end
                        end
                        array(1,si.size_row) = - Cmax*(sum(comb_mat(i,:)))-epsilon; %RHS (set negative to be on same side)

                        % Add constraint to the matrix constraints

                        array2 = - array - tempM(i,:);
                        rows = size(list3,1);
                        list3(rows+1,:) = array2;
                    end

                    %Define the last constraint to bound the binary combination

                    array = zeros(1,si.size_row);

                    % Define constraint matrix elements for binary terms
                    for counter=1:si.nb_downstream(k+1)

                        array(1,sum(si.sizes(1:bDem-1))+ countB + counter) = 2^(si.nb_downstream(k+1)-counter);

                    end

                    array(1,si.size_row) = (rowsM-1); %RHS
                    rows = size(list3,1);
                    list3(rows+1,:) = array;

                    countB = countB + si.nb_downstream(k+1);

                end  

            end
        end 
          
    end
end
            