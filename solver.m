%Sequence of codes: def_net --> network_struct --> solver
%>>> solver
tic
%% Constraint Equations
Aineq = zeros(0,0);
bineq = zeros(0,1);
Aeq = zeros(0,0);
beq = zeros(0,1);

for link = 1 : num_links
    linkStr = sprintf('link_%d',link);
        
    matrixA123 = [network_links.(linkStr).arrayineq;...
         network_links.(linkStr).DS_matrix];
     
     boundaryDS = [];
     boundaryDSe = [];
     %To incorporate network inputs and outputs
     if(network_links.(linkStr).junc_indicator(2) == 0 && network_links.(linkStr).junc_indicator(1) ~= 0)
         boundaryDS = network_links.(linkStr).input;
         boundaryDSe = network_links.(linkStr).input2;
     end
     if(network_links.(linkStr).junc_indicator(1) == 0 && network_links.(linkStr).junc_indicator(2) ~= 0)
         boundaryDS = network_links.(linkStr).output; 
         boundaryDSe = network_links.(linkStr).output2;
     end
     if(network_links.(linkStr).junc_indicator(1) == 0 && network_links.(linkStr).junc_indicator(2) == 0)
         boundaryDS = network_links.(linkStr).input;
         boundaryDSe = network_links.(linkStr).input2;
         boundaryDS = [boundaryDS; network_links.(linkStr).output]; 
         boundaryDSe = [boundaryDSe; network_links.(linkStr).output2];
     end
     
     if(~isempty(boundaryDS))
         matrixA123 = [matrixA123; boundaryDS];
     end
     
    pre_link_namestr = sprintf('link_%d',link-1);
    A123 = matrixA123(:,1:size(matrixA123,2)-1);
    A123 = [zeros(size(matrixA123,1),network_links.(linkStr).offset),...
        A123];
    A123 = [A123...
        zeros(size(matrixA123,1), network_links.totalLength-size(A123,2))];

    Aineq = [Aineq ; A123];

    rhA123 = matrixA123(:,size(matrixA123,2));
    bineq = [bineq ; rhA123];

    %Equality constraints
    if(any(link==VSL_link))
        matrixB123 = [network_links.(linkStr).arrayeq];
        
        if(~isempty(boundaryDSe))
             matrixB123 = [matrixB123; boundaryDSe];
        end
        
        B123 = matrixB123(:,1:size(matrixB123,2)-1);
        B123 = [zeros(size(matrixB123,1),network_links.(linkStr).offset),...
            B123];
        B123 = [B123...
            zeros(size(matrixB123,1), network_links.totalLength-size(B123,2))];

        Aeq = [Aeq ; B123];

        rhB123 = matrixB123(:,size(matrixB123,2));
        beq = [beq ; rhB123];
    else
        if(~isempty(boundaryDSe))
             matrixB123 = boundaryDSe;
             B123 = matrixB123(:,1:size(matrixB123,2)-1);
            B123 = [zeros(size(matrixB123,1),network_links.(linkStr).offset),...
                B123];
            B123 = [B123...
                zeros(size(matrixB123,1), network_links.totalLength-size(B123,2))];

            Aeq = [Aeq ; B123];

            rhB123 = matrixB123(:,size(matrixB123,2));
            beq = [beq ; rhB123];
        end
    end

end

% Add DS inequalities and Junction Equalities
ADS = Junc_ineq.DS_matrixEq(:,1:size(Junc_ineq.DS_matrixEq,2)-1);
rhADS = Junc_ineq.DS_matrixEq(:,size(Junc_ineq.DS_matrixEq,2));

Aineq = [Aineq ; ADS];
bineq = [bineq ; rhADS];

Aeq = [Aeq; Junc_ineq.Junc_matrix];
beq = [beq; zeros(size(Aeq,1)-size(beq,1),1)];

%% Lower and Upper Bounds and Variable Type
lb = zeros(network_links.totalLength,1);
for link = 1 : num_links
    
    linkStr = sprintf('link_%d',link);
    tmp_offset = network_links.(linkStr).offset;
%     set lower bound
    if(any(link==VSL_link))
        lb(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:qUS-1))+1:...
            network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:qUS))) = 0; 

        lb(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:kUS-1))+1:...
            network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:kUS))) = 0; 

        lb(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:qDS-1))+1:...
            network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:qDS))) = 0; 

        lb(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:ka-1))+1:...
            network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:ka))) = 0; 

        lb(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:bin-1))+1:...
            network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:bin))) = 0; 

%         lb(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:TCD-1))+1:...
%             network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:TCD))) = 0;

        lb(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:Dem-1))+1:...
        network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:Dem))) = -10000; 

        lb(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:Sup-1))+1:...
        network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:Sup))) = -10000; 

        lb(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:nbds-1))+1:...
        network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:nbus))) = 0; 

        lb(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:vslpenalty-1))+1:...
        network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:vslpenalty))) = -100;
    
        lb(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:new_Obj-1))+1:...
        network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:new_Obj))) = 0;
    else
        lb(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:qUS_n-1))+1:...
            network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:qUS_n))) = 0; 

        lb(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:qDS_n-1))+1:...
            network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:qDS_n))) = 0; 

        lb(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:Dem_n-1))+1:...
        network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:Dem_n))) = -10000; 

        lb(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:Sup_n-1))+1:...
        network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:Sup_n))) = -10000; 

        lb(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:nbds_n-1))+1:...
        network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:nbus_n))) = 0; 
    
        if(any(link==RM_link))
            lb(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:iRM-1))+1:...
                network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:iRM))) = 0.000001;
            
            lb(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:RMpenalty-1))+1:...
                 network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:RMpenalty))) = 0;

            
            lb(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:new_Obj_n-1))+1:...
                network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:new_Obj_n))) = 0;
        else
            lb(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:iRM-1))+1:...
                network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:iRM))) = 0;
        end
    end
end

ub = zeros(network_links.totalLength,1);
for link = 1 : num_links
    
    linkStr = sprintf('link_%d',link);
    %set upperower bound
    
    if(any(link==VSL_link))    
        ub(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:qUS-1))+1:...
            network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:qUS))) = 5; 

        ub(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:kUS-1))+1:...
            network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:kUS))) = network_links.(linkStr).km; 

        ub(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:qDS-1))+1:...
            network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:qDS))) = 5; 

        ub(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:ka-1))+1:...
            network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:ka))) = network_links.(linkStr).km; 

        ub(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:bin-1))+1:...
            network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:bin))) = 1; 

%         ub(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:TCD-1))+1:...
%             network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:TCD))) = network_links.(linkStr).km;

        ub(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:Dem-1))+1:...
        network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:Dem))) = 10000; 

        ub(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:Sup-1))+1:...
        network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:Sup))) = 10000; 

        ub(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:nbds-1))+1:...
        network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:nbus))) = 1;

        ub(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:vslpenalty-1))+1:...
        network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:vslpenalty))) = 100; 

        ub(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:new_Obj-1))+1:...
        network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:new_Obj))) = 10;

    else
        ub(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:qUS_n-1))+1:...
            network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:qUS_n))) = 5; 

        ub(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:qDS_n-1))+1:...
            network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:qDS_n))) = 5; 

        ub(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:Dem_n-1))+1:...
            network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:Dem_n))) = 10000; 

        ub(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:Sup_n-1))+1:...
            network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:Sup_n))) = 10000; 

        ub(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:nbds_n-1))+1:...
            network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:nbus_n))) = 1; 

                
        
        if(any(link==RM_link))
            ub(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:iRM-1))+1:...
                network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:iRM))) = 1;
            ub(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:RMpenalty-1))+1:...
                 network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:RMpenalty))) = 1;
            ub(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:new_Obj_n-1))+1:...
                network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:new_Obj_n))) = 10;
        else
            ub(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:iRM-1))+1:...
                network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:iRM))) = 10;
        end
    end
end

clear chtmp
for link = 1 : num_links
    %set ch options
    linkStr = sprintf('link_%d',link);

    if(any(link==VSL_link))
        chtmp(network_links.(linkStr).offset+1:...
            network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:bin-1))) = 'C';

        chtmp(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:bin-1))+1:...
            network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:bin))) = 'B';

        chtmp(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:bin))+1:...
            network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:Sup))) = 'C'; %change TCD to Sup when including Dem/Sup

        chtmp(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:Sup))+1:...
            network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:nbus))) = 'B'; %Binaries for Demand Supply
        % counter example - numsteps = 4, switch at 2, 3 links, 2 1-1 junctions,
        % weight on second link +3 rest -1, no density in objective
        chtmp(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:nbus))+1:...
            network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:vslpenalty))) = 'C';

        chtmp(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:vslpenalty))+1:...
            network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:new_Obj))) = 'C';
    else
        chtmp(network_links.(linkStr).offset+1:...
            network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:nbds_n-1))) = 'C';

        chtmp(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:nbds_n-1))+1:...
            network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:nbus_n))) = 'B'; %Binaries for Demand Supply
        

        % counter example - numsteps = 4, switch at 2, 3 links, 2 1-1 junctions,
        % weight on second link +3 rest -1, no density in objective
        if(any(link==RM_link))
            chtmp(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:iRM-1))+1:...
                network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:iRM))) = 'C';
            chtmp(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:RMpenalty-1))+1:...
                 network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:RMpenalty))) = 'C';
            chtmp(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:new_Obj_n-1))+1:...
                network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:new_Obj_n))) = 'C';
        else
            chtmp(network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:iRM-1))+1:...
                network_links.(linkStr).offset+sum(network_links.(linkStr).sizes(1:iRM))) = 'C';
        end
    end
end 
 
for junc = 1 : num_juncs
   
    juncStr = sprintf('junc_%d',junc);
    
    lb(network_junc.(juncStr).offset + 1:...
        network_junc.(juncStr).offset + network_junc.(juncStr).count) = 0;

    ub(network_junc.(juncStr).offset + 1:...
        network_junc.(juncStr).offset + network_junc.(juncStr).count) = 1;
    
    chtmp(network_junc.(juncStr).offset + 1:...
        network_junc.(juncStr).offset + network_junc.(juncStr).count) = 'B';
end

ch = chtmp;

%% Objective Function

% % Set objective
% % Flow maximization
% f = zeros(network_links.totalLength,1);
% for link = 1 : num_links
%     linkStr = sprintf('link_%d',link);
%     
%     if(any(link==VSL_link))
%         f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:qUS-1))+1:...
%             network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:qUS))) = 1;
%         f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:qDS-1))+1:...
%             network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:qDS))) = 1;
% %             Penalty for speed change
% %         f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:vslpenalty-1))+1:...
% %             network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:vslpenalty))) = -1;
%     else
%         f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:qUS_n-1))+1:...
%             network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:qUS_n))) = 1;
%         f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:qDS_n-1))+1:...
%             network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:qDS_n))) = 1;
%     end
% end
% f = -f;

% f(sum(sizes(1:v-1))+1:sum(sizes(1:v))) = -1;

%maximizing the below commented part is equivalent to just maximizing the
%inflow into the network
% TTS = 0;
% for link = 1 : size(InputID,2)
%     li = InputID(link);
%     linkStr = sprintf('link_%d',li);
%     for t = 1 : num_steps       
%         TTS = TTS + (NDemand(link,t+1) - T*(f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:qUS-1))+t)));    
%     end
% end 

%% TTT
% % Total Travel Time
% f = zeros(network_links.totalLength,1);
% we = 1; % weight on density minimization objective
% we2 = 0.2;
% for link = 1 : num_links
%     
%     linkStr = sprintf('link_%d',link);
%     for n = 1 : num_steps
% %         density
%         for j = 1 : n
%             if(any(link==VSL_link))
%                 f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:qDS-1))+j) =...
%                     f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:qDS-1))+j) - we*T;
%                 
%                 f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:qUS-1))+j) =...
%                     f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:qUS-1))+j) + we*T;
%                     
%                 %Penalty for speed change
%                 f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:vslpenalty-1))+1:...
%                     network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:vslpenalty))) = +we2;
% 
%             else
%                 f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:qDS_n-1))+j) =...
%                     f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:qDS_n-1))+j) - we*T;
%                 
%                 f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:qUS_n-1))+j) =...
%                     f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:qUS_n-1))+j) + we*T;    
%                 
%                 if(any(link==RM_link))
% 
%                     %Penalty for speed change
%                     f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:RMpenalty-1))+1:...
%                         network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:RMpenalty))) = +we2;
% 
%                 end
%             
%             end
%             
%         end
%         
%         if(any(InputID==link))
%             if(any(link==VSL_link)) %penalty for the demand that was not allowed to enter
% 
%                 f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:qUS-1))+n) =...
%                     f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:qUS-1))+n) -we*T;
%             else
% 
%                 f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:qUS_n-1))+n) =...
%                     f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:qUS_n-1))+n) -we*T;        
%             end
%         end
%        
%         
% %         if(any(link==VSL_link))
% %             f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:qDS-1))+n)=...
% %                 f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:qDS-1)+n)) + we2;
% %         else
% %             f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:qDS_n-1))+n)=...
% %                 f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:qDS_n-1)+n)) + we2;
% %         end
% %         if(any(link==VSL_link))
% %             f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:TCD-1))+(n-1)*8+4) = +we;
% %         else
% %             f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:qDS_n-1))+n) =...
% %                 f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:qDS_n-1))+n) + 1/network_links.(linkStr).v;
% %         end
%     end
% end

%% TCD
% Total Jam Density
f = zeros(network_links.totalLength,1);
we_add = 0.0;
we_mult = 0.07;
we = zeros(length(InputID),num_steps);
% we = we_base*ones(length(InputID),num_steps);
we2 = 10;
for link = 1 : num_links
    linkStr = sprintf('link_%d',link);
    for n = num_steps : num_steps
        if(any(link==VSL_link))
            f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:new_Obj-1))+n) = network_links.(linkStr).L;
        elseif(any(link==RM_link))
            f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:new_Obj_n-1))+n) = network_links.(linkStr).L;
        else
            f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:iRM-1))+n) = network_links.(linkStr).L;
        end

    end
    
    %Penalty for VSL and RM control signal variations
    if(any(link==VSL_link))
        %Penalty for speed change
        f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:vslpenalty-1))+1:...
            network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:vslpenalty))) = +we2;            
    elseif(any(link==RM_link))
        %Penalty for ramp conntrol change
        f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:RMpenalty-1))+1:...
            network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:RMpenalty))) = +we2;            
    end 
    
    
    if(any(InputID==link)) %penalty for the demand that was not allowed to enter
        input = find(InputID == link);

%         we(input,1:num_steps) = we_add + we_mult*(Qin_all(input,Simulation_tracker)-Qin_original(input,Simulation_tracker))...
%                           /(Qin_original(input,Simulation_tracker)*num_steps);
        we(input,1:num_steps) = we_add + we_mult*Qin_all(input,Simulation_tracker);
                      
        for n = 1 : num_steps
            
            if(any(link==VSL_link)) 
                f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:qUS-1))+n) =...
                    f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:qUS-1))+n) -we(input,n)*T;
            else
                f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:qUS_n-1))+n) =...
                    f(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:qUS_n-1))+n) -we(input,n)*T;        
            end
        end
    end
end

%% Run Solver

ctype = ch;

% opt = Cplex();
gap = 0.001; %0.05 = 5percent
% options = cplexoptimset('mip.strategy.probe',3,'preprocessing.reduce',3,'emphasis.mip',4,'mip.tolerances.mipgap',gap,...
%         'mip.tolerances.absmipgap',gap,'mip.display',4,'mip.tolerances.integrality',1e-08,'display','on');
% options = cplexoptimset('randomseed',100,'mip.strategy.search',1,'mip.strategy.probe',3,'preprocessing.reduce',3,'emphasis.mip',4,'mip.tolerances.mipgap',gap,...
%         'mip.display',4,'mip.tolerances.integrality',1e-08,'display','on');
% options = cplexoptimset('mip.strategy.probe',3,'preprocessing.reduce',3,'emphasis.mip',4,'mip.tolerances.mipgap',gap,...
%         'mip.display',4,'mip.tolerances.integrality',1e-08,'display','on');
toc

tic
[x, fval, exitflag, output] = ...
    cplexmilp(f, Aineq, bineq, Aeq, beq,[ ], [ ], [ ], lb, ub, ctype, [ ]);
   
% Output_times(run_no,solverRun) = toc
Output_times(Simulation_tracker,1) = toc;
Output_times(Simulation_tracker,1)

fprintf ('\nSolution status = %s \n', output.cplexstatusstring);
fprintf ('Solution value = %f \n', fval);

if(~(output.cplexstatusstring == "integer infeasible"))
   
    comparearray = [];
    for link = 1 : num_links
        linkStr = sprintf('link_%d',link);
        n = num_steps;
        if(any(link==VSL_link))
            comparearray(link,1) = x(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:new_Obj-1))+n);
        elseif(any(link==RM_link))
            comparearray(link,1) = x(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:new_Obj_n-1))+n);
        else
            comparearray(link,1) = x(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:iRM-1))+n);
        end

    end

    %% Arrange Final Objective Value
    % % Total Travel Time
    % FinObj = 0;
    % for link = 1 : num_links
    %     linkStr = sprintf('link_%d',link);
    %     FinObj = FinObj - network_links.(linkStr).ID_n(end)*num_steps*T;
    % end
    % 
    % for input = 1 : length(InputID)
    %     FinObj = FinObj + original.NSupply(input, end);
    % end
    % 
    % FinObj = FinObj + fval

    % Total Jam Density
    FinObj = 0;
    % for link = 1 : num_links
    %     linkStr = sprintf('link_%d',link);
    %     FinObj = FinObj - network_links.(linkStr).ID_n(end)*num_steps*T;
    % end
% 
%     for input = 1 : length(InputID)
%         for n = 1 : num_steps
%             FinObj = FinObj + we(input,n)*(original.NSupply(input,n+1)-original.NSupply(input,n));
%         end
%     end
% 
%     FinObj = FinObj + fval

end
