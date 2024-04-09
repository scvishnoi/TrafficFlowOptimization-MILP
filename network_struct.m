% This script:
% - Creates link structs
% - Assigns parameters to link structs
% - Creates junction structs
% - Assigns parameters to junction structs inclusing link connections


%% Creating network structs
%network structs
network_links = struct;
network_config
%assigning parameters to link structs : Storing the traffic flow
%parameters, information about link length and the switches

for junc = 1 : num_juncs
    juncStr = sprintf('junc_%d',junc);
    for i = 1 : length(network_junc.(juncStr).inlabel)
        linkStr = sprintf('link_%d',network_junc.(juncStr).inlabel(i));
        network_links.(linkStr).junc_indicator(1) = 1;
    end
    
    for i = 1 : length(network_junc.(juncStr).outlabel)
        linkStr = sprintf('link_%d',network_junc.(juncStr).outlabel(i));
        network_links.(linkStr).junc_indicator(2) = 1;
    end
end

%%
%Calculating number of junction binaries
total_bins = 0;
num_juncBinaries=0;
for junc = 1 : num_juncs
    
    juncStr = sprintf('junc_%d',junc);
    num_in = length(network_junc.(juncStr).inlabel);
    num_out = length(network_junc.(juncStr).outlabel);
    num_juncBinaries = zeros(1,num_in);
    
    for in=1:num_in
       
        opts = 1;
        
        for out=1:num_out
           
            if(network_junc.(juncStr).Alpha(out,in)>0)
                opts = opts+1;
            end
            
        end
        
    num_juncBinaries(in)= ceil(log2(opts));    
    end
    
  network_junc.(juncStr).num_juncBinaries = num_juncBinaries;
  total_bins = total_bins + sum(network_junc.(juncStr).num_juncBinaries);
end


%% Set Constraints
Aeq = zeros(0,0);
Aineq = zeros(0,0);

network_links.totalLength = 0;


for link = 1 : num_links
    
    linkStr = sprintf('link_%d',link);
    if(any(link==VSL_link)) %if link is a VSL controlled link
        si = setConstraints(network_links.(linkStr).L, network_links.(linkStr).chi,...
            network_links.(linkStr).xi, network_links.(linkStr).Sw, network_links.(linkStr).vel,...
            network_links.(linkStr).w, network_links.(linkStr).km, network_links.(linkStr).kcrit,...
            network_links.(linkStr).Q, network_links.(linkStr).sizes, network_links.(linkStr).size_row, network_junc, num_steps, T,...
            network_links.(linkStr).junc_indicator, network_links.(linkStr).ID_x, network_links.(linkStr).ID_k, network_links.(linkStr).ID_n);

        network_links.(linkStr).sizes = si.sizes;
        network_links.(linkStr).size_row = si.size_row;
        % Created for each link
        [network_links.(linkStr).arrayeq, network_links.(linkStr).arrayineq] = si.setAllConstraints(num_steps,T);
        [network_links.(linkStr).DS_matrix] = si.setDemandSupply(num_steps,T); %demand supply variables constraints

        %To incorporate network demands and supplies
        if(network_links.(linkStr).junc_indicator(2) == 0)
            flag = 1;
            [network_links.(linkStr).input] = si.setBoundaryDS(NSupply(InputID == link,:), num_steps, T, flag);
            [network_links.(linkStr).input2] = si.setMatrixDS(num_steps, T, flag);
        end
        if(network_links.(linkStr).junc_indicator(1) == 0)
            flag = 2;
            [network_links.(linkStr).output] = si.setBoundaryDS(NDemand(OutputID == link,:), num_steps, T, flag);
            [network_links.(linkStr).output2] = si.setMatrixDS(num_steps, T, flag);
        end
    else %if the link is not VSL controlled - will also include ramps (some of which are ramp controlled)
        si = setConstraints_noVSL(network_links.(linkStr).L, network_links.(linkStr).chi,...
            network_links.(linkStr).xi, network_links.(linkStr).v,...
            network_links.(linkStr).w, network_links.(linkStr).km, network_links.(linkStr).kc,...
            network_links.(linkStr).sizes, network_links.(linkStr).size_row, network_junc, num_steps, T,...
            network_links.(linkStr).junc_indicator, network_links.(linkStr).ID_x, network_links.(linkStr).ID_k, network_links.(linkStr).ID_n,network_links.(linkStr).RM);

        network_links.(linkStr).sizes = si.sizes;
        network_links.(linkStr).size_row = si.size_row;
        % Created for each link
        [network_links.(linkStr).arrayineq] = si.setAllConstraints(num_steps,T);
        [network_links.(linkStr).DS_matrix] = si.setDemandSupply(num_steps,T); %demand supply variables constraints


        %To incorporate network demands and supplies
        if(network_links.(linkStr).junc_indicator(2) == 0)
            flag = 1;
            [network_links.(linkStr).input] = si.setBoundaryDS(NSupply(InputID == link,:), num_steps, T, flag);
            [network_links.(linkStr).input2] = si.setMatrixDS(num_steps, T, flag);
        end
        if(network_links.(linkStr).junc_indicator(1) == 0)
            flag = 2;
            [network_links.(linkStr).output] = si.setBoundaryDS(NDemand(OutputID == link,:), num_steps, T, flag);
            [network_links.(linkStr).output2] = si.setMatrixDS(num_steps, T, flag);
        end

    end
    
    if (link == 1)    
        network_links.(linkStr).offset = 0;
    else
        pre_link_namestr = sprintf('link_%d',link-1);
        network_links.(linkStr).offset = network_links.(pre_link_namestr).offset + network_links.(pre_link_namestr).count;
    end
    network_links.(linkStr).count = sum(si.sizes);
    network_links.totalLength = network_links.totalLength + network_links.(linkStr).count; 
end

network_junc.offset = network_links.totalLength;

for junc = 1 : num_juncs
    juncStr = sprintf('junc_%d',junc);
    if(junc == 1)
        network_junc.(juncStr).offset = network_junc.offset;
    else
        pre_junc_namestr = sprintf('junc_%d',junc-1);
        network_junc.(juncStr).offset = network_junc.(pre_junc_namestr).offset + network_junc.(pre_junc_namestr).count;
    end
    network_junc.(juncStr).count = sum(network_junc.(juncStr).num_juncBinaries)*num_steps; 
end
network_links.totalLength = network_links.totalLength + total_bins*num_steps; %including binary variables

si = setJConstraints(num_juncs, network_links, network_junc, VSL_link);

% % Created Once for the network
Junc_ineq = struct;
[Junc_ineq.DS_matrixEq] = si.setMatrixDS(num_steps,T); %Demand Supply Inequality constraints
[Junc_ineq.Junc_matrix] = si.setMatrixJunc(num_steps); % Junction equality constraints