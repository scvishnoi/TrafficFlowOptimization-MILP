%% Configuration
num_ml = 18;
num_on = 8;
num_off = 8;

%Mainline links
ml = 1 : num_ml;
ml_links = struct;

link_lengths = [0 0.61 0.91 1.38 1.61 1.91 2.42 2.84 3.65 4.35 5.15 5.92 7.09 7.42 8.12 8.81 9.30 9.91 10.13]; %from maps
ml_links.Lengths = 1609.34*(link_lengths(2:end)-link_lengths(1:(end-1)));
km_ml = 0.5; %maximum density for mainline section
w = -5; v = 33;
kcritical = km_ml*w/(w-v);
% ml_links.V = []; %speed limit; 0 if variable
ml_links.ID_k = zeros(1,num_ml-1);

%Initial Densities on all mainline links
ml_links.ID_k([1 2 3 4]) = 0.02; 
ml_links.ID_k([5 6 7 8 9 10]) = kcritical; 
ml_links.ID_k([11 12 13 14]) = 0.01;
ml_links.ID_k([15 16]) = kcritical; 
ml_links.ID_k([17 18])= kcritical;

km_ramps = 0.5/4; %since the ramp will have one lane and the highway had 4 lanes
kcriticalramps = km_ramps*w/(w-v);
%on-ramps
onramp = [20 22 24 26 28 29 31 33];
on_ramps = struct;

on_ramps.Lengths = ones(1,num_on)*1000; %arbitrary 1km length
on_ramps.V = v*ones(1,num_on); 
on_ramps.ID_k = kcriticalramps*ones(1,num_on);

%off-ramps
offramp = [19 21 23 25 27 30 32 34];
off_ramps = struct;

off_ramps.Lengths = ones(1,num_on)*1000;
off_ramps.V = v*ones(1,num_off); 
off_ramps.ID_k = kcriticalramps*ones(1,num_off);

Lengths = [ml_links.Lengths on_ramps.Lengths off_ramps.Lengths];
ID_k = [ml_links.ID_k on_ramps.ID_k off_ramps.ID_k];

%% Storing data in links
for link = 1 : num_links
    
    linkStr = sprintf('link_%d',link);
    
    network_links.(linkStr).L = Lengths(link);
    network_links.(linkStr).chi = Lengths(link);
    network_links.(linkStr).xi = 0;
    network_links.(linkStr).w = w;
    network_links.(linkStr).Sw = Sw;
    if(any(link==ml))
        network_links.(linkStr).km = km_ml;
    else
        network_links.(linkStr).km = km_ramps;
    end
    if(any(link==VSL_link))
        network_links.(linkStr).vel = vel;
        network_links.(linkStr).Sw = Sw;
        network_links.(linkStr).kcrit = ((network_links.(linkStr).km).*network_links.(linkStr).w)./(network_links.(linkStr).w-network_links.(linkStr).vel); %critical density;
        network_links.(linkStr).Q = (network_links.(linkStr).kcrit).*network_links.(linkStr).vel;
        network_links.(linkStr).sizes = sizes;
        network_links.(linkStr).size_row = size_row;
    else
        network_links.(linkStr).v = v;
        network_links.(linkStr).kc = ((network_links.(linkStr).km).*network_links.(linkStr).w)./(network_links.(linkStr).w-network_links.(linkStr).v); %critical density;
        if(any(link==RM_link))
            network_links.(linkStr).sizes = sizes_RM;
            network_links.(linkStr).size_row = size_row_RM;
            network_links.(linkStr).RM = 1;
        else
            network_links.(linkStr).sizes = sizes_noVSL;
            network_links.(linkStr).size_row = size_row_noVSL;
            network_links.(linkStr).RM = 0;
        end
    end
    
   
    if(Simulation_tracker == 1)
        original.(linkStr).ID_x = [0 network_links.(linkStr).L];
        original.(linkStr).ID_k = [round(ID_k(link),5)];
        original.(linkStr).ID_n = [0 -cumsum(((original.(linkStr).ID_x(2:end))-(original.(linkStr).ID_x(1:end-1))).*(original.(linkStr).ID_k))];
 
        network_links.(linkStr).ID_x = [0 network_links.(linkStr).L];
        network_links.(linkStr).ID_k = [round(ID_k(link),5)];
        All_v = []; AllRMc = []; Allvsel = [];
    else
        original.(linkStr).ID_x = sim_links.(linkStr).newID_x;
        original.(linkStr).ID_k = round(sim_links.(linkStr).newID_k,5);
        original.(linkStr).ID_n = [0 -cumsum(((original.(linkStr).ID_x(2:end))-(original.(linkStr).ID_x(1:end-1))).*(original.(linkStr).ID_k))];
    
        network_links.(linkStr).ID_x = sim_links.(linkStr).newID_x;
        network_links.(linkStr).ID_k = round(sim_links.(linkStr).newID_k,5);
    end

    network_links.(linkStr).ID_n = [0 -cumsum(((network_links.(linkStr).ID_x(2:end))-(network_links.(linkStr).ID_x(1:end-1))).*(network_links.(linkStr).ID_k))]; %known Cumulative vehicles values;
    network_links.(linkStr).junc_indicator =[0 0];
  
end

%% Setting Demand and Supply

for input = 1 : length(InputID)
    original.NSupply(input,:) = [0 cumsum(Qin_all(input,Simulation_tracker:Simulation_tracker+num_steps-1)*T)];
    NSupply(input, :) = [0 cumsum(Qin_all(input,Simulation_tracker:Simulation_tracker+num_steps-1)*T)];
end

for output = 1 : length(OutputID)
    linkStr = sprintf('link_%d',OutputID(output));
    original.NDemand(output, :) = original.(linkStr).ID_n(end) + [0 cumsum(Qout_all(output,Simulation_tracker:Simulation_tracker+num_steps-1)*T)];
    NDemand(output, :) = network_links.(linkStr).ID_n(end) + [0 cumsum(Qout_all(output,Simulation_tracker:Simulation_tracker+num_steps-1)*T)];
end
