% This script defines all the necessary parameters for optimization

%% Control framework
% All_v = []; %Execute once and comment
% All_IDk = []; %Execute once and comment
%Simulation Configuration
num_steps = 10; %N_p = Prediction Horizon time-steps
T = 60; %Duration of each time-step; minimum value = (chi - xi)/(minimum v)

N_c = 1;
N_p = 10; %not used
%Network Compenents
num_links = 34;
VSL_link = [15 16];
% VSL_link = [];
RM_link = [29 31];
% RM_link=[]; 
num_juncs = 17; %1 one-one, merge, 0 diverge, 0 two-two: Defined in 'network_struct'

onramp = [20 22 24 26 28 29 31 33];
offramp = [19 21 23 25 27 30 32 34];

%Link Configuration
% xi =0; chi = 1000; %length of link in meters: Can be modified for each link in 'network_struct'
% L = chi - xi;

%Switch Configuration
Sw = [5]; %Speed Limit switch timings; minimum gap between switches = (chi - xi)/(w*T)
num_sec = size(Sw,2)+1; %number of sections separated by switches

%% Traffic Flow Diagram Parameters (Common to entire network - for now)
vel = [19 25 33]; %possible speed limit/free-flow speed settings
w = -5; %congestion velocity
km = 0.5; %maximum density
kcrit = (km.*w)./(w-vel); %critical density
Q = kcrit.*vel; %maximum/critical flow

v = 33; %for no VSL links
% kc = (km.*w)./(w-v); %for no VSL links

%Initial Density Conditions
% ID_x = [0 L]; %if L = 500

%Demand and Supply at network boundaries
% Network Inputs
InputID = [1 20 22 24 26 28 29 31 33];
% Network Outputs
OutputID = [18 19 21 23 25 27 30 32 34];

NSupply = zeros(length(InputID), num_steps+1);
NDemand = zeros(length(OutputID), num_steps+1);


%Limiting values for some variables
epsilon = 1e-5;
Mk = 100;
mk = -100;
Mint = 10000;
mint = -10000;
vmax = max(vel);

%% Variable counting for Projection and calculation
%array to store the number of different variables
sizes = [num_steps  num_steps  num_steps  num_steps*3  num_sec*3 num_sec-1 num_steps]; %does not include Dem and Sup or their binary variables
size_row = sum(sizes)+1; %last column to store constant value in constraints

sizes_noVSL = [num_steps num_steps num_steps];
size_row_noVSL = sum(sizes_noVSL)+1;

sizes_RM = [num_steps num_steps num_steps num_steps-1 num_steps];
size_row_RM = sum(sizes_noVSL)+1;
    
%% Variable indices
%position of variables in that array for VSL links
qUS = 1; %upstream flows
kUS = 2; %upstream density
qDS = 3; %downstream flows
ka = 4; %auxilliary continuoius variables for upstream density
bin = 5; %binary variables for velocity selection
%     TCD = 6; % to store flow/v for Total Congestion Delay
Dem = 6; %demand variables
Sup = 7; %supply variables
nbds = 8; %num bin downstream
nbus = 9; %num bin upstream
vslpenalty = 10; %stores absolute difference between speeds in consecutive sections
new_Obj = 11;

%position of variables in that array for non-VSL links
qUS_n = 1; %upstream flows
qDS_n = 2; %downstream flows
Dem_n = 3; %demand variables
Sup_n = 4; %supply variables
nbds_n = 5; %num bin downstream
nbus_n = 6; %num bin upstream
iRM = 7; %in RM controlled links
RMpenalty = 8;
new_Obj_n = 9;

%% Set In-flow and Out-flow for entire simulation period

%% Solver Runs
rng(8);
num_runs = 1;
Initial_density = zeros(4, num_runs);
runs_per_condition = 1;

Simulation_tracker = 0;

Simulation_tracker = Simulation_tracker + 1;
Demand_Supply_All
for run_no = 1 : num_runs
    Junctions
    network_struct
    solver
    Result_interpretation
end