% This file defines the Junction structs
% Each junction struct contains the information about the input and output
% links, and split ratios.
% Junctions can be of four types: 1-1, 2-1, 1-2, 2-2

network_junc = struct;

%% One-to-one
junc =14;
juncStr = sprintf('junc_%d',junc);

network_junc.(juncStr).type = 1; % one-one
network_junc.(juncStr).inlabel = 14; % input links
network_junc.(juncStr).outlabel = 15; %output links
network_junc.(juncStr).Alpha = 1; %allocation to merge links (sum = 1)
network_junc.(juncStr).Alpha_supply = 1; 

%% merge 2-1

junc = 2;

juncStr = sprintf('junc_%d',junc);

network_junc.(juncStr).type = 2; % merge
network_junc.(juncStr).inlabel = [2 20]; % input links
network_junc.(juncStr).outlabel = 3; %output links
network_junc.(juncStr).Alpha = [ 1  1  0;                              
                             0   0  0];
network_junc.(juncStr).Alpha_supply = [ 0.9 0.1  0;                              
                                      0   0   0]; 

junc = 4;

juncStr = sprintf('junc_%d',junc);

network_junc.(juncStr).type = 2; % merge
network_junc.(juncStr).inlabel = [4 22]; % input links
network_junc.(juncStr).outlabel = 5; %output links
network_junc.(juncStr).Alpha = [ 1  1  0;                              
                             0   0  0];
network_junc.(juncStr).Alpha_supply = [ 0.9 0.1  0;                              
                                      0   0   0]; 
                                  
junc = 6;

juncStr = sprintf('junc_%d',junc);

network_junc.(juncStr).type = 2; % merge
network_junc.(juncStr).inlabel = [6 24]; % input links
network_junc.(juncStr).outlabel = 7; %output links
network_junc.(juncStr).Alpha = [ 1  1  0;                              
                             0   0  0];
network_junc.(juncStr).Alpha_supply = [ 0.9 0.1  0;                              
                                      0   0   0]; 

junc = 8;

juncStr = sprintf('junc_%d',junc);

network_junc.(juncStr).type = 2; % merge
network_junc.(juncStr).inlabel = [8 26]; % input links
network_junc.(juncStr).outlabel = 9; %output links
network_junc.(juncStr).Alpha = [ 1  1  0;                              
                             0   0  0];
network_junc.(juncStr).Alpha_supply = [ 0.9 0.1  0;                              
                                      0   0   0];

junc = 10;

juncStr = sprintf('junc_%d',junc);

network_junc.(juncStr).type = 2; % merge
network_junc.(juncStr).inlabel = [10 28]; % input links
network_junc.(juncStr).outlabel = 11; %output links
network_junc.(juncStr).Alpha = [ 1  1  0;                              
                             0   0  0];
network_junc.(juncStr).Alpha_supply = [ 0.9 0.1  0;                              
                                      0   0   0]; 
                                  
junc = 11;

juncStr = sprintf('junc_%d',junc);

network_junc.(juncStr).type = 2; % merge
network_junc.(juncStr).inlabel = [11 29]; % input links
network_junc.(juncStr).outlabel = 12; %output links
network_junc.(juncStr).Alpha = [ 1  1  0;                              
                             0   0  0];
network_junc.(juncStr).Alpha_supply = [ 0.9 0.1  0;                              
                                      0   0   0]; 

junc = 13;

juncStr = sprintf('junc_%d',junc);

network_junc.(juncStr).type = 2; % merge
network_junc.(juncStr).inlabel = [13 31]; % input links
network_junc.(juncStr).outlabel = 14; %output links
network_junc.(juncStr).Alpha = [ 1  1  0;                              
                             0   0  0];
network_junc.(juncStr).Alpha_supply = [ 0.9 0.1  0;                              
                                      0   0   0];
                                  
junc = 16;

juncStr = sprintf('junc_%d',junc);

network_junc.(juncStr).type = 2; % merge
network_junc.(juncStr).inlabel = [16 33]; % input links
network_junc.(juncStr).outlabel = 17; %output links
network_junc.(juncStr).Alpha = [ 1  1  0;                              
                             0   0  0];
network_junc.(juncStr).Alpha_supply = [ 0.8 0.2  0;                              
                                      0   0   0];

%% diverge 1-2

junc = 1;

juncStr = sprintf('junc_%d',junc);

network_junc.(juncStr).type = 3; % diverge
network_junc.(juncStr).inlabel = 1; % input links
network_junc.(juncStr).outlabel = [2 19]; %output links
network_junc.(juncStr).Alpha = [ 0.9    0; 
                              0.1    0;
                               0     0];

network_junc.(juncStr).Alpha_supply = [ 1   0;
                                     1   0;
                                     0   0];     
                                
                                 
junc = 3;

juncStr = sprintf('junc_%d',junc);

network_junc.(juncStr).type = 3; % diverge
network_junc.(juncStr).inlabel = 3; % input links
network_junc.(juncStr).outlabel = [4 21]; %output links
network_junc.(juncStr).Alpha = [ 0.9    0; 
                              0.1    0;
                               0     0];

network_junc.(juncStr).Alpha_supply = [ 1   0;
                                     1   0;
                                     0   0]; 
                                 
junc = 5;

juncStr = sprintf('junc_%d',junc);

network_junc.(juncStr).type = 3; % diverge
network_junc.(juncStr).inlabel = 5; % input links
network_junc.(juncStr).outlabel = [6 23]; %output links
network_junc.(juncStr).Alpha = [ 0.9    0; 
                              0.1    0;
                               0     0];

network_junc.(juncStr).Alpha_supply = [ 1   0;
                                     1   0;
                                     0   0];     
                                
                                 
junc = 7;

juncStr = sprintf('junc_%d',junc);

network_junc.(juncStr).type = 3; % diverge
network_junc.(juncStr).inlabel = 7; % input links
network_junc.(juncStr).outlabel = [8 25]; %output links
network_junc.(juncStr).Alpha = [ 0.8    0; 
                              0.2    0;
                               0     0];

network_junc.(juncStr).Alpha_supply = [ 1   0;
                                     1   0;
                                     0   0];     
                                 
                                 
junc = 9;

juncStr = sprintf('junc_%d',junc);

network_junc.(juncStr).type = 3; % diverge
network_junc.(juncStr).inlabel = 9; % input links
network_junc.(juncStr).outlabel = [10 27]; %output links
network_junc.(juncStr).Alpha = [ 0.9    0; 
                              0.1    0;
                               0     0];

network_junc.(juncStr).Alpha_supply = [ 1   0;
                                     1   0;
                                     0   0]; 
                                 
junc = 12;

juncStr = sprintf('junc_%d',junc);

network_junc.(juncStr).type = 3; % diverge
network_junc.(juncStr).inlabel = 12; % input links
network_junc.(juncStr).outlabel = [13 30]; %output links
network_junc.(juncStr).Alpha = [ 0.9    0; 
                              0.1    0;
                               0     0];

network_junc.(juncStr).Alpha_supply = [ 1   0;
                                     1   0;
                                     0   0];     
                                
                                 
junc = 15;

juncStr = sprintf('junc_%d',junc);

network_junc.(juncStr).type = 3; % diverge
network_junc.(juncStr).inlabel = 15; % input links
network_junc.(juncStr).outlabel = [16 32]; %output links
network_junc.(juncStr).Alpha = [ 0.8    0; 
                              0.2    0;
                               0     0];

network_junc.(juncStr).Alpha_supply = [ 1   0;
                                     1   0;
                                     0   0];  
                                 
junc = 17;

juncStr = sprintf('junc_%d',junc);

network_junc.(juncStr).type = 3; % diverge
network_junc.(juncStr).inlabel = 17; % input links
network_junc.(juncStr).outlabel = [18 34]; %output links
network_junc.(juncStr).Alpha = [ 0.9    0; 
                              0.1    0;
                               0     0];

network_junc.(juncStr).Alpha_supply = [ 1   0;
                                     1   0;
                                     0   0];  
    