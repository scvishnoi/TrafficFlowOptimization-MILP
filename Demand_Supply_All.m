total_time = 7*num_steps;
if(Simulation_tracker == 1)
    Qin_all=[];
    % InputID =
    %      1    20    22    24    26    28    29    31    33
    Qin_all(1,:) = [1.8*ones(1,num_steps) 1.9*ones(1,num_steps) 2.0*ones(1,num_steps) 2.1*ones(1,num_steps) 2.0*ones(1,num_steps) 1.9*ones(1,num_steps) 1.9*ones(1,num_steps)];        
    Qin_all(2,:) = [0.25*ones(1,num_steps) 0.28*ones(1,num_steps) 0.28*ones(1,num_steps) 0.3*ones(1,num_steps) 0.29*ones(1,num_steps) 0.25*ones(1,num_steps) 0.25*ones(1,num_steps)];
    Qin_all(3,:) = [0.25*ones(1,num_steps) 0.28*ones(1,num_steps) 0.28*ones(1,num_steps) 0.3*ones(1,num_steps) 0.29*ones(1,num_steps) 0.25*ones(1,num_steps) 0.25*ones(1,num_steps)];
    Qin_all(4,:) = [0.28*ones(1,num_steps) 0.3*ones(1,num_steps) 0.28*ones(1,num_steps) 0.3*ones(1,num_steps) 0.28*ones(1,num_steps) 0.26*ones(1,num_steps) 0.26*ones(1,num_steps)];
    Qin_all(5,:) = [0.28*ones(1,num_steps) 0.3*ones(1,num_steps) 0.28*ones(1,num_steps) 0.3*ones(1,num_steps) 0.28*ones(1,num_steps) 0.26*ones(1,num_steps) 0.26*ones(1,num_steps)];
    Qin_all(6,:) = [0.2*ones(1,num_steps) 0.24*ones(1,num_steps) 0.25*ones(1,num_steps) 0.25*ones(1,num_steps) 0.24*ones(1,num_steps) 0.22*ones(1,num_steps) 0.22*ones(1,num_steps)];
    Qin_all(7,:) = [0.25*ones(1,num_steps) 0.26*ones(1,num_steps) 0.25*ones(1,num_steps) 0.28*ones(1,num_steps) 0.26*ones(1,num_steps) 0.28*ones(1,num_steps) 0.28*ones(1,num_steps)];
    Qin_all(8,:) = [0.25*ones(1,num_steps) 0.26*ones(1,num_steps) 0.25*ones(1,num_steps) 0.28*ones(1,num_steps) 0.26*ones(1,num_steps) 0.28*ones(1,num_steps) 0.28*ones(1,num_steps)];
    Qin_all(9,:) = [0.3*ones(1,num_steps) 0.28*ones(1,num_steps) 0.26*ones(1,num_steps) 0.28*ones(1,num_steps) 0.3*ones(1,num_steps) 0.3*ones(1,num_steps) 0.3*ones(1,num_steps)];            
    Qin_original = Qin_all;
end

Qout_all=[];
for output = 1 : length(OutputID)
    if(any(OutputID(output)==offramp))
        if(OutputID(output)==32)
            Qout =2*ones(1,total_time);            
        else
            Qout = 2*ones(1, total_time); %no restriction - always above the maximum flow            
        end
    else
        Qout = [2.0*ones(1,num_steps) 2.2*ones(1,num_steps) 1.9*ones(1,num_steps) 2.0*ones(1,num_steps) 1.9*ones(1,num_steps) 2.0*ones(1,num_steps) 2.0*ones(1,num_steps)];           
    end
    
    Qout_all(output, :) = Qout;
end