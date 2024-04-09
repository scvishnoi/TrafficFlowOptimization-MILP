vsel = zeros(num_links,size(Sw,2)+1);
vseltemp = zeros(1,size(Sw,2)+1);

RMc = zeros(num_steps, length(RM_link));

for link = 1 : num_links

    linkStr = sprintf('link_%d',link);

    if(any(link == VSL_link))

        for sec = 1 : size(Sw,2)+1
            temp = x(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:bin-1))+(sec-1)*3+1 :network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:bin-1))+(sec-1)*3+3);
            b = find( temp >= 1-0.00001 & temp <= 1+0.00001);

            vseltemp(sec) = network_links.(linkStr).vel(b);

        end
        
        vsel(link,1:(size(Sw,2)+1)) = vseltemp;
    else
        vsel(link,1:(size(Sw,2)+1)) = (network_links.(linkStr).v)*ones(1,size(Sw,2)+1);
    end
    
    
    if(any(link==RM_link))
        RMc(:, find(link==RM_link)) = x(network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:iRM-1)) + 1:...
           network_links.(linkStr).offset + sum(network_links.(linkStr).sizes(1:iRM)));
    end
end

% All_v = [All_v vsel(:,1)];
% All_v = vsel(:,1);

Allvsel = [Allvsel vsel(:,1)];

AllRMc = [AllRMc; RMc(1,:)];

main_mod1
ini_den_split
%this last line (code) stores the new initial densities in newID_x and
%newID_k elements of the sim_links.(linkStr) structure