classdef setJConstraints
    properties
        num_juncs
        network_links
        network_junc
        
        VSL_link
    end
    
    methods
        function si = setJConstraints(num_juncs, network_links, network_junc, VSL_link)
            si.num_juncs = num_juncs;
            si.network_links = network_links;
            si.network_junc = network_junc;
            
            si.VSL_link = VSL_link;
        end
        
        
        function [list] = setMatrixDS(si, num_steps, T)    
            
            list = zeros(0,si.network_links.totalLength+1);
%             num_junc = sum(structfun(@numel,si.network_junc));
%             on_count = 0;
            
            
            %Binary variables useful counter as a reference
            countB=0;
            Cmax = 500000; 
            
            % nb: number of binary variables per point THIS CAN BE REDUCED
            nb = ceil(log2(10));
            
            % comb_mat: combination matrix for possible binary combinations
            comb_mat = zeros(2^nb,nb);
%             epsilon = 0.00001;
            epsilon = 0;
            for i = 1:2^nb
                comb_mat(i,:) = double(dec2bin(i-1,nb)=='1');
            end
            
            
            for junc = 1 : si.num_juncs
                
                
                juncStr = sprintf('junc_%d',junc);
                num_in = length(si.network_junc.(juncStr).inlabel);
                num_out = length(si.network_junc.(juncStr).outlabel);
                
%                 on_count = on_count + si.network_junc.(juncStr).on(1,1); %0 as there are no ramps
                                                                               
                     for in = 1 : num_in
                        
                        for ti = 1 : num_steps
                            
                            tempM = zeros(0,0);
                            tempArray = zeros(1,si.network_links.totalLength+1);
                            
                            linkStr = sprintf('link_%d',si.network_junc.(juncStr).inlabel(1,in));
                            
                            link = si.network_junc.(juncStr).inlabel(1,in);
                            if(any(link == si.VSL_link))
                                qDS = 3;
                                Dem = 6;
                            else
                                qDS = 2;
                                Dem = 3;
                            end
                            
                            tmp_offset = si.network_links.(linkStr).offset+...
                                sum( si.network_links.(linkStr).sizes(1:qDS-1) );
                            
                            tempArray(1, tmp_offset + ti)...
                                = 1;
                             tempArray(1, si.network_links.totalLength+1) = epsilon;
                             
                            tmp_offset = si.network_links.(linkStr).offset+...
                                sum( si.network_links.(linkStr).sizes(1:Dem-1) );
                            
                            arraym = zeros(1,si.network_links.totalLength+1);
                            
                            arraym(1, tmp_offset + 1)= -1/T;
                            %arraym(1, tmp_offset + ti)= -1/T;
                            arraym(1, tmp_offset + ti+1)= 1/T; 
                            
                            if(ti > 1)
                                tmp_offset = si.network_links.(linkStr).offset+...
                                sum( si.network_links.(linkStr).sizes(1:qDS-1) );
                                arraym(1, tmp_offset + 1: tmp_offset + ti-1)= -1;
                            end
                           
                            array2 = tempArray - arraym;
                            if(~isempty(array2))
                                rows = size(list,1);
                                list(rows+1,:) = array2;
                                rows = size(tempM,1);
                                tempM(rows+1,:) = arraym;
                            end
                           
                           for out = 1 : num_out
                               
                             if(si.network_junc.(juncStr).Alpha(out,in) > 0)
                           
                               linkStr = sprintf('link_%d',si.network_junc.(juncStr).outlabel(1,out));
                               
                               link = si.network_junc.(juncStr).outlabel(1,out);
                               if(any(link == si.VSL_link))
                                    qUS = 1;
                                    Sup = 7;
                                else
                                    qUS = 1;
                                    Sup = 4;
                               end
                               
                               tmp_offset = si.network_links.(linkStr).offset +...
                               sum( si.network_links.(linkStr).sizes(1:Sup-1) );
                         
                               arraym = zeros(1,si.network_links.totalLength+1);
                               n2flow = 1/T;
                             
                              arraym(1, tmp_offset +1)= -si.network_junc.(juncStr).Alpha_supply(out,in)*n2flow/si.network_junc.(juncStr).Alpha(out,in);  
                              arraym(1, tmp_offset + ti+1)= si.network_junc.(juncStr).Alpha_supply(out,in)*n2flow/si.network_junc.(juncStr).Alpha(out,in);  
                              
                              if(ti>1)
                                
                                tmp_offset = si.network_links.(linkStr).offset+...
                                sum( si.network_links.(linkStr).sizes(1:qUS-1) );  
                                 
                                arraym(1, tmp_offset +1:tmp_offset + ti-1)= -si.network_junc.(juncStr).Alpha_supply(out,in)*1/si.network_junc.(juncStr).Alpha(out,in);  
                                  
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
                           
                           rowsM = size(tempM,1);
                           
                      for i = 1:rowsM
                        
                        array = tempArray;
                      
                        %Decode the binary combinations
                        
                        for counter = 1 : si.network_junc.(juncStr).num_juncBinaries(in)
                            if comb_mat(i,nb-si.network_junc.(juncStr).num_juncBinaries(in)+counter) == 1
                                tmp_offset = si.network_junc.offset;
                               
                                array(1,tmp_offset + countB+ counter) = -Cmax;
                                
                            elseif comb_mat(i,nb-si.network_junc.(juncStr).num_juncBinaries(in)+counter) == 0
                                tmp_offset = si.network_junc.offset;
                                array(1,tmp_offset + countB+ counter) = Cmax;
                                
                            end
                        end
                        array(1,si.network_links.totalLength+1) = -Cmax*(sum(comb_mat(i,:)))-epsilon; %RHS (set negative to be on same side)
                        
                        % Add constraint to the MILP matrix
                        
                        
                        array2 = + tempM(i,:) - array;
                        rows = size(list,1);
                        list(rows+1,:) = array2;
                     end
                    
                     %Define the last constraint to bound the binary combination
                    
                    array = zeros(1,si.network_links.totalLength+1);
                    
                    for counter=1:si.network_junc.(juncStr).num_juncBinaries(in)
                        tmp_offset = si.network_junc.offset;
                        array(1, tmp_offset + countB + counter) = 2^(si.network_junc.(juncStr).num_juncBinaries(in)-counter);
                        
                    end
                    
                    array(1,si.network_links.totalLength+1) = (rowsM-1); %RHS (maximum possible value of the binary comb)
                    rows = size(list,1);
                    list(rows+1,:) = array;
                    
                    countB = countB + si.network_junc.(juncStr).num_juncBinaries(in);
                           
                        end
                        
                     end             
                                                          
            end
            
        end
            
        
        function [list]= setMatrixJunc(si, num_steps)
%             qUS = 1; %upstream flows
%             qDS = 3; %downstream flows
                        
            list = zeros(0,si.network_links.totalLength);
%             num_junc = sum(structfun(@numel,si.network_junc));
%             on_count = 0;
            for junc = 1 : si.num_juncs
                juncStr = sprintf('junc_%d',junc);
                num_in = length(si.network_junc.(juncStr).inlabel);
                num_out = length(si.network_junc.(juncStr).outlabel);
                
%                 on_count = on_count + si.network_junc.(juncStr).on(1,1);
                
                for out = 1 : num_out
                    
                    for ti = 1 : num_steps
                        tempArray = zeros(1,si.network_links.totalLength);
                        
                        % set Incoming road
                        for in = 1 : num_in
                            linkStr = sprintf('link_%d',si.network_junc.(juncStr).inlabel(1,in));
                            
                            link = si.network_junc.(juncStr).inlabel(1,in);
                           if(any(link == si.VSL_link))
                                qDS = 3;
                            else
                                qDS = 2;
                           end
                                
                            tmp_offset = si.network_links.(linkStr).offset+...
                                sum( si.network_links.(linkStr).sizes(1:qDS-1));
                            
                            tempArray(1, tmp_offset + ti)...
                                = si.network_junc.(juncStr).Alpha(out,in);
                            
                        end
                        
                        % set outgoing road
                        linkStr = sprintf('link_%d',si.network_junc.(juncStr).outlabel(1,out));
                        
                        link = si.network_junc.(juncStr).outlabel(1,out);
                           if(any(link == si.VSL_link))
                                qUS = 1;
                            else
                                qUS = 1;
                           end
                        
                        tmp_offset = si.network_links.(linkStr).offset +...
                            sum( si.network_links.(linkStr).sizes(1,1:qUS-1) );
                        
                        tempArray(1, tmp_offset + ti) = -1;
                        
                        rows = size(list,1);
                        list(rows+1,:) = tempArray;
                        
                    end
                end
            end
            
        end
        
    end
    
end