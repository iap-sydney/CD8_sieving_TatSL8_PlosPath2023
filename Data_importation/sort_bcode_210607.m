%Steffen Docken
%14-4-21
%Code to reorder barcodes based on prevelance over time
%(barcodes seen at first time point are ordered by size
%at that time point, then barcodes seen later, then those
%seen in LNMC but not Plasma, and then those seen only in PBMC

%Code edited on 22-4-21 to remove sequences with "Unique" barcodes or SL8 misreads
%from data

function [Data_Plasma_new, Data_LNMC_new, Data_PBMC_new] = ...
    sort_bcode_210607(Data_Plasma_old, Data_LNMC_old, Data_PBMC_old)

Data_Plasma_new = Data_Plasma_old;
Data_LNMC_new = Data_LNMC_old;
Data_PBMC_new = Data_PBMC_old;

last_bcode_ind = 1;%will track the index + 1 of the lowest detected barcode
% at the last time point



%reording barcodes
for ii = 1:3
    switch ii
        case 1 %order based on Plasma
            Primary_Data = Data_Plasma_new;
        case 2 %order based on LNMC
            Primary_Data = Data_LNMC_new;
        case 3 %order based on PBMC
            Primary_Data = Data_PBMC_new;
    end
    
    if isempty(Primary_Data.dpi) %checking that there is actually data
        %for this sample type
        continue
        
    end
    
    num_timepoints = length(Primary_Data.dpi); %number of time points
    
    
    for kk = 1:num_timepoints
        
        switch ii %getting barcode count above LOD data for ordering
            case 1 %order based on Plasma
                barcode_data_inform_order = ...
                    max(Data_Plasma_new.B.count(:,1,kk) - Data_Plasma_new.B.LOD(:,1,kk), 0);
            case 2 %order based on LNMC
                barcode_data_inform_order = ...
                    max(Data_LNMC_new.B.count(:,1,kk) - Data_LNMC_new.B.LOD(:,1,kk),0);
            case 3 %order based on PBMC
                barcode_data_inform_order = ...
                    max(Data_PBMC_new.B.count(:,1,kk) - Data_PBMC_new.B.LOD(:,1,kk),0);
        end
        
        %barcode ranking
        [bcode_tots_it, I_bcode_it] = sort(barcode_data_inform_order(last_bcode_ind:end), 'descend');
        %sorting barcodes seen at this time point and getting the indices
        %keyboard
        
        
        for ll = 1:3 %looping over Plamsa, LNMC, and PBMC to be reordered
            
            switch ll
                case 1 %reorder Plasma
                    reorder_Data = Data_Plasma_new;
                case 2 %reorder LNMC
                    reorder_Data = Data_LNMC_new;
                case 3 %reorder PBMC
                    reorder_Data = Data_PBMC_new;
            end
            
            if isempty(reorder_Data.dpi) %checking that there is actually data
                %for this sample type
                
                continue
            end
            
            reorder_Data.barcode_ID(last_bcode_ind:end) = ...
                reorder_Data.barcode_ID(last_bcode_ind-1 + I_bcode_it);
            %update order of barcode_ID
            
            num_timepoints_reorder = length(reorder_Data.dpi); %number of time points
            
            for mm = 1:num_timepoints_reorder
                
                for nn = 1:3
                    switch nn
                        case 1 %reordering counts
                            B_data = reorder_Data.B.count;
                            W_data = reorder_Data.W.count;
                            
                        case 2 %reordering LOD
                            B_data = reorder_Data.B.LOD;
                            W_data = reorder_Data.W.LOD;
                            
                        case 3 %reordering LOD_ind
                            B_data = reorder_Data.B.LOD_ind;
                            W_data = reorder_Data.W.LOD_ind;
                    end
                    
                    %Total Barcodes
                    B_data(last_bcode_ind:end, 1, mm) = ...
                        B_data(last_bcode_ind-1 + I_bcode_it, 1, mm);
                
                    % Barcode-SL8 var reordered by barcode
                    W_data(last_bcode_ind:end, :, mm) = ...
                        W_data(last_bcode_ind-1 + I_bcode_it, :, mm);
                    
                    switch nn
                        case 1 %reordering counts
                            reorder_Data.B.count = B_data;
                            reorder_Data.W.count = W_data;
                            
                        case 2 %reordering LOD
                            reorder_Data.B.LOD = B_data;
                            reorder_Data.W.LOD = W_data;
                            
                        case 3 %reordering LOD_ind
                            reorder_Data.B.LOD_ind = B_data;
                            reorder_Data.W.LOD_ind = W_data;
                    end
                    
                end
                
            end
            
            
            switch ll
                case 1 %reorder Plasma
                    Data_Plasma_new = reorder_Data;
                case 2 %reorder LNMC
                    Data_LNMC_new = reorder_Data;
                case 3 %reorder PBMC
                    Data_PBMC_new = reorder_Data;
            end
            
        end
        
        
        if (sum(bcode_tots_it) > 0)
            last_bcode_ind = find(bcode_tots_it, 1, 'last') + last_bcode_ind;
        end %updating last_bcode_ind if anymore barcodes were detected,
        %otherwise leaving it what it was before.
            
        
    end
    
end


