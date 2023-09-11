%Steffen Docken
%07-06-21
%This code generates arrays of counts for barcode, SL8 variant, and
%barcode-SL8 variant combos

function [Data_SL8var_Plasma, Data_SL8var_LNMC, Data_SL8var_PBMC] = ...
    SL8_variant_count_210607(Data_orig_Plasma, Data_orig_LNMC, ...
    Data_orig_PBMC)

if ((~isempty(Data_orig_Plasma.dpi))&&(~isempty(Data_orig_LNMC.dpi))&&...
        (sum(sum((Data_orig_Plasma.SL8_ID - Data_orig_LNMC.SL8_ID).^2)) ||...
        sum((Data_orig_Plasma.barcode_ID - Data_orig_LNMC.barcode_ID).^2)))
    disp('ERROR: barcode or SL8 info not the same across Plasma and LNMC');
end

if ((~isempty(Data_orig_Plasma.dpi))&&(~isempty(Data_orig_PBMC.dpi))&&...
        (sum(sum((Data_orig_Plasma.SL8_ID - Data_orig_PBMC.SL8_ID).^2)) ||...
        sum((Data_orig_Plasma.barcode_ID - Data_orig_PBMC.barcode_ID).^2)))
    disp('ERROR: barcode or SL8 info not the same across Plasma and PBMC');
end

num_barcode = length(Data_orig_Plasma.barcode_ID);

%collecting SL8 info from original data set, which will be combined into
%counts for each individual variant
SL8_ID_array = Data_orig_Plasma.SL8_ID;

SL8_var = unique(SL8_ID_array(1,:));
num_SL8_var = length(SL8_var);

SL8_var_array = zeros(9, num_SL8_var); %will contain the variant ID (first
%entry) and the code for the amino acid sequence (entries 2:9)

%copying over data from full data set that is unchanged for variant
%specific data set
for ii = 1:3
    switch ii
        case 1 %create new Plasma Data struct
            Data_orig = Data_orig_Plasma;
        case 2 %create new LNMC Data struct
            Data_orig = Data_orig_LNMC;
        case 3 %create new PBMC Data struct
            Data_orig = Data_orig_PBMC;
    end
    
    Data_SL8var.T = Data_orig.T;
    Data_SL8var.S = Data_orig.S;
    Data_SL8var.A = Data_orig.A;
    Data_SL8var.dpi = Data_orig.dpi;
    Data_SL8var.phase = Data_orig.phase;
    Data_SL8var.timing = Data_orig.timing;
    Data_SL8var.sample = Data_orig.sample;
    Data_SL8var.seq = Data_orig.seq;
    Data_SL8var.barcode_ID = Data_orig.barcode_ID;
    Data_SL8var.B = Data_orig.B;
    
    V.count = zeros(1,num_SL8_var, length(Data_orig.dpi));
    V.LOD = zeros(1,num_SL8_var, length(Data_orig.dpi));
    V.LOD_ind = zeros(1,num_SL8_var, length(Data_orig.dpi));
    Data_SL8var.V = V; %structure to hold count, LOD, and above detection
    %indicator for variants
    
    W_var.count = zeros(num_barcode, num_SL8_var, length(Data_orig.dpi));
    W_var.LOD = zeros(num_barcode, num_SL8_var, length(Data_orig.dpi));
    W_var.LOD_ind = zeros(num_barcode, num_SL8_var, length(Data_orig.dpi));
    Data_SL8var.W_var = W_var; %structure to hold count, LOD, and above detection
    %indicator for barcode-variant combos

    Data_SL8var.animal_ID = Data_orig.animal_ID;
    
    switch ii
        case 1
            Data_SL8var_Plasma_init = Data_SL8var;
        case 2
            Data_SL8var_LNMC_init = Data_SL8var;
        case 3
            Data_SL8var_PBMC_init = Data_SL8var;
    end
end


%% summing columns to get totals for each variant
for ii = 1:3
    switch ii
        case 1 %create new Plasma Data struct
            Data_orig = Data_orig_Plasma;
            Data_SL8var = Data_SL8var_Plasma_init;
        case 2 %create new LNMC Data struct
            Data_orig = Data_orig_LNMC;
            Data_SL8var = Data_SL8var_LNMC_init;
        case 3 %create new PBMC Data struct
            Data_orig = Data_orig_PBMC;
            Data_SL8var = Data_SL8var_PBMC_init;
    end
    
    if isempty(Data_orig.dpi)
        continue
    end
    
    for jj = 1:num_SL8_var
        var_ind = find(~(SL8_ID_array(1,:) - SL8_var(jj)).^2); %columns of
        %M_orig and W_orig that pertain to current variant

        if (sum(sum((SL8_ID_array(2:9,var_ind) - ...
                SL8_ID_array(2:9, var_ind(1))).^2)))
            disp(strcat("ERROR: variant not identical aa in animal ", ...
                Data_SL8var.animal_ID));
        end %checking that amino acid sequence is same for all columns registered
        %as this variant

        SL8_var_array(:,jj) = SL8_ID_array(1:9, var_ind(1)); %recording SL8 
        %variant and amino acid sequence
        
        Data_SL8var.V.count(1,jj,:) = sum(Data_orig.M.count(1,var_ind, :), 2);
        Data_SL8var.W_var.count(:, jj, :) = sum(Data_orig.W.count(:, var_ind, :), 2);
        %summing columns for current variant
        
    end
    
    num_timepoint = length(Data_orig.dpi);
    for jj = 1:num_timepoint
        [~, Data_SL8var.V.LOD(1,:,jj), Data_SL8var.W_var.LOD(:,:,jj)] = ...
            bcode_SL8_LOD_210604(Data_SL8var.B.count(:,1,jj), ...
            Data_SL8var.V.count(1,:,jj), Data_SL8var.S(jj), Data_SL8var.A(jj));
    end
    
    Data_SL8var.V.LOD_ind = Data_SL8var.V.count > Data_SL8var.V.LOD;
    Data_SL8var.W_var.LOD_ind = Data_SL8var.W_var.count > Data_SL8var.W_var.LOD;
    
    Data_SL8var.SL8_ID = SL8_var_array; %storing the array with SL8 ID info
    
    switch ii
        case 1 %updating Plasma Data struct
            Data_SL8var_Plasma_init = Data_SL8var;
        case 2 %updating LNMC Data struct
            Data_SL8var_LNMC_init = Data_SL8var;
        case 3 %updating PBMC Data struct
            Data_SL8var_PBMC_init = Data_SL8var;
    end
    
end



%% reordering variants based on number of mutations (0, then 1, etc). First,
%ordering by size and time of first detected. Then this will be subordering
%for when order by number of mutations 
%(variants seen at first time point are ordered by size
%at that time point, then variants seen later, then those
%seen in LNMC but not Plasma, and then those seen only in PBMC
        
Data_SL8var_Plasma = Data_SL8var_Plasma_init;
Data_SL8var_LNMC = Data_SL8var_LNMC_init;
Data_SL8var_PBMC = Data_SL8var_PBMC_init;

last_var_ind = 1;%will track the index + 1 of the lowest detected variant
% at the last time point

%reording variants
for ii = 1:3
    switch ii
        case 1 %order based on Plasma
            Primary_Data = Data_SL8var_Plasma;
        case 2 %order based on LNMC
            Primary_Data = Data_SL8var_LNMC;
        case 3 %order based on PBMC
            Primary_Data = Data_SL8var_PBMC;
    end
    
    if isempty(Primary_Data.dpi) %checking that there is actually data
        %for this sample type
        continue
        
    end
    
    num_timepoints = length(Primary_Data.dpi); %number of time points
    
    
    for kk = 1:num_timepoints
        
        switch ii %getting variant count above LOD data for ordering
            case 1 %order based on Plasma
                variant_data_inform_order = ...
                    max(Data_SL8var_Plasma.V.count(1,:,kk) - Data_SL8var_Plasma.V.LOD(1,:,kk), 0);
            case 2 %order based on LNMC
                variant_data_inform_order = ...
                    max(Data_SL8var_LNMC.V.count(1,:,kk) - Data_SL8var_LNMC.V.LOD(1,:,kk),0);
            case 3 %order based on PBMC
                variant_data_inform_order = ...
                    max(Data_SL8var_PBMC.V.count(1,:,kk) - Data_SL8var_PBMC.V.LOD(1,:,kk),0);
        end
        
        %barcode ranking
        [var_tots_it, I_var_it] = sort(variant_data_inform_order(last_var_ind:end), 'descend');
        %sorting variants seen at this time point and getting the indices
        %keyboard
        
        
        for ll = 1:3 %looping over Plamsa, LNMC, and PBMC to be reordered
            
            switch ll
                case 1 %reorder Plasma
                    reorder_Data = Data_SL8var_Plasma;
                case 2 %reorder LNMC
                    reorder_Data = Data_SL8var_LNMC;
                case 3 %reorder PBMC
                    reorder_Data = Data_SL8var_PBMC;
            end
            
            if isempty(reorder_Data.dpi) %checking that there is actually data
                %for this sample type
                
                continue
            end
            
            reorder_Data.SL8_ID(:,last_var_ind:end) = ...
                reorder_Data.SL8_ID(:,last_var_ind-1 + I_var_it);
            %update order of SL8_ID
            
            num_timepoints_reorder = length(reorder_Data.dpi); %number of time points
            
            for mm = 1:num_timepoints_reorder
                
                for nn = 1:3
                    switch nn
                        case 1 %reordering counts
                            V_data = reorder_Data.V.count;
                            W_var_data = reorder_Data.W_var.count;
                            
                        case 2 %reordering LOD
                            V_data = reorder_Data.V.LOD;
                            W_var_data = reorder_Data.W_var.LOD;
                            
                        case 3 %reordering LOD_ind
                            V_data = reorder_Data.V.LOD_ind;
                            W_var_data = reorder_Data.W_var.LOD_ind;
                    end
                    
                    %Total Variants
                    V_data(1, last_var_ind:end, mm) = ...
                        V_data(1, last_var_ind-1 + I_var_it, mm);
                
                    % Barcode-SL8 var reordered by variant
                    W_var_data(:, last_var_ind:end, mm) = ...
                        W_var_data(:, last_var_ind-1 + I_var_it, mm);
                    
                    switch nn
                        case 1 %reordering counts
                            reorder_Data.V.count = V_data;
                            reorder_Data.W_var.count = W_var_data;
                            
                        case 2 %reordering LOD
                            reorder_Data.V.LOD = V_data;
                            reorder_Data.W_var.LOD = W_var_data;
                            
                        case 3 %reordering LOD_ind
                            reorder_Data.V.LOD_ind = V_data;
                            reorder_Data.W_var.LOD_ind = W_var_data;
                    end
                    
                end
                
            end
            
            
            switch ll
                case 1 %reorder Plasma
                    Data_SL8var_Plasma = reorder_Data;
                case 2 %reorder LNMC
                    Data_SL8var_LNMC = reorder_Data;
                case 3 %reorder PBMC
                    Data_SL8var_PBMC = reorder_Data;
            end
            
        end
        
        
        if (sum(var_tots_it) > 0)
            last_var_ind = find(var_tots_it, 1, 'last') + last_var_ind;
        end %updating last_var_ind if anymore variants were detected,
        %otherwise leaving it what it was before.

    end
    
end

%% getting number of mutations in each variant and then ordering by that
num_mut = sum(Data_SL8var_Plasma.SL8_ID(2:end,:) > 0);
if isempty(num_mut)
    disp('ERROR: number of mutations in variant not calculated');
end
last_var_ind = 1;%will track the index + 1 of the lowest variant reordered

for jj = 0:8
    var_num_mut_it = find(num_mut == jj); %getting indices for variants
    %with number of mutations as current iteration
    
    if isempty(var_num_mut_it)
        continue
    end
    
    var_reorder_ind_it = [1:(last_var_ind-1), var_num_mut_it,...
        setdiff(last_var_ind:num_SL8_var, var_num_mut_it)];
    %reordering variants. variants with fewer than the current number of
    %mutations stay the same (1:(last_var_ind-1)), then the variants with
    %the current number of mutations in the same order as previously 
    %(var_mut_it), and then all variants not reordered yet, but in their
    %previous order (setdiff(last_var_ind:num_SL8_var, var_num_mut_it))
    
    
    for ll = 1:3 %looping over Plamsa, LNMC, and PBMC to be reordered
            
        switch ll
            case 1 %reorder Plasma
                reorder_Data = Data_SL8var_Plasma;
            case 2 %reorder LNMC
                reorder_Data = Data_SL8var_LNMC;
            case 3 %reorder PBMC
                reorder_Data = Data_SL8var_PBMC;
        end
    
        if isempty(reorder_Data.dpi)
            continue
        end
    
        for nn = 1:3
            switch nn
                case 1 %reordering counts
                    V_data = reorder_Data.V.count;
                    W_var_data = reorder_Data.W_var.count;
                    
                case 2 %reordering LOD
                    V_data = reorder_Data.V.LOD;
                    W_var_data = reorder_Data.W_var.LOD;
                    
                case 3 %reordering LOD_ind
                    V_data = reorder_Data.V.LOD_ind;
                    W_var_data = reorder_Data.W_var.LOD_ind;
            end
            
            V_data(1, :,:) = V_data(1, var_reorder_ind_it,:);
            W_var_data(:, :,:) = W_var_data(:, var_reorder_ind_it,:);
            
            switch nn
                case 1 %reordering counts
                    reorder_Data.V.count = V_data;
                    reorder_Data.W_var.count = W_var_data;

                case 2 %reordering LOD
                    reorder_Data.V.LOD = V_data;
                    reorder_Data.W_var.LOD = W_var_data;

                case 3 %reordering LOD_ind
                    reorder_Data.V.LOD_ind = V_data;
                    reorder_Data.W_var.LOD_ind = W_var_data;
            end
        end
        
        reorder_Data.SL8_ID = reorder_Data.SL8_ID(:,var_reorder_ind_it);
        
        switch ll
            case 1 %reorder Plasma
                Data_SL8var_Plasma = reorder_Data;
            case 2 %reorder LNMC
                Data_SL8var_LNMC = reorder_Data;
            case 3 %reorder PBMC
                Data_SL8var_PBMC = reorder_Data;
        end
    end
    
    num_mut = num_mut(var_reorder_ind_it); %reordering num_mut to account 
    %for reordering of this iteration
    
    last_var_ind = last_var_ind + length(var_num_mut_it); %This updates the
    %tracking of the index + 1 of the lowest variant reordered
    
    
end

if (last_var_ind ~= length(num_mut)+1)
    disp('ERROR: number of variants reordered by number of mutations does not match')
end
    
    
