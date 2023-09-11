%Steffen Docken
%01-07-21
%This code gathers information about VL WT composition around CD8 depletion
% and the barcodes present before and after

clear
close all

load '../UPenn_FTY_barcode_SL8_Mutation_data_210607.mat'

addpath('../Data_importation/');
addpath('../General_Functions/');

save_R_data = 1; %indicator of if data saved for R (0 = no, 1 = yes)

Num_Animals = length(Data_SL8_Plasma);

%% plotting return of WT after CD8 depletion (only plotting animals not 
%necropsied before CD8 depletion)

VL_frac_WT_post_CD8depl.old_react = [];
VL_frac_WT_post_CD8depl.new_react = []; %this struct will hold the vector 
%of fraction WT in total viral load AFTER (not on same day) CD8 depletion
%for animals that had already reactivated (first vector) and the animal
%that had not previously reactivated (second vector)


PrePost_CD8depl_R_data_export_arrays = {'Animal', 'ATI2 timing', 'Barcode ID', 'Pre vs Post',...
    'Last Detection', 'WT frac at Last Det'};

num_Rexport_data = 1; %will hold current length of Rexport data


for ii = 1:Num_Animals
    VL_data_it = VL_data(ii);
    
    if (strcmp(VL_data_it.Group, 'A'))
        ART_days = Group_A_ART_days;
    else
        ART_days = Group_B_ART_days;
    end
    
    if (~strcmp(VL_data_it.animal_ID, Data_SL8_Plasma(ii).animal_ID))
        disp('ERROR: VL and SL8 animals do not match up');
        keyboard
    end
    
    [Data_SL8_Plasma_animal, ~, ~] = SL8_variant_count_210607(Data_SL8_Plasma(ii),...
        Data_SL8_LNMC(ii), Data_SL8_PBMC(ii)); %Data strucs for current
    %animal. will only contain SL8 amino acid variants detected in Plasma 
    % above LOD in at least one time point
    
    seq_ind_ATI2 = find(Data_SL8_Plasma_animal.dpi >= ART_days(4));
    %sequencing time points after second ATI

    if (isempty(VL_data_it.nec_dpi))
        cens_time = VL_data_it.CD8depl_dpi; %if not necropsied, censor at 
        %CD8 depletion

        seq_ind_early_ATI2 = seq_ind_ATI2(Data_SL8_Plasma_animal.dpi(seq_ind_ATI2) <= ...
            VL_data_it.CD8depl_dpi);
        %indices of sequencing timepoints during ATI2 before CD8
        %depletion

        seq_ind_postCD8depl = setdiff(seq_ind_ATI2, seq_ind_early_ATI2);
        %indices of postCD8 depletion sequencing

    else
        cens_time = VL_data_it.nec_dpi; %if necropsied, censor at necropsy

        seq_ind_early_ATI2 = seq_ind_ATI2; %animal necropsied, so all 
        % sequencings are 'pre-CD8 depletion'
        seq_ind_postCD8depl = [];

    end

    [react_dpi, react_ind, react_det] = find_react_dpi(VL_data_it, ...
        ART_days(4), cens_time); %finding time of reactivation in ATI-2

    %% recording fraction of VL that is WT after CD8 depletion
    if (~isempty(seq_ind_postCD8depl)) %only if sequencing happened after 
        % CD8 depletion
        first_postCD8_seq_ind = seq_ind_postCD8depl(1);

        [~, frac_VL_WT_postCD8depl] = ...
            background_subtraction(Data_SL8_Plasma_animal.V.count(1,1,first_postCD8_seq_ind),...
            Data_SL8_Plasma_animal.V.LOD(1,1,first_postCD8_seq_ind), Data_SL8_Plasma_animal.S(first_postCD8_seq_ind)); 
        %recording fraction WT above background right after CD8 depletion.

        if (react_det == 1) %reactivation detected before CD8 depletion
            VL_frac_WT_post_CD8depl.old_react = ...
                [VL_frac_WT_post_CD8depl.old_react; ...
                [Data_SL8_Plasma_animal.dpi(first_postCD8_seq_ind) - VL_data_it.CD8depl_dpi,...
                max(0, frac_VL_WT_postCD8depl)]];
                %days post CD8 depletion and fraction WT in VL in animal 
                % previously reactivated
        else
            VL_frac_WT_post_CD8depl.new_react = ...
                [VL_frac_WT_post_CD8depl.new_react; ...
                [Data_SL8_Plasma_animal.dpi(first_postCD8_seq_ind) - VL_data_it.CD8depl_dpi,...
                max(0, frac_VL_WT_postCD8depl)]];
                %days post CD8 depletion and fraction WT in VL in animal 
                % that did not reactivate pre-CD8 depletion

        end

    end

    %% getting last detection time point and corresponding percent WT for
    % barcodes detected during ATI-2

    %obtaining list of barcodes first detected before or after CD8
    %depletion during ATI-2

    if (isempty(seq_ind_postCD8depl)) %either no sequencing post CD8 
        % depletion or animal was necropsied pre-depletion

        Early_ATI2_bcode_ind = find(sum(Data_SL8_Plasma_animal.B.LOD_ind(:,...
            1,seq_ind_early_ATI2), 3));
        %indices for barcodes detected during ATI-2 before CD8 depletion
        
        Post_CD8depl_bcode_ind = []; 
        %indices for barcodes detected after CD8 depletion

        if (isempty(VL_data_it.CD8depl_dpi))
            ATI_timing_des = "necropsy"; %to record animal was necropsied
            %before CD8 depletion
        else
            ATI_timing_des = "no post CD8depl seq";
            %to record no sequencing after CD8 depletion
        end
    elseif ((~isempty(seq_ind_postCD8depl))&&... %sequencing after CD8-depl
            (react_det == 1)&&... %reactivation detected before CD8-depl
            (~isempty(seq_ind_early_ATI2))) %sequencing was done 
        % pre-CD8 depletion

        Early_ATI2_bcode_ind = find(sum(Data_SL8_Plasma_animal.B.LOD_ind(:,...
            1,seq_ind_early_ATI2), 3));
        %indices for barcodes detected during ATI-2 before CD8 depletion

        All_bcode_post_CD8depl_ind = find(sum(Data_SL8_Plasma_animal.B.LOD_ind(:,...
            1,seq_ind_postCD8depl), 3));
        %indices for all barcodes detected after CD8 depletion

        Post_CD8depl_bcode_ind = setdiff(All_bcode_post_CD8depl_ind, ...
            Early_ATI2_bcode_ind); %barcodes detected post CD8 depletion, 
        %but not in sequencing before
        disp(strcat("Animal ", Data_SL8_Plasma_animal.animal_ID, ...
            ": last seq ", num2str(VL_data_it.CD8depl_dpi - ...
            Data_SL8_Plasma_animal.dpi(seq_ind_early_ATI2(end))), ...
            " days before CD8 depl"));

        ATI_timing_des = "React pre CD8depl; Seq pre and post";
    elseif ((~isempty(seq_ind_postCD8depl))&&... %sequencing after CD8-depl
            (react_det == 1)&&... %reactivation detected before CD8-depl
            (isempty(seq_ind_early_ATI2))) %sequencing not done 
        % pre-CD8 depletion

        Early_ATI2_bcode_ind = [];
        %indices for barcodes detected during ATI-2 before CD8 depletion
        %(none)

        Post_CD8depl_bcode_ind = find(sum(Data_SL8_Plasma_animal.B.LOD_ind(:,...
            1,seq_ind_postCD8depl), 3));
        %indices for all barcodes detected after CD8 depletion (all
        %barcodes detected in ATI-2)

        CD8depl_time_since_last_VL_det = VL_data_it.CD8depl_dpi -...
            max(VL_data_it.t((VL_data_it.VL > 60)&...
            (VL_data_it.t < VL_data_it.CD8depl_dpi)));
        %time from last VL > detection to CD8 depletion
        ATI_timing_des = strcat("No seq pre CD8depl; ", ...
            num2str(CD8depl_time_since_last_VL_det), " days since VL det");

        disp(strcat("Animal ", Data_SL8_Plasma_animal.animal_ID, ": ", ...
            ATI_timing_des));

    elseif ((~isempty(seq_ind_postCD8depl))&&... %sequencing after CD8-depl
            (react_det == 0)) %reactivation detected after CD8-depl
        
        Early_ATI2_bcode_ind = [];
        %indices for barcodes detected during ATI-2 before CD8 depletion
        %(none)

        Post_CD8depl_bcode_ind = find(sum(Data_SL8_Plasma_animal.B.LOD_ind(:,...
            1,seq_ind_postCD8depl), 3));
        %indices for all barcodes detected after CD8 depletion (all
        %barcodes detected in ATI-2)

        ATI_timing_des = "React post CD8depl";
        disp(strcat("Animal ", Data_SL8_Plasma_animal.animal_ID,...
            ": React post CD8depl"));

    else
        disp(strcat(Data_SL8_Plasma_animal.animal_ID, " not chategorized"));
    end


    for kk = 1:2 %pre and post CD8 depletion barcodes
        switch kk
            case 1 %pre-CD8 depletion barcodes
                bcode_ind_it = Early_ATI2_bcode_ind;
            case 2 %post-CD8 depletion barcodes
                bcode_ind_it = Post_CD8depl_bcode_ind;
        end

        if (isempty(bcode_ind_it))
            continue %skip iteration if there are no barcodes for this 
            % category
        end

        for jj = 1:length(bcode_ind_it)
            last_det_ind_it = find(Data_SL8_Plasma_animal.B.LOD_ind(...
                bcode_ind_it(jj), 1, 1:(seq_ind_ATI2(1)-1)), 1, ...
                'last');%index of timing of last detection of barcode 
            % before ATI-2
            
            if (isempty(last_det_ind_it))
                last_det_var_it = 0; %setting indicator for timing of last 
                %detection of this barcode to 0 to indicate it was never seen
                %before
                last_det_str_it = 'Never';
            elseif (Data_SL8_Plasma_animal.dpi(last_det_ind_it) < 15)
                last_det_var_it = 1;%setting indicator for timing of last 
                %detection of this barcode to 1 to indicate it was last seen in
                %Pre-ART
                last_det_str_it = 'Primary';
            elseif ((Data_SL8_Plasma_animal.dpi(last_det_ind_it) >= ART_days(2))&&...
                    (Data_SL8_Plasma_animal.dpi(last_det_ind_it) <= 250))
                last_det_var_it = 2;%setting indicator for timing of last 
                %detection of this barcode to 2 to indicate it was last seen in
                %Early ATI-1
                last_det_str_it = 'Early ATI-1';
            elseif ((Data_SL8_Plasma_animal.dpi(last_det_ind_it) > 250)&&...
                    (Data_SL8_Plasma_animal.dpi(last_det_ind_it) <= ART_days(3)))
                last_det_var_it = 3;%setting indicator for timing of last 
                %detection of this barcode to 3 to indicate it was last seen in
                %Late ATI-1
                last_det_str_it = 'Late ATI-1';
            else
                disp('ERROR: last detection of barcode indeterminant');
                keyboard
            end
            
            if (last_det_var_it == 0)
                bcode_WT_last_det_it = nan; %recording last detection of WT
                %as not a number if never detected
            else
                
                [~, bcode_WT_last_det_it] = background_subtraction(Data_SL8_Plasma_animal.W_var.count(bcode_ind_it(jj), 1, last_det_ind_it),...
                    Data_SL8_Plasma_animal.W_var.LOD(bcode_ind_it(jj), 1, last_det_ind_it),...
                    Data_SL8_Plasma_animal.B.count(bcode_ind_it(jj), 1, last_det_ind_it));

                bcode_WT_last_det_it = max(0, bcode_WT_last_det_it); %shifting 
            % negative values (after background subtraction) to 0.
            end
            
            
            %adding to data for R export
            PrePost_CD8depl_R_data_export_arrays{num_Rexport_data + 1, 1} = ...
                Data_SL8_Plasma_animal.animal_ID;
            PrePost_CD8depl_R_data_export_arrays{num_Rexport_data + 1, 2} = ...
                ATI_timing_des;
            PrePost_CD8depl_R_data_export_arrays{num_Rexport_data + 1, 3} = ...
                Data_SL8_Plasma_animal.barcode_ID(bcode_ind_it(jj));
            switch kk
                case 1
                    PrePost_CD8depl_R_data_export_arrays{num_Rexport_data + ...
                        1, 4} = 'Pre'; %barcode detected pre-CD8 depletion
                case 2
                    PrePost_CD8depl_R_data_export_arrays{num_Rexport_data + ...
                        1, 4} = 'Post'; %barcode detected post-CD8 depletion
            end

            PrePost_CD8depl_R_data_export_arrays{num_Rexport_data + 1, 5} = ...
                last_det_str_it;
            PrePost_CD8depl_R_data_export_arrays{num_Rexport_data + 1, 6} = ...
                bcode_WT_last_det_it;
            num_Rexport_data = num_Rexport_data + 1; %increasing count of 
            %number data points recorded
            
        end
        
    end

end

if (save_R_data == 1)
    writecell(PrePost_CD8depl_R_data_export_arrays, 'PrePost_CD8depl_barcodes_data_for_R_230120.csv');
end