%Steffen Docken
%22-2-22
%This code generates summary statistics used in the paper

clear
close all

load UPenn_FTY_barcode_SL8_Mutation_data_210607.mat

addpath('General_Functions/');
addpath('Data_importation/');

CD8_depl_start = 465; %CD8 depletion was started on this day post 
%inoculation or soon after. All VL measurements after this were after CD8 
%depletion

day_early_ATI = 10; %maximum number of days since reactivation for sequencing
% to count as as early ATI


VL_LOD = 60; %copies/ml. limit of detection for viral load

Num_animals = length(Data_SL8_Plasma);

%array to hold data:
VL_summary_table = table; %will hold animal 
VL_summary_names = {'animal', 'FTY_info', 'max_PreART_log10VL', 'day_14_log10VL', 'ART_1_suppression',...
    'ATI_1_react_timing', 'ATI_1_log10peak', 'ATI_1_setpoint',...
    'ART_2_suppression_time', 'ATI_2_react_timing', 'ATI_2_log10peak',...
    'ATI_2_setpoint', 'log10VL_preCD8depl', 'log10VL_postCD8depl', ...
    'log10VL_increase_postCD8depl', 'preCD8depl_timing', 'postCD8depl_timing'};
VL_summary_types = {'cellstr', 'cellstr', 'double', 'double', 'int8', 'int8',...
    'double', 'double', 'int8', 'int8',...
    'double', 'double', 'double', 'double', 'double', 'int8', 'int8'};


for ii = 1:Num_animals
    VL_summary_table_it = table('Size', [1, 17],...
        'VariableTypes', VL_summary_types,...
        'VariableNames', VL_summary_names);
    
    VL_data_it = VL_data(ii); %VL data for this animals
    
    switch VL_data_it.Group
        case 'A'
            ART_days_it = Group_A_ART_days;
        case 'B'
            ART_days_it = Group_B_ART_days;
    end
    
    VL_summary_table_it.animal = {VL_data_it.animal_ID};
    VL_summary_table_it.FTY_info = {VL_data_it.FTY_info};
    
    VL_summary_table_it.max_PreART_log10VL = log10(max(VL_data_it.VL(VL_data_it.t <= 14)));

    VL_summary_table_it.day_14_log10VL = log10(VL_data_it.VL(VL_data_it.t == 14));
    
    %% 2 days below detection threshold
    VL_undet_2days_ind = ((VL_data_it.VL(1:end-1) == 60)&(VL_data_it.VL(2:end)==60));
    %indicator of if VL is below detection on each day and the subsequent
    %measurement
    ART_1_suppression_ind = find((VL_undet_2days_ind)&...
        (VL_data_it.t(1:end-1) > ART_days_it(1)), 1, 'first');
    
    VL_summary_table_it.ART_1_suppression = VL_data_it.t(ART_1_suppression_ind);
    %this is dpi
    %% time to rebound (react on first day off ART counts
    %as react in 0 days)

    [react_ATI1_dpi, react_ATI1_ind, react_ATI1_det] = ...
        find_react_dpi(VL_data_it, ART_days_it(2), ART_days_it(3));

    if (react_ATI1_det == 0)
        VL_summary_table_it.ATI_1_react_timing = nan;
        VL_summary_table_it.ATI_1_log10peak = nan;
        VL_summary_table_it.ATI_1_setpoint = nan;
    else
    
        VL_summary_table_it.ATI_1_react_timing = react_ATI1_dpi -...
            ART_days_it(2);
    
        
        %% ATI-1 peak 
        ATI_1_peak_ind_range2 = find(VL_data_it.t <= VL_data_it.t(react_ATI1_ind) + 30,...
            1, 'last');
        
        VL_summary_table_it.ATI_1_log10peak = ...
            log10(max(VL_data_it.VL(react_ATI1_ind:ATI_1_peak_ind_range2)));
        
        %% ATI-1 set point: time weighted AUC of log10 VL from 30 to 60 days after rebound
        ATI_1_setpoint_ind1 = ATI_1_peak_ind_range2 + 1;
        ATI_1_setpoint_time2 = min(ART_days_it(3), VL_data_it.t(react_ATI1_ind) + 60);
        ATI_1_setpoint_ind2 = find(VL_data_it.t <= ATI_1_setpoint_time2, 1, 'last');
        
        if (ATI_1_setpoint_ind1 > ATI_1_setpoint_ind2)
            VL_summary_table_it.ATI_1_setpoint = NaN;
        elseif (ATI_1_setpoint_ind1 == ATI_1_setpoint_ind2)
            VL_summary_table_it.ATI_1_setpoint = log10(VL_data_it.VL(ATI_1_setpoint_ind1));
        else
            VL_summary_table_it.ATI_1_setpoint =...
                trapz(VL_data_it.t(ATI_1_setpoint_ind1: ATI_1_setpoint_ind2), ...
                log10(VL_data_it.VL(ATI_1_setpoint_ind1: ATI_1_setpoint_ind2)))/...
                (VL_data_it.t(ATI_1_setpoint_ind2) - VL_data_it.t(ATI_1_setpoint_ind1));
        end
    end
    
    %% 2 days below detection threshold
    ART_2_suppression_ind = find((VL_undet_2days_ind)&...
        (VL_data_it.t(1:end-1) > ART_days_it(3)), 1, 'first');
    
    VL_summary_table_it.ART_2_suppression_time = VL_data_it.t(ART_2_suppression_ind) -...
        ART_days_it(3); %number of days since reinitiation of ART to suppression
    
    %% time to rebound (react on first day off ART counts
    %as react in 0 days)
    if(isempty(VL_data_it.CD8depl_dpi))
        ATI2_end_dpi = VL_data_it.nec_dpi;
    else
        ATI2_end_dpi = VL_data_it.CD8depl_dpi;
    end

    [react_ATI2_dpi, react_ATI2_ind, react_ATI2_det] = ...
        find_react_dpi(VL_data_it, ART_days_it(4), ATI2_end_dpi);

    if (react_ATI2_det == 0)

        VL_summary_table_it.ATI_2_react_timing = nan;
        VL_summary_table_it.ATI_2_log10peak = nan;
        VL_summary_table_it.ATI_2_setpoint = nan;

    else
    
        VL_summary_table_it.ATI_2_react_timing = VL_data_it.t(react_ATI2_ind) -...
            ART_days_it(4);
        
        %% ATI-2 peak
        ATI_2_peak_ind_range2 = find(VL_data_it.t <= VL_data_it.t(react_ATI2_ind) + 30,...
            1, 'last');
        
        VL_summary_table_it.ATI_2_log10peak = log10(max(VL_data_it.VL(react_ATI2_ind:ATI_2_peak_ind_range2)));
        
        %% ATI-2 set point: time weighted AUC of log10 VL from 30 to 60 days after rebound
        ATI_2_setpoint_ind1 = ATI_2_peak_ind_range2 + 1;
            
        ATI_2_setpoint_ind2 = find(VL_data_it.t <= min(ATI2_end_dpi,...
            VL_data_it.t(react_ATI2_ind) + 60), 1, 'last');
        %finding last time point for set point calculation. Either minimum
        %of necropsy day, CD8 depletion day, or 60 days after detected
        %reactivation
    
        if (ATI_2_setpoint_ind1 > ATI_2_setpoint_ind2)
            VL_summary_table_it.ATI_2_setpoint = NaN; %no relevant data >30 
            % days after reactivation
        elseif (ATI_2_setpoint_ind1 == ATI_2_setpoint_ind2)
            VL_summary_table_it.ATI_2_setpoint = log10(VL_data_it.VL(ATI_2_setpoint_ind1));
            %only 1 relevant day between 30 and 60 days post reactivation
        else
            VL_summary_table_it.ATI_2_setpoint =...
                trapz(VL_data_it.t(ATI_2_setpoint_ind1: ATI_2_setpoint_ind2), ...
                log10(VL_data_it.VL(ATI_2_setpoint_ind1: ATI_2_setpoint_ind2)))/...
                (VL_data_it.t(ATI_2_setpoint_ind2) - VL_data_it.t(ATI_2_setpoint_ind1));
        end

    end
    
    %% VL pre-CD8 depletion
    if (isempty(VL_data_it.CD8depl_dpi))
        VL_summary_table_it.log10VL_preCD8depl = nan;
        VL_summary_table_it.log10VL_postCD8depl = nan;
        VL_summary_table_it.log10VL_increase_postCD8depl = nan;

        VL_summary_table_it.preCD8depl_timing = nan;
        VL_summary_table_it.postCD8depl_timing = nan;
        
    else
        
        preCD8depl_ind = find(VL_data_it.t <= VL_data_it.CD8depl_dpi, 1, 'last');
        VL_summary_table_it.log10VL_preCD8depl = log10(VL_data_it.VL(preCD8depl_ind));
        VL_summary_table_it.preCD8depl_timing = VL_data_it.t(preCD8depl_ind)-...
            VL_data_it.CD8depl_dpi;
        
        postCD8depl_ind = find(VL_data_it.t > VL_data_it.CD8depl_dpi, 1, 'first');
        VL_summary_table_it.log10VL_postCD8depl = log10(VL_data_it.VL(postCD8depl_ind));
        VL_summary_table_it.postCD8depl_timing = VL_data_it.t(postCD8depl_ind)-...
            VL_data_it.CD8depl_dpi;
        
        VL_summary_table_it.log10VL_increase_postCD8depl = ...
            VL_summary_table_it.log10VL_postCD8depl-...
            VL_summary_table_it.log10VL_preCD8depl;
    end
    
    %% adding to table
    if (ii == 1)
        VL_summary_table = VL_summary_table_it;
    else
        VL_summary_table = [VL_summary_table; VL_summary_table_it];
    end
end


ATI_peak_ind = (~isnan(VL_summary_table.ATI_1_log10peak)&...
    ~isnan(VL_summary_table.ATI_2_log10peak));
[p_ATI_peak, h_ATI_peak] = signrank(VL_summary_table.ATI_1_log10peak(ATI_peak_ind),...
    VL_summary_table.ATI_2_log10peak(ATI_peak_ind));

disp(strcat("Comparing log10 peak VL of ATI-1 vs. ATI-2: p-value = ", ...
    num2str(p_ATI_peak), "; paired Wilcoxon signed rank test"));

ATI_setpoint_ind = (~isnan(VL_summary_table.ATI_1_setpoint)&...
    ~isnan(VL_summary_table.ATI_2_setpoint));
[p_ATI_setpoint, h_ATI_setpoint] = signrank(VL_summary_table.ATI_1_setpoint(ATI_setpoint_ind),...
    VL_summary_table.ATI_2_setpoint(ATI_setpoint_ind));
disp(strcat("Comparing time weighted log10 set point VL of ATI-1 vs. ATI-2: p-value = ", ...
    num2str(p_ATI_setpoint), "; paired Wilcoxon signed rank test"));

%% Protocol

%array to hold data:
Protocol_summary_table = table; %will hold animal 
Protocol_summary_names = {'animal', 'ART1_duration', 'ATI1_duration', ...
    'ART2_duration', 'ATI2_to_CD8depl', 'LNMC_ART1', 'LNMC_ATI1', 'LNMC_ART2', 'LNMC_ATI2', ...
    'PBMC_ART1', 'PBMC_ATI1', 'PBMC_ART2', 'PBMC_ATI2'};
Protocol_summary_types = {'cellstr', 'int8', 'int8', 'int8', 'int8',...
    'int8', 'int8', 'int8', 'int8', 'int8', 'int8', 'int8', 'int8'};


for ii = 1:Num_animals
    Protocol_summary_table_it = table('Size', [1, length(Protocol_summary_names)],...
        'VariableTypes', Protocol_summary_types,...
        'VariableNames', Protocol_summary_names);
    
    VL_data_it = VL_data(ii); %VL data for this animals

    switch VL_data_it.Group
        case 'A'
            ART_days_it = Group_A_ART_days;
        case 'B'
            ART_days_it = Group_B_ART_days;
    end

    Protocol_summary_table_it.animal = {VL_data_it.animal_ID};

    Protocol_summary_table_it.ART1_duration = ART_days_it(2)-ART_days_it(1);
    %Number of days of first round of ART: first day off ART
    %(ART_days_it(2)) - first day of ART (ART_days_it(1))

    Protocol_summary_table_it.ATI1_duration = ART_days_it(3)-ART_days_it(2);
    %Number of days of first ATI: first day of second round of ART
    %(ART_days_it(3)) - first day without ART (ART_days_it(2))

    Protocol_summary_table_it.ART2_duration = ART_days_it(4)-ART_days_it(3);
    %Number of days of second round of ART: first day off second round of ART
    %(ART_days_it(4)) - first day of second round of ART (ART_days_it(3))

    if (isempty(VL_data_it.CD8depl_dpi))
        Protocol_summary_table_it.ATI2_to_CD8depl = NaN;
    else
        Protocol_summary_table_it.ATI2_to_CD8depl = VL_data_it.CD8depl_dpi - ...
            ART_days_it(4);
    end

    LNMC_SL8_Data_it = Data_SL8_LNMC(ii);
    PBMC_SL8_Data_it = Data_SL8_PBMC(ii);
    if ((~isempty(LNMC_SL8_Data_it.animal_ID) && (~strcmp(VL_data_it.animal_ID, LNMC_SL8_Data_it.animal_ID)))||...
            (~isempty(PBMC_SL8_Data_it.animal_ID) && (~strcmp(VL_data_it.animal_ID, PBMC_SL8_Data_it.animal_ID))))
        disp('ERROR: not all same animals')
        keyboard
    end

    for jj = 1:4
        switch jj
            case 1
                phase_it = 'ART-1';
            case 2
                phase_it = 'ATI-1';
            case 3
                phase_it = 'ART-2';
            case 4
                phase_it = 'ATI-2';
        end

        %% LNMC
        LNMC_phase_it_ind = find(strcmp(LNMC_SL8_Data_it.phase, phase_it));

        if (isempty(LNMC_phase_it_ind))
            LNMC_dpi_it = NaN;
        else
            LNMC_dpi_it = LNMC_SL8_Data_it.dpi(LNMC_phase_it_ind(1));
        end
        if (length(LNMC_phase_it_ind) > 1)
            disp(strcat("ERROR: multiple time points for LNMC in animal ", ...
                Protocol_summary_table_it.animal));
        end

        %% PBMC
        PBMC_phase_it_ind = find(strcmp(PBMC_SL8_Data_it.phase, phase_it));

        if (isempty(PBMC_phase_it_ind))
            PBMC_dpi_it = NaN;
        else
            PBMC_dpi_it = PBMC_SL8_Data_it.dpi(PBMC_phase_it_ind(1));
        end
        if (length(PBMC_phase_it_ind) > 1)
            disp(strcat("ERROR: multiple time points for PBMC in animal ", ...
                Protocol_summary_table_it.animal));
        end


        switch jj
            case 1
                Protocol_summary_table_it.LNMC_ART1 = LNMC_dpi_it;
                Protocol_summary_table_it.PBMC_ART1 = PBMC_dpi_it;

            case 2
                Protocol_summary_table_it.LNMC_ATI1 = LNMC_dpi_it;
                Protocol_summary_table_it.PBMC_ATI1 = PBMC_dpi_it;

            case 3
                Protocol_summary_table_it.LNMC_ART2 = LNMC_dpi_it;
                Protocol_summary_table_it.PBMC_ART2 = PBMC_dpi_it;

            case 4
                Protocol_summary_table_it.LNMC_ATI2 = LNMC_dpi_it;
                Protocol_summary_table_it.PBMC_ATI2 = PBMC_dpi_it;
        end
    end

    %% adding to table
    if (ii == 1)
        Protocol_summary_table = Protocol_summary_table_it;
    else
        Protocol_summary_table = [Protocol_summary_table; Protocol_summary_table_it];
    end
end

%% WT kinetics and largest non-WT variant

%array to hold data:
WT_summary_table = table; 
WT_summary_names = {'animal', 'FTY_info', 'WT_day14', 'WT_LNMC_ART1', 'WT_PBMC_ART1', ...
    'WT_early_ATI1', 'early_ATI1_timing', 'WT_late_ATI1', 'late_ATI1_dpi',...
    'late_ATI1_timing', 'WT_LNMC_ART2', 'WT_PBMC_ART2', 'WT_early_ATI2',...
    'early_ATI2_timing', 'WT_postCD8depl', 'days_postCD8depl'};
WT_summary_types = {'cellstr', 'cellstr', 'double', 'double', 'double', 'double',...
    'int8', 'double', 'int8', 'int8', 'double', 'double', 'double',...
    'int8', 'double', 'int8'};

%WT sequence
WT_SL8_seq = 'STPESANL';

%array to hold data:
nonWT_vars_geq20percent = cell(1,2);
nonWT_vars_geq20percent{1,1} = 'SL8ID';
nonWT_vars_geq20percent{1,2} = 'SL8 aa seq';
nonWT_vars_geq20percent_vec = [];


for ii = 1:Num_animals
    WT_summary_table_it = table('Size', [1, length(WT_summary_names)],...
        'VariableTypes', WT_summary_types,...
        'VariableNames', WT_summary_names);
    
    VL_data_it = VL_data(ii); %VL data for this animals

    [Plasma_SL8_Data_it, LNMC_SL8_Data_it, PBMC_SL8_Data_it] = ...
        SL8_variant_count_210607(Data_SL8_Plasma(ii), Data_SL8_LNMC(ii), ...
        Data_SL8_PBMC(ii)); %SL8 data for this animal

    switch VL_data_it.Group
        case 'A'
            ART_days_it = Group_A_ART_days;
        case 'B'
            ART_days_it = Group_B_ART_days;
    end

    if ((~isempty(LNMC_SL8_Data_it.animal_ID) && (~strcmp(VL_data_it.animal_ID, LNMC_SL8_Data_it.animal_ID)))||...
            (~isempty(PBMC_SL8_Data_it.animal_ID) && (~strcmp(VL_data_it.animal_ID, PBMC_SL8_Data_it.animal_ID))))
        disp('ERROR: not all same animals')
        keyboard
    end

    WT_summary_table_it.animal = {VL_data_it.animal_ID};
    
    WT_summary_table_it.FTY_info = {VL_data_it.FTY_info};

    %% day 14 WT
    day_14_ind = find(Plasma_SL8_Data_it.dpi == 14);
    [~, WT_summary_table_it.WT_day14] = background_subtraction(Plasma_SL8_Data_it.V.count(1, 1, day_14_ind),...
        Plasma_SL8_Data_it.V.LOD(1, 1, day_14_ind), ...
        Plasma_SL8_Data_it.S(day_14_ind));

    %% ART DNA
    for kk = 1:2
        for jj = 1:2
            switch jj
                case 1
                    SL8_Data_it = LNMC_SL8_Data_it;
                case 2
                    SL8_Data_it = PBMC_SL8_Data_it;
            end
            
            phase_it_ind = find(strcmp(SL8_Data_it.phase, ...
                strcat('ART-', num2str(kk))));
    
            if (isempty(phase_it_ind))
                WT_prop_it = NaN;
            else
                [~, WT_prop_it] = background_subtraction(SL8_Data_it.V.count(1, 1, phase_it_ind),...
                    SL8_Data_it.V.LOD(1, 1, phase_it_ind), ...
                    SL8_Data_it.S(phase_it_ind));

                WT_prop_it = max(WT_prop_it, 0); %making sure WT prop is >0
            end %getting proportion WT after background subtraction

            if (length(phase_it_ind) > 1)
                disp(strcat("ERROR: multiple time points for DNA in animal ", ...
                    Protocol_summary_table_it.animal));
            end
    
            switch kk
                case 1 %ART-1
                    switch jj
                        case 1
                            WT_summary_table_it.WT_LNMC_ART1 = WT_prop_it;
                        case 2
                            WT_summary_table_it.WT_PBMC_ART1 = WT_prop_it;
                    end
                case 2 %ART-2
                    switch jj
                        case 1
                            WT_summary_table_it.WT_LNMC_ART2 = WT_prop_it;
                        case 2
                            WT_summary_table_it.WT_PBMC_ART2 = WT_prop_it;
                    end
            end
        end

    end

    %% early ATI
    
    for jj = 1:2

        if (isempty(VL_data_it.CD8depl_dpi))
            phase_it_ind = find(strcmp(Plasma_SL8_Data_it.phase, ...
                strcat('ATI-', num2str(jj))), 1);
        else
            phase_it_ind = find(strcmp(Plasma_SL8_Data_it.phase, ...
                strcat('ATI-', num2str(jj)))&...
                Plasma_SL8_Data_it.dpi <= VL_data_it.CD8depl_dpi, 1);
        end
        
        if (isempty(phase_it_ind))
            WT_prop_it = NaN;
            timing_it = NaN;
            num_new_var_it = NaN;
            double_mut_on_dom_it = NaN;

        else
            [~, WT_prop_it] = background_subtraction(Plasma_SL8_Data_it.V.count(1, 1, phase_it_ind),...
                Plasma_SL8_Data_it.V.LOD(1, 1, phase_it_ind), ...
                Plasma_SL8_Data_it.S(phase_it_ind));

            WT_prop_it = max(WT_prop_it, 0);

            timing_it = Plasma_SL8_Data_it.timing(phase_it_ind);
            %days since first detection that sequencing occured (if
            %sequencing is same day as detection, this is written as 1)

        end

        switch jj
            case 1
                WT_summary_table_it.WT_early_ATI1 = WT_prop_it;
                WT_summary_table_it.early_ATI1_timing = timing_it;

            case 2
                WT_summary_table_it.WT_early_ATI2 = WT_prop_it;
                WT_summary_table_it.early_ATI2_timing = timing_it;

        end


    end

    %% late ATI-1 info
    phase_it_ind = find(Plasma_SL8_Data_it.dpi == ART_days_it(3));% animals
        % with sequencing time point at resumption of ART

    if (isempty(phase_it_ind)) 
        WT_summary_table_it.WT_late_ATI1 = NaN;
        WT_summary_table_it.late_ATI1_dpi = NaN;
        WT_summary_table_it.late_ATI1_timing = NaN;

    else
        [~, WT_prop_it] = background_subtraction(Plasma_SL8_Data_it.V.count(1, 1, phase_it_ind),...
            Plasma_SL8_Data_it.V.LOD(1, 1, phase_it_ind), ...
            Plasma_SL8_Data_it.S(phase_it_ind));

        WT_prop_it = max(WT_prop_it, 0); %ensuring WT prop is >0

        WT_summary_table_it.WT_late_ATI1 = WT_prop_it;
        WT_summary_table_it.late_ATI1_timing = Plasma_SL8_Data_it.timing(phase_it_ind);
        %days since first detection that sequencing occured (if
        %sequencing is same day as detection, this is written as 1)

        WT_summary_table_it.late_ATI1_dpi = Plasma_SL8_Data_it.dpi(phase_it_ind);
    end

    %% post CD8-delpetion WT percent
    if(~isempty(VL_data_it.CD8depl_dpi))
        postCD8_depl_ind = find(Plasma_SL8_Data_it.dpi >VL_data_it.CD8depl_dpi, ...
            1, 'first');

        if(isempty(postCD8_depl_ind))
            WT_summary_table_it.WT_postCD8depl = nan;
            WT_summary_table_it.days_postCD8depl = nan;
        else
            [~, WT_prop_it] = background_subtraction(Plasma_SL8_Data_it.V.count(1, 1, postCD8_depl_ind),...
                Plasma_SL8_Data_it.V.LOD(1, 1, postCD8_depl_ind), ...
                Plasma_SL8_Data_it.S(postCD8_depl_ind));
    
            WT_summary_table_it.WT_postCD8depl = max(WT_prop_it, 0);
            WT_summary_table_it.days_postCD8depl = ...
                Plasma_SL8_Data_it.dpi(postCD8_depl_ind) - ...
                VL_data_it.CD8depl_dpi;
            %days since first detection that sequencing occured (if
            %sequencing is same day as detection, this is written as 1)
        end
    end

     %% adding to table
    if (ii == 1)
        WT_summary_table = WT_summary_table_it;

    else
        WT_summary_table = [WT_summary_table; WT_summary_table_it];

    end

    %% non-WT vars greater than 20 percent at some point.
    for jj = 1:length(Plasma_SL8_Data_it.dpi)
        [~, var_frac_it] = background_subtraction(Plasma_SL8_Data_it.V.count(1, :, jj),...
            Plasma_SL8_Data_it.V.LOD(1, :, jj), ...
            Plasma_SL8_Data_it.S(jj));

        var_geq20percent_ind_all_it = find(var_frac_it >= 0.2);
        var_geq20percent_SL8_all_it = Plasma_SL8_Data_it.SL8_ID(1, ...
                    var_geq20percent_ind_all_it);
        var_geq20percent_SL8_it = setdiff(var_geq20percent_SL8_all_it, 0); 
        %removing WT from set of variants

        if(~isempty(var_geq20percent_SL8_it))

            for kk = 1:length(var_geq20percent_SL8_it)
                SL8_ID_it = var_geq20percent_SL8_it(kk); %current SL8
                
                if(ismember(SL8_ID_it, nonWT_vars_geq20percent_vec))
                    continue %skip iteration if current variant already 
                    % recorded
                end

                nonWT_vars_geq20percent_vec = [nonWT_vars_geq20percent_vec;...
                    SL8_ID_it];
                %tracking non-WT variants that have been recorded.

                SL8_seq_it = WT_SL8_seq; %redefining current SL8 sequence to WT 
                %(which will be edited)
                
                SL8_aa_seq_code_it = Plasma_SL8_Data_it.SL8_ID(2:9, ...
                    Plasma_SL8_Data_it.SL8_ID(1,:) == SL8_ID_it);
                
                for ll = 1:8
                    if(SL8_aa_seq_code_it(ll) == 0)
                        continue %this aa is wild-type
                    end
        
                    SL8_seq_it(ll) = aa_translation_210603(SL8_aa_seq_code_it(ll));
                end

                nonWT_vars_geq20percent{end+1, 1} = SL8_ID_it;
                nonWT_vars_geq20percent{end, 2} = SL8_seq_it;

            end
        end

    end
end


early_ATI1_seq_ind = (~isnan(WT_summary_table.early_ATI1_timing)&...
    (WT_summary_table.early_ATI1_timing <= day_early_ATI+1)); %only including
%sequencing done "early" after detected rebound. The '+1' is because
%rebound timing is recoreded with day of detection as day 1.
[p_WT_early_ATI1_preART, h_WT_early_ATI1_preART] = signrank(WT_summary_table.WT_day14(early_ATI1_seq_ind),...
    WT_summary_table.WT_early_ATI1(early_ATI1_seq_ind));
disp(strcat("Testing if percent WT is different at day 14 vs. early ATI-1: p-value = ",...
    num2str(p_WT_early_ATI1_preART), "; paired Wilcoxon signed rank test"));

[p_FTY_lateATI1_WT,h_FTY_lateATI1_WT] = ranksum(WT_summary_table.WT_late_ATI1(strcmp(WT_summary_table.FTY_info, "FTY")),...
    WT_summary_table.WT_late_ATI1(strcmp(WT_summary_table.FTY_info, "Control")));
disp(strcat("Testing if percent WT late in ATI-1 is different for FTY and Control animals: p-value = ",...
    num2str(p_FTY_lateATI1_WT), "; Wilcoxon rank sum test"));

early_ATI1_ATI2_seq_ind = (~isnan(WT_summary_table.early_ATI1_timing)&...
    (WT_summary_table.early_ATI1_timing <= day_early_ATI+1)&...
    (~isnan(WT_summary_table.early_ATI2_timing))&...
    (WT_summary_table.early_ATI2_timing <= day_early_ATI+1));%only including
%sequencing done "early" after detected rebound. The '+1' is because
%rebound timing is recoreded with day of detection as day 1.
[p_WT_early_ATI1_ATI2, h_WT_early_ATI1_ATI2] = signrank(WT_summary_table.WT_early_ATI2(early_ATI1_ATI2_seq_ind),...
    WT_summary_table.WT_early_ATI1(early_ATI1_ATI2_seq_ind));
disp(strcat("Testing if percent WT is different at early ATI-1 vs. early ATI-2: p-value = ",...
    num2str(p_WT_early_ATI1_ATI2), "; paired Wilcoxon signed rank test"));


WT_Comparisons_table = table;
WT_Comparisons_table.Comparisons = {'diff_LNMC_ART1_day14'; 'diff_PBMC_ART1_day14';...
    'diff_LNMC_ART2_lateATI1'; 'diff_PBMC_ART2_lateATI1'};
WT_Comparisons_table.min = [min(WT_summary_table.WT_LNMC_ART1 - WT_summary_table.WT_day14);...
    min(WT_summary_table.WT_PBMC_ART1 - WT_summary_table.WT_day14);...
    min(WT_summary_table.WT_LNMC_ART2 - WT_summary_table.WT_late_ATI1);...
    min(WT_summary_table.WT_PBMC_ART2 - WT_summary_table.WT_late_ATI1)];

WT_Comparisons_table.max = [max(WT_summary_table.WT_LNMC_ART1 - WT_summary_table.WT_day14);...
    max(WT_summary_table.WT_PBMC_ART1 - WT_summary_table.WT_day14);...
    max(WT_summary_table.WT_LNMC_ART2 - WT_summary_table.WT_late_ATI1);...
    max(WT_summary_table.WT_PBMC_ART2 - WT_summary_table.WT_late_ATI1)];

%% CD8 levels

%array to hold data:
TL8CD8_summary_table = table; 
TL8CD8_summary_names = {'animal', 'Treatment', 'PBMC_day14_TL8', 'PBMC_day14_day',...
    'PBMC_ART1_TL8', 'PBMC_ART1_day', 'PBMC_ATI1_start_TL8', 'PBMC_ATI1_start_day', ...
    'PBMC_ATI1_early_TL8', 'PBMC_ATI1_early_day', 'PBMC_ATI1_end_TL8', 'PBMC_ATI1_end_day',...
    'PBMC_ART2_early_TL8', 'PBMC_ART2_early_day', 'PBMC_ART2_end_TL8', 'PBMC_ART2_end_day', ...
    'LNMC_ART1_TL8', 'LNMC_ART1_day', 'LNMC_ATI1_start_TL8', 'LNMC_ATI1_start_day',...
    'LNMC_ATI1_TL8', 'LNMC_ATI1_day', 'LNMC_ART2_early_TL8', 'LNMC_ART2_early_day'};
TL8CD8_summary_types = {'cellstr', 'cellstr', 'double', 'double', 'double', 'double',...
     'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', ...
     'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', ...
     'double', 'double'};

PBMC_TL8_totCD8_data = readtable('Data/3-9-22 TL8-PBMC-LNMC_SSD_names_corrected.xlsx', ...
    'Sheet', 'TL8 - PBMC', 'Range', 'B3');
LNMC_TL8_totCD8_data = readtable('Data/3-9-22 TL8-PBMC-LNMC_SSD_names_corrected.xlsx', ...
    'Sheet', 'TL8 - LNMC', 'Range', 'B3');

for ii = 1:Num_animals
    TL8CD8_summary_table_it = table('Size', [1, length(TL8CD8_summary_names)],...
        'VariableTypes', TL8CD8_summary_types,...
        'VariableNames', TL8CD8_summary_names);
    
    VL_data_it = VL_data(ii); %VL data for this animals

    TL8CD8_summary_table_it.animal = {VL_data_it.animal_ID};

    TL8CD8_summary_table_it.Treatment = {VL_data_it.FTY_info};

    %PBMC data
    PBMC_TL8_totCD8_data_it = PBMC_TL8_totCD8_data(strcmp(PBMC_TL8_totCD8_data.MonkeyID,...
                    VL_data_it.animal_ID),:);

    day_14_ind_it = find(PBMC_TL8_totCD8_data_it.DaysP_i_ == 14);
    TL8CD8_summary_table_it.PBMC_day14_day = PBMC_TL8_totCD8_data_it.DaysP_i_(day_14_ind_it);
    TL8CD8_summary_table_it.PBMC_day14_TL8 = PBMC_TL8_totCD8_data_it.x_TL8_TotalCD8(day_14_ind_it);
    

    ART1_ind_it = find((PBMC_TL8_totCD8_data_it.DaysP_i_ >=185)&(PBMC_TL8_totCD8_data_it.DaysP_i_ <=195));
    TL8CD8_summary_table_it.PBMC_ART1_day = PBMC_TL8_totCD8_data_it.DaysP_i_(ART1_ind_it);
    TL8CD8_summary_table_it.PBMC_ART1_TL8 = PBMC_TL8_totCD8_data_it.x_TL8_TotalCD8(ART1_ind_it);
    

    ATI1_start_ind_it = find((PBMC_TL8_totCD8_data_it.DaysP_i_ >=215)&(PBMC_TL8_totCD8_data_it.DaysP_i_ <=221));
    TL8CD8_summary_table_it.PBMC_ATI1_start_day = PBMC_TL8_totCD8_data_it.DaysP_i_(ATI1_start_ind_it);
    TL8CD8_summary_table_it.PBMC_ATI1_start_TL8 = PBMC_TL8_totCD8_data_it.x_TL8_TotalCD8(ATI1_start_ind_it);
    

    ATI1_early_ind_it = find((PBMC_TL8_totCD8_data_it.DaysP_i_ >= 225)&(PBMC_TL8_totCD8_data_it.DaysP_i_ <=270));
    TL8CD8_summary_table_it.PBMC_ATI1_early_day = PBMC_TL8_totCD8_data_it.DaysP_i_(ATI1_early_ind_it);
    TL8CD8_summary_table_it.PBMC_ATI1_early_TL8 = PBMC_TL8_totCD8_data_it.x_TL8_TotalCD8(ATI1_early_ind_it);
    

    ATI1_end_ind_it = find((PBMC_TL8_totCD8_data_it.DaysP_i_ >=281)&(PBMC_TL8_totCD8_data_it.DaysP_i_ <= 282));
    TL8CD8_summary_table_it.PBMC_ATI1_end_day = PBMC_TL8_totCD8_data_it.DaysP_i_(ATI1_end_ind_it);
    TL8CD8_summary_table_it.PBMC_ATI1_end_TL8 = PBMC_TL8_totCD8_data_it.x_TL8_TotalCD8(ATI1_end_ind_it);


    ART2_early_ind_it = find((PBMC_TL8_totCD8_data_it.DaysP_i_ >=305)&(PBMC_TL8_totCD8_data_it.DaysP_i_ <= 315));
    TL8CD8_summary_table_it.PBMC_ART2_early_day = PBMC_TL8_totCD8_data_it.DaysP_i_(ART2_early_ind_it);
    TL8CD8_summary_table_it.PBMC_ART2_early_TL8 = PBMC_TL8_totCD8_data_it.x_TL8_TotalCD8(ART2_early_ind_it);


    ART2_end_ind_it = find((PBMC_TL8_totCD8_data_it.DaysP_i_ >=370)&(PBMC_TL8_totCD8_data_it.DaysP_i_ <= 375));
    TL8CD8_summary_table_it.PBMC_ART2_end_day = PBMC_TL8_totCD8_data_it.DaysP_i_(ART2_end_ind_it);
    TL8CD8_summary_table_it.PBMC_ART2_end_TL8 = PBMC_TL8_totCD8_data_it.x_TL8_TotalCD8(ART2_end_ind_it);


    %LNMC data
    LNMC_TL8_totCD8_data_it = LNMC_TL8_totCD8_data(strcmp(LNMC_TL8_totCD8_data.MonkeyID,...
                    VL_data_it.animal_ID),:);

    ART1_ind_it = find((LNMC_TL8_totCD8_data_it.DaysP_i_ >=185)&(LNMC_TL8_totCD8_data_it.DaysP_i_ <=195));
    if (~isempty(ART1_ind_it))
        TL8CD8_summary_table_it.LNMC_ART1_day = LNMC_TL8_totCD8_data_it.DaysP_i_(ART1_ind_it);
        TL8CD8_summary_table_it.LNMC_ART1_TL8 = LNMC_TL8_totCD8_data_it.x_TL8_TotalCD8(ART1_ind_it);
    else
        TL8CD8_summary_table_it.LNMC_ART1_day = nan;
        TL8CD8_summary_table_it.LNMC_ART1_TL8 = nan;
    end


    ATI1_start_ind_it = find((LNMC_TL8_totCD8_data_it.DaysP_i_ >=215)&(LNMC_TL8_totCD8_data_it.DaysP_i_ <=221));
    if (~isempty(ATI1_start_ind_it))
        TL8CD8_summary_table_it.LNMC_ATI1_start_day = LNMC_TL8_totCD8_data_it.DaysP_i_(ATI1_start_ind_it);
        TL8CD8_summary_table_it.LNMC_ATI1_start_TL8 = LNMC_TL8_totCD8_data_it.x_TL8_TotalCD8(ATI1_start_ind_it);
    else

        TL8CD8_summary_table_it.LNMC_ATI1_start_day = nan;
        TL8CD8_summary_table_it.LNMC_ATI1_start_TL8 = nan;
    end


    ATI1_ind_it = find((LNMC_TL8_totCD8_data_it.DaysP_i_ >=222)&(LNMC_TL8_totCD8_data_it.DaysP_i_ <=282));
    if (~isempty(ATI1_ind_it))
        TL8CD8_summary_table_it.LNMC_ATI1_day = LNMC_TL8_totCD8_data_it.DaysP_i_(ATI1_ind_it);
        TL8CD8_summary_table_it.LNMC_ATI1_TL8 = LNMC_TL8_totCD8_data_it.x_TL8_TotalCD8(ATI1_ind_it);
    else

        TL8CD8_summary_table_it.LNMC_ATI1_day = nan;
        TL8CD8_summary_table_it.LNMC_ATI1_TL8 = nan;
    end


    ART2_early_ind_it = find((LNMC_TL8_totCD8_data_it.DaysP_i_ >=305)&(LNMC_TL8_totCD8_data_it.DaysP_i_ <=315));
    if (~isempty(ART2_early_ind_it))
        TL8CD8_summary_table_it.LNMC_ART2_early_day = LNMC_TL8_totCD8_data_it.DaysP_i_(ART2_early_ind_it);
        TL8CD8_summary_table_it.LNMC_ART2_early_TL8 = LNMC_TL8_totCD8_data_it.x_TL8_TotalCD8(ART2_early_ind_it);
    else
        TL8CD8_summary_table_it.LNMC_ART2_early_day = nan;
        TL8CD8_summary_table_it.LNMC_ART2_early_TL8 = nan;
    end
    
    %% adding to table
    if (ii == 1)
        TL8CD8_summary_table = TL8CD8_summary_table_it;

    else
        TL8CD8_summary_table = [TL8CD8_summary_table; TL8CD8_summary_table_it];

    end
    
end

for ii = 1:3
    switch ii
        case 1 %day 14 vs. ART-1, PBMC, all animals
            t1_data = TL8CD8_summary_table.PBMC_day14_TL8;
            t2_data = TL8CD8_summary_table.PBMC_ART1_TL8;
            comparison_name = '%TL8 specific in PBMC on day14 vs. ART1';
            
%         case 2 %ART-1 vs. ART-2 end, PBMC, all animals; #computed in R
%         code
%             t1_data = TL8CD8_summary_table.PBMC_ART1_TL8;
%             t2_data = TL8CD8_summary_table.PBMC_ART2_end_TL8;
%             comparison_name = '%TL8 specific in PBMC at end of ART1 vs. ART2';
            
        case 2 %day 14 vs. ATI-1 end, PBMC, Control animals
            t1_data = TL8CD8_summary_table.PBMC_day14_TL8(strcmp(TL8CD8_summary_table.Treatment, 'Control'));
            t2_data = TL8CD8_summary_table.PBMC_ATI1_end_TL8(strcmp(TL8CD8_summary_table.Treatment, 'Control'));
            comparison_name = '%TL8 specific in PBMC in Controls on day14 vs. end of ATI1';
            
        case 3 %ART-1 vs. ATI-1 end, PBMC, Control animals
            t1_data = TL8CD8_summary_table.PBMC_ART1_TL8(strcmp(TL8CD8_summary_table.Treatment, 'Control'));
            t2_data = TL8CD8_summary_table.PBMC_ATI1_end_TL8(strcmp(TL8CD8_summary_table.Treatment, 'Control'));
            comparison_name = '%TL8 specific in PBMC in Controls at end of ART1 vs. end of ATI1';

    end

    ind_it = find(~isnan(t1_data)&~isnan(t2_data));
    n_it = length(ind_it);

    [p_it, h_it] = signrank(t1_data(ind_it), t2_data(ind_it));

    disp(strcat(comparison_name, ": p-value = ", num2str(p_it), ...
        ", n = ", num2str(n_it), "; paired Wilcoxon signed rank test"))
    
end

%% summary stats
TL8CD8_stats_table_it = table('Size', [8, length(TL8CD8_summary_names)],...
        'VariableTypes', TL8CD8_summary_types,...
        'VariableNames', TL8CD8_summary_names);

TL8CD8_stats_table_it.animal = {'min'; 'max'; 'mean'; 'median'; ...
    'Control min'; 'Control max'; 'Control mean'; 'Control median'};
TL8CD8_stats_table_it{1,3:end} = min(TL8CD8_summary_table{:, 3:end}, [],'omitnan');
TL8CD8_stats_table_it{2,3:end} = max(TL8CD8_summary_table{:, 3:end}, [],'omitnan');
TL8CD8_stats_table_it{3,3:end} = mean(TL8CD8_summary_table{:, 3:end}, 'omitnan');
TL8CD8_stats_table_it{4,3:end} = median(TL8CD8_summary_table{:, 3:end}, 'omitnan');

TL8CD8_stats_table_it{5,3:end} = min(TL8CD8_summary_table{strcmp(TL8CD8_summary_table.Treatment, 'Control'), 3:end}, [],'omitnan');
TL8CD8_stats_table_it{6,3:end} = max(TL8CD8_summary_table{strcmp(TL8CD8_summary_table.Treatment, 'Control'), 3:end}, [],'omitnan');
TL8CD8_stats_table_it{7,3:end} = mean(TL8CD8_summary_table{strcmp(TL8CD8_summary_table.Treatment, 'Control'), 3:end}, 'omitnan');
TL8CD8_stats_table_it{8,3:end} = median(TL8CD8_summary_table{strcmp(TL8CD8_summary_table.Treatment, 'Control'), 3:end}, 'omitnan');

TL8CD8_summary_table = [TL8CD8_summary_table; ...
    TL8CD8_stats_table_it];

save("Summary_statistics.mat", "VL_summary_table", "Protocol_summary_table",...
    "WT_summary_table", "WT_Comparisons_table", "nonWT_vars_geq20percent",...
    "TL8CD8_summary_table");