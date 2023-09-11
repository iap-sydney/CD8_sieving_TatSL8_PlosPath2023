%Steffen Docken
%21-4-21
%This code imports data and saves it as a .csv for use in R
%% AUC calculations are currently done from the time point before VL 
%detected (either in primary infection or during ATIs) and the level of VL
%on the day before detected VL is listed as the limit of detection (60
%copies/ml)

clear
close all
load UPenn_FTY_barcode_SL8_Mutation_data_210607.mat

addpath('General_Functions/');
addpath('Data_importation/');

VL_LOD = 60; %copies/ml. limit of detection for viral load

Num_animals = length(Data_SL8_Plasma);


Seq_Animal_ID_vec = [];
Seq_Animal_Group_vec = [];
Seq_FTY_vec = [];
Sample_vec = [];
WT_count = [];
WT_LOD = [];
WT_det = [];
WT_frac_above_background = [];
Seq_dpi_vec = [];
phase_vec = [];
timing_vec = [];
AUC_VL_post_det_vec = [];
AUC_log10VL_post_det_vec = [];%vectors that will hold final large arrays to be exported
%to R

CD8_depl_Animal_ID_vec = strings(Num_animals, 1);
CD8_depl_vec = zeros(Num_animals,1); %vectors to hold animal name and day
%of CD8 depletion (if it occured)


VL_vec = [];
VL_det_vec = [];
VL_Animal_ID_vec = [];
VL_Animal_Group_vec = [];
VL_FTY_vec = [];
VL_dpi_vec = [];
VL_CD8_depl_vec = [];
%will hold large vectors for VL info exported to R

for ii = 1:Num_animals
    %% get data for Plasma, LNMC, and PBMC overall time course
    [Plasma_data_it, LNMC_data_it, PBMC_data_it] = ...
        SL8_variant_count_210607(Data_SL8_Plasma(ii), Data_SL8_LNMC(ii), ...
        Data_SL8_PBMC(ii));%data for this animals

    VL_data_it = VL_data(ii); %VL data for this animals

    CD8depl_dpi = VL_data_it.CD8depl_dpi;

    CD8_depl_Animal_ID_vec(ii) = VL_data_it.animal_ID;
    if (isempty(CD8depl_dpi))
        CD8_depl_vec(ii) = [];
    else
        CD8_depl_vec(ii) = CD8depl_dpi;
    end
    
    for jj = 1:3 %looping Plasma, LNMC, and PBMC
    
        switch jj
            case 1 %Plasma
                Seq_data_it = Plasma_data_it;
            case 2 %LNMC
                Seq_data_it = LNMC_data_it;
            case 3 %PBMC
                Seq_data_it = PBMC_data_it;
        end
        
        Num_timepoints_it = length(Seq_data_it.dpi);
        
        if ((isempty(Num_timepoints_it))||...
                (~strcmp(VL_data_it.animal_ID,Seq_data_it.animal_ID)))
            continue
        end
        
        
        Seq_Animal_ID_vec_it = strings(Num_timepoints_it, 1);
        Seq_Animal_Group_vec_it = strings(Num_timepoints_it, 1);
        Seq_FTY_vec_it = strings(Num_timepoints_it,1);
        Sample_vec_it = strings(Num_timepoints_it,1);
        WT_count_it = zeros(Num_timepoints_it,1);
        WT_LOD_it = zeros(Num_timepoints_it, 1);
        WT_det_it = zeros(Num_timepoints_it, 1);
        WT_frac_above_background_it = zeros(Num_timepoints_it, 1);
        Seq_dpi_vec_it = Seq_data_it.dpi';%dpi of each sample
        phase_vec_it = strings(Num_timepoints_it,1);%phase of protocol for each sample
        timing_vec_it = Seq_data_it.timing'-1; %timing of within phase of 
        %protocol for each sample (need to subtract 1 from all non-Pre-ART
        %time points to get days since reactivation or ART initiation or 
        %reactivation, since Kevin defines day 1 of ART-1, ATI-1, ART-2,
        %and ATI-2 as day of ART initiation or day of reactivation). Will
        %add 1 back on for Pre-ART time points
        AUC_VL_post_det_vec_it = zeros(Num_timepoints_it, 1);
        AUC_log10VL_post_det_vec_it = zeros(Num_timepoints_it, 1);
        
        
        for kk = 1:Num_timepoints_it
            Seq_Animal_ID_vec_it(kk) = Seq_data_it.animal_ID; %recording animal
            Seq_Animal_Group_vec_it(kk) = VL_data_it.Group; %recording the
            %collection group of the current animal
            Seq_FTY_vec_it(kk) = VL_data_it.FTY_info; %recording if treated with
            %FTY
            
            Sample_vec_it(kk) = Seq_data_it.sample{kk};
            phase_vec_it(kk) = Seq_data_it.phase{kk};%phase of protocol for each sample
            
            WT_count_it(kk) = Seq_data_it.V.count(1,1,kk); %number of sequences
            %with WT
            WT_LOD_it(kk) = Seq_data_it.V.LOD(1,1,kk); %LOD of WT at current
            %time point
            WT_det_it(kk) = Seq_data_it.V.LOD_ind(1,1,kk); %if WT detected (1)
            %or not (0) at current time
            if (WT_det_it(kk) == 1)
                [~, WT_frac_above_background_it(kk)] =...
                    background_subtraction(Seq_data_it.V.count(1,1,kk),...
                    Seq_data_it.V.LOD(1,1,kk), Seq_data_it.S(kk)); %percent WT 
                %above LOD
            else
                WT_frac_above_background_it(kk) = 0; %setting percent WT to 
                %0 if below detection threshold
            end
            
            
            if (jj == 1) %only recording AUCs for Plasma samples
                
                dpi_ind = find(Seq_data_it.dpi(kk) == VL_data_it.t);
            
                if isempty(dpi_ind)
                    disp(strcat("ERROR: sequencing and VL timing do not match for animal ",...
                        VL_data_it.animal_ID, " sequencing at ", ...
                        num2str(Seq_data_it.dpi(kk)), " dpi"));
                    AUC_VL_post_det_vec_it(kk) = NaN;
                    AUC_log10VL_post_det_vec_it(kk) = NaN;
                    continue
                end
                
                %creating vectors of VL measurments and associated timing
                %from either infection or last time before detection (for
                %Pre-ART and ATIs, respectively) up to current sequencing
                %time
                if (strcmp(Seq_data_it.phase{kk}, 'Pre-ART'))
                    timing_vec_it(kk) = timing_vec_it(kk) + 1; %increasing
                    %value for within phase timing during Pre-ART, because
                    %the original timing value was days since infection
                    %with infection day defined as 0 (different from ART
                    %and ATIs)
                    
                    VL_AUC_comp_array = [[0;VL_LOD],...
                        [VL_data_it.t(1:dpi_ind);VL_data_it.VL(1:dpi_ind)]];
                    %array of time and VL level for computing AUC VL and
                    %AUC log10(VL) from infection to current point (use LOD
                    %as VL at first time point)
                else
                    React_dpi = Seq_data_it.dpi(kk) - timing_vec_it(kk);
                    %dpi of reactivation (may be current time if sequenced
                    %first day of reactivation

                    Last_predet_VL_ind = find(VL_data_it.t < React_dpi, 1, 'last');
                    %finding index of last VL measurement before
                    %Reactivation detected
                    
                    VL_AUC_comp_array = [VL_data_it.t(Last_predet_VL_ind:dpi_ind);...
                        VL_data_it.VL(Last_predet_VL_ind:dpi_ind)];
                end
                
                AUC_VL_post_det_vec_it(kk) = ...
                    AUC_exp_between_measure_210507(VL_AUC_comp_array(1,:),...
                    VL_AUC_comp_array(2,:)); %calculating AUC of VL assuming 
                %exponential changes in VL between measurements
                
                AUC_log10VL_post_det_vec_it(kk) = trapz(VL_AUC_comp_array(1,:),...
                    log10(VL_AUC_comp_array(2,:))); %calculating AUC of 
                %log10(VL) when assuming VL changes exponentially between
                %measurements (this results in just a trapezoidal
                %calculation of the integral)
                
            else
                AUC_VL_post_det_vec_it(kk)= NaN;
                AUC_log10VL_post_det_vec_it(kk) = NaN;
            end
                
        end
        
        Seq_Animal_ID_vec = [Seq_Animal_ID_vec; Seq_Animal_ID_vec_it];
        Seq_Animal_Group_vec = [Seq_Animal_Group_vec; Seq_Animal_Group_vec_it];
        Seq_FTY_vec = [Seq_FTY_vec; Seq_FTY_vec_it];
        Sample_vec = [Sample_vec; Sample_vec_it];
        WT_count = [WT_count; WT_count_it];
        WT_LOD = [WT_LOD; WT_LOD_it];
        WT_det = [WT_det; WT_det_it];
        WT_frac_above_background = [WT_frac_above_background; ...
            WT_frac_above_background_it];
        Seq_dpi_vec = [Seq_dpi_vec; Seq_dpi_vec_it];
        phase_vec = [phase_vec; phase_vec_it];
        timing_vec = [timing_vec; timing_vec_it];
        AUC_VL_post_det_vec = [AUC_VL_post_det_vec; AUC_VL_post_det_vec_it];
        AUC_log10VL_post_det_vec = [AUC_log10VL_post_det_vec; ...
            AUC_log10VL_post_det_vec_it];
        %Seq_CD8_depl_vec = [Seq_CD8_depl_vec; Seq_CD8_depl_vec_it];
        
    end
    
    %% recording VL info for VL data frame in R
    VL_vec_it = VL_data_it.VL';
    VL_det_vec_it = VL_data_it.det_ind';
    VL_dpi_vec_it = VL_data_it.t';
    
    num_VL = length(VL_vec_it);
    
    
    VL_Animal_ID_vec_it = strings(num_VL,1);
    VL_Animal_Group_vec_it = strings(num_VL,1);
    VL_FTY_vec_it = strings(num_VL,1);
    %VL_CD8_depl_vec_it = strings(num_VL, 1);
    
    for kk = 1:num_VL
        VL_Animal_ID_vec_it(kk) = VL_data_it.animal_ID;
        VL_Animal_Group_vec_it(kk) = VL_data_it.Group;
        VL_FTY_vec_it(kk) = VL_data_it.FTY_info;
        
    end
    
    VL_vec = [VL_vec; VL_vec_it];
    VL_det_vec = [VL_det_vec; VL_det_vec_it];
    VL_Animal_ID_vec = [VL_Animal_ID_vec; VL_Animal_ID_vec_it];
    VL_Animal_Group_vec = [VL_Animal_Group_vec; VL_Animal_Group_vec_it];
    VL_FTY_vec = [VL_FTY_vec; VL_FTY_vec_it];
    VL_dpi_vec = [VL_dpi_vec; VL_dpi_vec_it];
    %VL_CD8_depl_vec = [VL_CD8_depl_vec; VL_CD8_depl_vec_it];
    
    
end



save('WT_timecourse_data_for_R.mat', 'Group_A_ART_days', 'Group_B_ART_days',...
    'Seq_Animal_ID_vec', 'Seq_Animal_Group_vec', 'Seq_FTY_vec', 'Sample_vec', 'WT_count', ...
    'WT_LOD', 'WT_det', 'WT_frac_above_background', 'Seq_dpi_vec',...
    'phase_vec', 'timing_vec', 'AUC_VL_post_det_vec', 'AUC_log10VL_post_det_vec',...%'Seq_CD8_depl_vec', 
    'VL_vec', 'VL_det_vec', 'VL_Animal_ID_vec', 'VL_Animal_Group_vec',...
    'VL_FTY_vec', 'VL_dpi_vec', 'VL_CD8_depl_vec', 'CD8_depl_Animal_ID_vec',...
    'CD8_depl_vec');

%% saving data in spreadsheets
num_WT_data_points = length(WT_frac_above_background);

WT_frac_data_array = cell(num_WT_data_points+1, 13);

WT_frac_data_array{1,1} = "WT Fraction above Background";
WT_frac_data_array{1,2} = "WT detection";
WT_frac_data_array{1,3} = "WT count";
WT_frac_data_array{1,4} = "WT LOD";
WT_frac_data_array{1,5} = "Animal";
WT_frac_data_array{1,6} = "FTY Treatment";
WT_frac_data_array{1,7} = "Sample";
WT_frac_data_array{1,8} = "dpi";
WT_frac_data_array{1,9} = "phase";
WT_frac_data_array{1,10} = "timing";
WT_frac_data_array{1,11} = "AUC VL post detection";
WT_frac_data_array{1,12} = "AUC log10VL post detection";
WT_frac_data_array{1,13} = "Collection Group";

for ii = 1:num_WT_data_points
    WT_frac_data_array{ii+1,1} = WT_frac_above_background(ii);
    WT_frac_data_array{ii+1,2} = WT_det(ii);
    WT_frac_data_array{ii+1,3} = WT_count(ii);
    WT_frac_data_array{ii+1,4} = WT_LOD(ii);
    WT_frac_data_array{ii+1,5} = Seq_Animal_ID_vec(ii);
    WT_frac_data_array{ii+1,6} = Seq_FTY_vec(ii);
    WT_frac_data_array{ii+1,7} = Sample_vec(ii);
    WT_frac_data_array{ii+1,8} = Seq_dpi_vec(ii);
    WT_frac_data_array{ii+1,9} = phase_vec(ii);
    WT_frac_data_array{ii+1,10} = timing_vec(ii);
    WT_frac_data_array{ii+1,11} = AUC_VL_post_det_vec(ii);
    WT_frac_data_array{ii+1,12} = AUC_log10VL_post_det_vec(ii);
    WT_frac_data_array{ii+1,13} = Seq_Animal_Group_vec(ii);
    
end

writecell(WT_frac_data_array, 'WT_timecourse_data_for_R.csv');


num_VL_data_points = length(VL_vec);

VL_data_array = cell(num_VL_data_points+ 1, 6);

VL_data_array{1,1} = "Viral Load (copies/ml)";
VL_data_array{1,2} = "Viral Load detection";
VL_data_array{1,3} = "Animal";
VL_data_array{1,4} = "FTY Treatment";
VL_data_array{1,5} = "dpi";
%VL_data_array{1,6} = "CD8 depletion";
VL_data_array{1,6} = "Collection Group";

for ii = 1:num_VL_data_points
    VL_data_array{ii+1,1} = VL_vec(ii);
    VL_data_array{ii+1,2} = VL_det_vec(ii);
    VL_data_array{ii+1,3} = VL_Animal_ID_vec(ii);
    VL_data_array{ii+1,4} = VL_FTY_vec(ii);
    VL_data_array{ii+1,5} = VL_dpi_vec(ii);
    %VL_data_array{ii+1,6} = VL_CD8_depl_vec(ii);
    VL_data_array{ii+1,6} = VL_Animal_Group_vec(ii);
end


writecell(VL_data_array, 'VL_timecourse_data_for_R.csv');

ART_days_array = cell(3, 5);
ART_days_array{1,2} = "First Day of ART1";
ART_days_array{1,3} = "First Day without ART1";
ART_days_array{1,4} = "First Day of ART2";
ART_days_array{1,5} = "First Day without ART2";
ART_days_array{2,1} = "Group A";
ART_days_array{2,2} = Group_A_ART_days(1);
ART_days_array{2,3} = Group_A_ART_days(2);
ART_days_array{2,4} = Group_A_ART_days(3);
ART_days_array{2,5} = Group_A_ART_days(4);
ART_days_array{3,1} = "Group B";
ART_days_array{3,2} = Group_B_ART_days(1);
ART_days_array{3,3} = Group_B_ART_days(2);
ART_days_array{3,4} = Group_B_ART_days(3);
ART_days_array{3,5} = Group_B_ART_days(4);

writecell(ART_days_array, 'ART_days_data_for_R.csv');


CD8_depl_days_array = cell(Num_animals+1, 2);
CD8_depl_days_array{1,1} = "Animal";
CD8_depl_days_array{1,2} = "CD8depl_dpi";
for ii = 1:Num_animals
    CD8_depl_days_array{ii+1, 1} = CD8_depl_Animal_ID_vec(ii);
    CD8_depl_days_array{ii + 1, 2} = CD8_depl_vec(ii);
end

writecell(CD8_depl_days_array, 'CD8depl_days_data_for_R.csv');
