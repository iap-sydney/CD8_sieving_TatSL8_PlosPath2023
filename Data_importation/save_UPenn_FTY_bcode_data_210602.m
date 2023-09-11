%Steffen Docken
%09-04-21 (Original code)
%This code imports and saves the FTY barcode and SL8 mutation data sent on
%07-06-21

%_raw means including all sequencing runs
%_short, and _SL8 means only including sequencing runs above Template and
%amplification cut off, and when there are repeated runs above cut offs,
%only the one with the most templates is used, but only looking at _short
%or _SL8 sequencing

clear

Template_Plasma_cutoff = 1e3;
Template_DNA_cutoff = 200; %at least 1000 templates of Plasma and 200 of 
%DNA to consider data
Sequence_Plasma_cutoff = 200;
Sequence_DNA_cutoff = 100; %at least 200 sequences of Plasma and 100 of 
%DNA to consider data

cutoff_vec = [Template_Plasma_cutoff, Template_DNA_cutoff,...
    Sequence_Plasma_cutoff, Sequence_DNA_cutoff];

UPenn_FTY_barcode_data = '../Data/Penn SIV barcode data summary.SD_20230208.xlsx';

data_sheets = sheetnames(UPenn_FTY_barcode_data);

animal_list = data_sheets(2:end);

Num_animals = length(animal_list);

%initializing data structure for sequencing output
Data_struct_init.T = []; %number of templates at each time point
Data_struct_init.S = []; %number of sequence reads at each time point
Data_struct_init.A = []; %amplification at each time point
Data_struct_init.dpi = []; %day post infection of each time point
Data_struct_init.phase = {}; %protocol phase (e.g., ART-1, ATI-1, etc.) of
% each time point
Data_struct_init.timing = []; %timing within phase of each time point
Data_struct_init.sample = {}; %sample sequenced (Plasma, PBMC, or LNMC) at
% each time point
Data_struct_init.seq = {}; %short or SL8 sequencing done at each time point
Data_struct_init.barcode_ID = []; %list of all barcode IDs in sequencing 
% output (same order as in count data)
Data_struct_init.SL8_ID = []; %list of IDs for all SL8 variants in 
% sequencing output (same order as in count data)

%barcode count data (3-D arrays: each row a different barcode, 1 column,
%each element of 3rd dimension a different sequencing run/time point
B_struct_init.count = []; %count
B_struct_init.LOD = []; %limit of detection
B_struct_init.LOD_ind = []; %indicates if corresponding barcode is above 
% the LOD at the corresponding time point

%SL8 nucleotide variant count data (3-D arrays: 1 row, each column a 
% different variant, each element of 3rd dimension a different sequencing 
% run/time point
M_struct_init.count = [];%count
M_struct_init.LOD = [];%limit of detection
M_struct_init.LOD_ind = []; %indicates if corresponding variant is above 
% the LOD at the corresponding time point

%barcode-SL8 nucleotide variant combination count data (3-D arrays: each 
% row a different barcode, each column a different variant, each element of
% 3rd dimension a different sequencing run/time point
W_struct_init.count = [];%count
W_struct_init.LOD =[];%limit of detection
W_struct_init.LOD_ind = []; %indicates if corresponding barcode-variant
% combination is above the LOD at the corresponding time point

Data_struct_init.B = B_struct_init;
Data_struct_init.M = M_struct_init;
Data_struct_init.W = W_struct_init;

Data_struct_init.animal_ID = []; %ID for animal

Data_short_Plasma(Num_animals) = Data_struct_init;
Data_short_LNMC(Num_animals) = Data_struct_init;
Data_short_PBMC(Num_animals) = Data_struct_init;
Data_SL8_Plasma(Num_animals) = Data_struct_init;
Data_SL8_LNMC(Num_animals) = Data_struct_init;
Data_SL8_PBMC(Num_animals) = Data_struct_init;
Data_all_runs(Num_animals) = Data_struct_init;
%initializing data structures that will hold all data runs (short, SL8, 
% %and runs with bad template number or low amplification)

for ii = 1:Num_animals
    Mutation_data_headings_raw = readcell(UPenn_FTY_barcode_data,...
        'Sheet', animal_list{ii}, 'Range', '3:7'); %collecting the column
    %heading info.  Type of cells or virus (Plasma, LN, or PBMC) is in the
    %first row, type of sequencing is the first item in the 4th row, and day
    %post inoculation is at the end of the last row
    
    Num_seq_runs_raw = size(Mutation_data_headings_raw,2)-1; %number of sequencing
    %runs for this animal (first column is barcode/mutation info)
    
    %% recording timing of sequencing, type, #templates, and #sequences
    T_raw = zeros(1, Num_seq_runs_raw); %number of templates for each data time point
    dpi_raw = zeros(1, Num_seq_runs_raw); %days post infection of each data time point
    phase_raw = cell(1, Num_seq_runs_raw); %phase of protocol of each time point 
    %(e.g., Pre-ART, ART-1, ATI-1, etc.)
    timing_raw = zeros(1, Num_seq_runs_raw); %timing within phase of protocol for 
    %each time point (e.g., for ATI-1.007, it is 7, because it is the 7th
    %day since first detection of VL in ATI-1; first day of detection during 
    % ATI is denoted 1)
    sample_raw = cell(1, Num_seq_runs_raw); %Sample type (PBMC, LNMC, Plasma)
    seq_raw = cell(1, Num_seq_runs_raw); %Sequence type (short or SL8)
    
    for jj = 1:Num_seq_runs_raw
        
        if (~isnumeric(Mutation_data_headings_raw{2, jj +1}))
            disp(strcat("ERROR: Templates not numeric in animal ",...
                animal_list{ii}, ", entry ", num2str(jj)));
            continue
        end
        
        T_raw(jj) = Mutation_data_headings_raw{2, jj +1}; %number of templates
        
        dpi_info1 = split(Mutation_data_headings_raw{end, jj + 1}, 'day ');
        dpi_info2 = split(dpi_info1{2}, ' (');
        
        dpi_raw(jj) = str2double(dpi_info2(1)); %dpi
        
        phase_info1 = split(Mutation_data_headings_raw{end, jj + 1}, " ");
        phase_info2 = split(phase_info1{1}, '.');
        phase_raw{jj} = phase_info2{1};%phase of protocol
        timing_raw(jj) = str2double(phase_info2(2));%timing within phase
        
        sample_raw{jj} = Mutation_data_headings_raw{1,jj + 1}; %PBMC, LNMC, or Plasma
        
        seq_type = split(Mutation_data_headings_raw{4,jj+1}, '; ');
        seq_raw{jj} = seq_type{1}; %short or SL8
        
    end
    
    %% importing barcode-mutation names and raw count data
    Barcode_mut_names_raw_raw = readcell(UPenn_FTY_barcode_data, ...
        'Sheet', animal_list{ii}, 'Range', 'A:A'); %barcode and mutation names
    %start on row 9.
    
    Barcode_mut_names_raw_raw = Barcode_mut_names_raw_raw(9:end); %removing 
    %initial rows
    
    Seq_count_raw_raw = readmatrix(UPenn_FTY_barcode_data, ...
        'Sheet', animal_list{ii}, 'Range', 'B9'); %raw sequence numbers
    
    %% getting list of all barcodes and SL8 IDs ever seen in this animal
    %including ones that may be below LOD
    [barcode_ID_raw, SL8_ID_raw] = bcode_SL8ID_list_identification_210604(Barcode_mut_names_raw_raw, animal_list{ii});
    
    num_bcode_raw = length(barcode_ID_raw); %number of barcodes detected in this
    %animal
    num_var_raw = length(SL8_ID_raw); %number of SL8 variants detected in this
    %animal
    
    %% recording counts of each barcode, each SL8 variant, and each barcode
    %-SL8 variant combo
    
    [B_raw, M_raw, W_raw, S_raw, A_raw] = bcode_SL8ID_total_count_210604(Barcode_mut_names_raw_raw,...
        Seq_count_raw_raw, barcode_ID_raw, SL8_ID_raw, T_raw, animal_list{ii});
    %B_raw struct that will hold count of each 
    %barcode at each time point (rows and 3rd Dimension), and LOD, and
    %indicator of if barcode is above LOD (0 = below LOD, 1 = above LOD)
    
    %M_raw struct that will hold count of each 
    %SL8 variant at each time point (columns and 3rd Dimension), and LOD, and
    %indicator of if barcode is above LOD (0 = below LOD, 1 = above LOD)
    
    %W_raw struct that will hold count of each barcode with each
    %SL8 variant at each time point (rows, columns, and 3rd Dimension), and LOD, and
    %indicator of if barcode is above LOD (0 = below LOD, 1 = above LOD)
    
    %S_raw is number of sequences for each data time point with barcode and
    %SL8 variant identified
    
    %A_raw is amplification factor for each data time point (S/T)
    
    
    
    %% Seperate out Short, Good SL8, and old runs (and seperate by sample type)
    
    [final_short_Plasma_ind, final_short_LNMC_ind, final_short_PBMC_ind,...
    final_SL8_Plasma_ind, final_SL8_LNMC_ind, final_SL8_PBMC_ind] = Best_seq_run_indices_210604(dpi_raw,...
        sample_raw, seq_raw, T_raw, S_raw, cutoff_vec); %obtaining indices for best 
    %sequencing runs (most templates) of short and SL8 sequencing runs of
    %each sample type (also ensuring enough templates and sequences
    %were attained
    
    for jj = 1:7
        switch jj
            case 1 %short Plasma
                ind_list = final_short_Plasma_ind;
            case 2 %short LNMC
                ind_list = final_short_LNMC_ind;
            case 3 %short PBMC
                ind_list = final_short_PBMC_ind;
            case 4 %SL8 Plasma
                ind_list = final_SL8_Plasma_ind;
            case 5 %SL8 LNMC
                ind_list = final_SL8_LNMC_ind;
            case 6 %SL8 PBMC
                ind_list = final_SL8_PBMC_ind;
            case 7 %all of raw data
                ind_list = 1:Num_seq_runs_raw;
        end
        
        if (isempty(ind_list)) %skipping iteration if not data for this 
            %sequencing type and sample type combo
            continue
        end
        
        Data_it.T = T_raw(ind_list);
        Data_it.S = S_raw(ind_list);
        Data_it.A = A_raw(ind_list);
        Data_it.dpi = dpi_raw(ind_list);
        Data_it.phase = phase_raw(ind_list);
        Data_it.timing = timing_raw(ind_list);
        Data_it.sample = sample_raw(ind_list);
        Data_it.seq = seq_raw(ind_list);
        Data_it.barcode_ID = barcode_ID_raw;
        Data_it.SL8_ID = SL8_ID_raw;
        
        B_it.count = B_raw.count(:,1, ind_list);
        B_it.LOD = B_raw.LOD(:,1, ind_list);
        B_it.LOD_ind = B_raw.LOD_ind(:, 1, ind_list);
        
        M_it.count = M_raw.count(1, :, ind_list);
        M_it.LOD = M_raw.LOD(1, :, ind_list);
        M_it.LOD_ind = M_raw.LOD_ind(1, :, ind_list);
        
        W_it.count = W_raw.count(:,:, ind_list);
        W_it.LOD = W_raw.LOD(:,:, ind_list);
        W_it.LOD_ind = W_raw.LOD_ind(:, :, ind_list);
        
        Data_it.B = B_it;
        Data_it.M = M_it;
        Data_it.W = W_it;
        
        Data_it.animal_ID = animal_list{ii};
        
        
        switch jj
            case 1 %short Plasma
                Data_short_Plasma(ii) = Data_it;
            case 2 %short LNMC
                Data_short_LNMC(ii) = Data_it;
            case 3 %short PBMC
                Data_short_PBMC(ii) = Data_it;
            case 4 %SL8 Plasma
                Data_SL8_Plasma(ii) = Data_it;
            case 5 %SL8 LNMC
                Data_SL8_LNMC(ii) = Data_it;
            case 6 %SL8 PBMC
                Data_SL8_PBMC(ii) = Data_it;
            case 7 %all of raw data
                Data_all_runs(ii) = Data_it;
        end
    end
            
    %% reordering short and SL8 sequencing data by size of barcodes
    [Data_short_Plasma(ii), Data_short_LNMC(ii), Data_short_PBMC(ii)] = ...
        sort_bcode_210607(Data_short_Plasma(ii), Data_short_LNMC(ii), Data_short_PBMC(ii));
    
    [Data_SL8_Plasma(ii), Data_SL8_LNMC(ii), Data_SL8_PBMC(ii)] = ...
        sort_bcode_210607(Data_SL8_Plasma(ii), Data_SL8_LNMC(ii), Data_SL8_PBMC(ii));
     
    
    animal_list{ii} %to track as animal data finishes being processed
end
   
    
%% import VL info
%% VL data

UPenn_FTY_VL_data = '../Data/Copy of 1st Betts FTY study 1 VL_day106 post MT807R1_SSD_necropsy_dpi_added_dates_corrected_230207.xlsx';

VL_thresh_val = 60; %value to be used for data points below threshold

data_sheets = sheetnames(UPenn_FTY_VL_data);

VL_animals_GroupA = readcell(UPenn_FTY_VL_data, 'Sheet', data_sheets{2},...
    'Range', 'K4:Q4');
num_GroupA = length(VL_animals_GroupA);
VL_animals_GroupB = readcell(UPenn_FTY_VL_data, 'Sheet', data_sheets{3},...
    'Range', 'I4:O4'); %saving the animals that are in Collection Group A 
%and B and the order they are in.
num_GroupB = length(VL_animals_GroupB);

VL_animals_FTY_GroupA = readcell(UPenn_FTY_VL_data, 'Sheet', data_sheets{2},...
    'Range', 'K5:Q5');
VL_animals_FTY_GroupB = readcell(UPenn_FTY_VL_data, 'Sheet', data_sheets{3},...
    'Range', 'I5:O5'); %saving the info on if animals were treated with FTY
%or control

VL_data_GroupA_num = readmatrix(UPenn_FTY_VL_data, 'Sheet', data_sheets{2},...
    'Range', 'F7:Q71');
VL_data_GroupB_num = readmatrix(UPenn_FTY_VL_data, 'Sheet', data_sheets{3},...
    'Range', 'D7:O64');
VL_data_GroupA_txt = readcell(UPenn_FTY_VL_data, 'Sheet', data_sheets{2},...
    'Range', 'F7:Q71');
VL_data_GroupB_txt = readcell(UPenn_FTY_VL_data, 'Sheet', data_sheets{3},...
    'Range', 'D7:O64');

%initializing data structure for viral load data
VL_data_struct_init.animal_ID =[]; %will hold animal ID
VL_data_struct_init.t = [];%will hold list of dpi for VL of animal
VL_data_struct_init.VL = [];%will hold measured VL values
VL_data_struct_init.det_ind = [];%will hold list of indices indicating if 
%values in .VL are above LOD (1) or below (0);
VL_data_struct_init.Group = [];%will hold 'A' or 'B' depending on the 
%collection group
VL_data_struct_init.FTY_info = []; %will hold if animal was treated with 
%FTY (FTY) or not (Control);
VL_data_struct_init.nec_dpi = []; %will hold the day post infection the 
%animal was necropsied if it was necropsied earlier than the other animals
VL_data_struct_init.CD8depl_dpi = []; %will hold the day post infection the
%animal was CD8 depleted if it was CD8 depleted.

VL_data(Num_animals) = VL_data_struct_init; 

%dates 4 animals were necropsied (from "Copy of 1st Betts FTY study 1
%VL_day106 post MT807R1_SSD_necropsy_dpi_added_210426.xlsx")
nec_14_10 = 392;
nec_13D064 = 396;
nec_RJy13 = 396;
nec_REe16 = 403;

%dates of CD8 depletion (from "Copy of 1st Betts FTY study 1
%VL_day106 post MT807R1_SSD_necropsy_dpi_added_210426.xlsx")
Collect_B_CD8depl = 465; %collection group B depleted on day 465 post infection
Collect_A1_CD8depl = 468; %CD8 depleted for all Collection group A animals
%except RYk16 and RYm15
Collect_A_RYk16_RYm15_CD8depl = 471; %CD8 depleted for animals RYk16 and RYm15

for ii = 1:Num_animals
    VL_data_it = VL_data_struct_init; %initializing VL_data_it
    
    VL_data_it.animal_ID = animal_list{ii};
    for jj = 1:num_GroupA
        if strcmp(VL_data_it.animal_ID, VL_animals_GroupA{jj})
            VL_data_num_it = VL_data_GroupA_num;
            VL_data_txt_it = VL_data_GroupA_txt;
            Group_ind = jj; %index for where this animal is in Group A
            VL_data_it.Group = 'A';
            if strcmp(VL_animals_FTY_GroupA{jj}, 'FTY')
                VL_data_it.FTY_info = 'FTY';
            elseif strcmp(VL_animals_FTY_GroupA{jj}, 'Control')
                VL_data_it.FTY_info = 'Control';
            else
                disp('ERROR: FTY status not recorded');
            end
            break; %ending loop searching for what collection group this 
            %animal is in
        elseif strcmp(VL_data_it.animal_ID, VL_animals_GroupB{jj})
            VL_data_num_it = VL_data_GroupB_num;
            VL_data_txt_it = VL_data_GroupB_txt;
            Group_ind = jj; %index for where this animal is in Group B
            VL_data_it.Group = 'B';
            if strcmp(VL_animals_FTY_GroupB{jj}, 'FTY')
                VL_data_it.FTY_info = 'FTY';
            elseif strcmp(VL_animals_FTY_GroupB{jj}, 'Control')
                VL_data_it.FTY_info = 'Control';
            else
                disp('ERROR: FTY status not recorded');
            end
            break; %ending loop searching for what collection group this 
            %animal is in
        end
        
        if (jj == num_GroupA)
            disp(strcat("Animal ", VL_data_it.animal_ID, " not found in VL data")); 
            %this means code got to 
            %end of loop and didn't find the group for this animal
        end
        
    end
    
    for jj = 1:length(VL_data_num_it(:,1))
        if ~isfinite(VL_data_num_it(jj, 2))
            continue %skipping this row because no day post infection is listed
        end
        
        if isfinite(VL_data_num_it(jj, 5+Group_ind))
            VL_data_it.t = [VL_data_it.t, VL_data_num_it(jj, 2)];
            VL_data_it.VL = [VL_data_it.VL, VL_data_num_it(jj, 5+Group_ind)];
            VL_data_it.det_ind = [VL_data_it.det_ind, 1]; %entering that VL
            % detected
        elseif strcmp(VL_data_txt_it{jj, 5+Group_ind}, '<60')
            VL_data_it.t = [VL_data_it.t, VL_data_num_it(jj, 2)];
            VL_data_it.VL = [VL_data_it.VL, VL_thresh_val]; %entering LOD for
            %points below LOW
            VL_data_it.det_ind = [VL_data_it.det_ind, 0]; %entering that VL
            %not detected
        elseif (strcmp(VL_data_txt_it(jj, 5+Group_ind), 'n/a')||...
                ismissing(VL_data_txt_it{jj, 5+Group_ind}))
            continue; %not entering an empty row in the data
        else
            disp('VL data value not entered');
        end
    end
    
    %% recording day of necropsy or CD8 depletion
    if (strcmp(VL_data_it.animal_ID, '14_10'))
        VL_data_it.nec_dpi = nec_14_10;
    elseif (strcmp(VL_data_it.animal_ID, '13D064'))
        VL_data_it.nec_dpi = nec_13D064;
    elseif (strcmp(VL_data_it.animal_ID, 'RJy13'))
        VL_data_it.nec_dpi = nec_RJy13;
    elseif (strcmp(VL_data_it.animal_ID, 'REe16'))
        VL_data_it.nec_dpi = nec_REe16;
    elseif (strcmp(VL_data_it.animal_ID, 'RYk16')||...
            strcmp(VL_data_it.animal_ID, 'RYm15'))
        VL_data_it.CD8depl_dpi = Collect_A_RYk16_RYm15_CD8depl;
    elseif (strcmp(VL_data_it.Group, 'A'))
        VL_data_it.CD8depl_dpi = Collect_A1_CD8depl;
    elseif (strcmp(VL_data_it.Group, 'B'))
        VL_data_it.CD8depl_dpi = Collect_B_CD8depl;
    end
    
    VL_data(ii) = VL_data_it;%saving current VL data struct for current animal
end


Group_A_ART_days = zeros(4,1);
Group_B_ART_days = zeros(4,1); %will hold day post infection of ART 
%initiation, 1st ART cessation, ART reinitiation, and 2nd ART cessation

Group_A_ART_days(1) = 14; %first day of ART
Group_A_ART_days(2) = 220; %first day without ART after first phase
Group_A_ART_days(3) = 281; %first day of second phase of ART
Group_A_ART_days(4) = 374; %first day without ART after second phase

Group_B_ART_days(1) = 14; %first day of ART
Group_B_ART_days(2) = 219; %first day without ART after first phase
Group_B_ART_days(3) = 282; %first day of second phase of ART
Group_B_ART_days(4) = 375; %first day without ART after second phase
    
save('../UPenn_FTY_barcode_SL8_Mutation_data_210607.mat', 'Data_all_runs',...
    'Data_short_Plasma', 'Data_short_LNMC', 'Data_short_PBMC',...
    'Data_SL8_Plasma', 'Data_SL8_LNMC', 'Data_SL8_PBMC',...
    'Group_A_ART_days', 'Group_B_ART_days', 'VL_data');