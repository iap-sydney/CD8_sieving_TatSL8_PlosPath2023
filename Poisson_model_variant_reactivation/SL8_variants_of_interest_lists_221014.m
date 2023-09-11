%Steffen Docken
%8-11-21
%Code to generate the lists of aa mutations that are 1 nt change from WT
%and common variants from the literature

clear
close all

load ../UPenn_FTY_barcode_SL8_Mutation_data_210607.mat

addpath('../Data_importation/')

Num_animals = length(Data_SL8_Plasma);

Single_aa_mut_list = []; %will contain the list of all SL8 aa mutations 
% that are (or can be) a single nt mutation from WT. SL8 variants listed by 
% their 1-D SL8 ID

Single_aa_mut_info_array = []; %will contain the list of all SL8 aa mutations 
% SL8 variants listed by 
% their 1-D SL8 ID, aa sequence, and nt sequence

for ii = 1:Num_animals
    for jj = 1:3
        switch jj
            case 1
                Data_it = Data_SL8_Plasma;
            case 2
                Data_it = Data_SL8_LNMC;
            case 3
                Data_it = Data_SL8_PBMC;
        end
        

        if (isempty(Data_it(ii).SL8_ID)) %skipping animals-sample combinations
            %with no data
            continue
        end

        aa_mut_ind_array = Data_it(ii).SL8_ID(2:9, :) > 0; %getting 
        % indicator for if each aa (rows) is mutated in each SL8 variant
        % (columns)

        num_aa_mut_it = sum(aa_mut_ind_array, 1); %getting total number of
        % aa mutations in each SL8 variant

        single_aa_mut_ind_it = find(num_aa_mut_it ==1);

        Single_aa_mut_list = [Single_aa_mut_list, Data_it(ii).SL8_ID(1, single_aa_mut_ind_it)];
        %adding current list of SL8 variants with single aa mutations

        Single_aa_mut_info_array = [Single_aa_mut_info_array, ...
            Data_it(ii).SL8_ID(:, single_aa_mut_ind_it)];

    end
end

Single_aa_mut_list = unique(Single_aa_mut_list); %removing duplicates in 
% list

%% analysis of single amino acid mutations

WT_SL8_aa_seq = 'STPESANL'; %WT amino acid sequence for Tat-SL8
WT_SL8_nt_seq = 'TCCACTCCAGAATCGGCCAACCTG'; %WT nucleotide sequence for 
% Tat-SL8

num_single_aa_mut = length(Single_aa_mut_list);

Single_aa_mut_info_array = unique(Single_aa_mut_info_array', "rows")';
    %removing repeats from multiple animals

Single_aa_mut_struct_init.SL8_ID = 0; %single numeric code for SL8 variant
Single_aa_mut_struct_init.aa_mut = {}; %will contain the amino acid sequence 
% for the SL8 variant
Single_aa_mut_struct_init.aa_mut_numeric = zeros(8, 1); %8-D vector numeric 
% code for amino acid SL8 variant
Single_aa_mut_struct_init.num_nt_mut_var = 0; %number of corresponding
% nucleotide variants
Single_aa_mut_struct_init.nt_mut = {}; %will contain all nt sequences
% corresponding to the current SL8 variant
Single_aa_mut_struct_init.nt_mut_numeric = zeros(8, 1); %8-D vector numeric 
% code for nucleotide SL8 variants
Single_aa_mut_struct_init.nt_mut_path_dist = 0; %number of nucleotide 
% mutations for each corresponding nt sequence
Single_aa_mut_struct_init.min_path_dist = 0; %minimum number of mutations
%across all nt variants for this amino acid sequence
Single_aa_mut_struct_init.stop_codon = 0; %did mutations add a stop codon? 
%(1 = yes, 0 = no)
Single_aa_mut_struct_init.min_path_count = 0; %number of nt sequences that 
% contain the minimum number of mutations
Single_aa_mut_struct_init.min_path_nt_mut = {}; %will contain all nt sequences
% corresponding to the current SL8 variant that are of the minimum length
Single_aa_mut_struct_init.min_path_nt_mut_numeric = zeros(8,1);% will 
% contain the 8-D vector numeric code for all nucleotide SL8 variants of
% the minimum length
Single_aa_mut_struct_init.min_path_nt_mut_contained = 0; %for each nt 
% sequence, this will contain the index of the corresponding 
% min_path_nt_mut(_numeric) variant that contains the minimal mutation
% present in the given nt sequence
Single_aa_mut_struct_init.num_w_min_path = 0; %will contain the number of nt
% variants that contain each minimal mutation listed in
% min_path_nt_mut(_numeric)

analysis_aa_mut_info_init.SL8ID = 0; %single numeric code for SL8 variant
analysis_aa_mut_info_init.SL8aaID = zeros(8,1);%8-D vector numeric 
% code for amino acid SL8 variant
analysis_aa_mut_info_init.orig_SL8ntID_set = zeros(8,1);%8-D vector numeric 
% code for nucleotide SL8 variants

for ii = 1:num_single_aa_mut
    aa_mut_SL8_it = Single_aa_mut_list(ii);

    Single_aa_mut_info_array_it = ...
        Single_aa_mut_info_array(:, Single_aa_mut_info_array(1,:) == aa_mut_SL8_it);

    Single_aa_mut_struct_it = Single_aa_mut_struct_init;

    Single_aa_mut_struct_it.SL8_ID = aa_mut_SL8_it; %SL8 variant ID
    Single_aa_mut_struct_it.aa_mut_numeric = Single_aa_mut_info_array_it(2:9,1); 
    %SL8 variant aa sequence

    SL8_aa_seq_it = WT_SL8_aa_seq; %redefining current SL8 aa sequence to WT 
    %(which will be edited)
    for kk = 1:8 %obtaining the aa sequence for the current variant
        if(Single_aa_mut_struct_it.aa_mut_numeric(kk) ~= 0)
            aa_it = aa_translation_210603(Single_aa_mut_struct_it.aa_mut_numeric(kk)); 
            if (strcmp(aa_it, 'Stop'))
                SL8_aa_seq_it(kk) = '.';
                Single_aa_mut_struct_it.stop_codon = 1; %indicating the mutation
                %gives a stop codon
            else
                SL8_aa_seq_it(kk) = aa_it;
            end
        end
    end

    Single_aa_mut_struct_it.aa_mut = SL8_aa_seq_it; %recording the aa sequence
    %for the current variant

    %nt sequence
    num_nt_mut_var_it = length(Single_aa_mut_info_array_it(1,:));
    Single_aa_mut_struct_it.num_nt_mut_var = num_nt_mut_var_it; %number of nt
    %variants of current aa variant

    Single_aa_mut_struct_it.nt_mut_numeric = Single_aa_mut_info_array_it(10:17,:);
    %SL8 variant nt sequence

    Single_aa_mut_struct_it.nt_mut = cell(1,num_nt_mut_var_it);
    Single_aa_mut_struct_it.nt_mut_path_dist = zeros(1, num_nt_mut_var_it);

    for jj = 1:num_nt_mut_var_it
        SL8_nt_seq_it = WT_SL8_nt_seq; %redefining current SL8 nt sequence to WT 
        %(which will be edited)
        
        nt_mut_numeric_it = Single_aa_mut_struct_it.nt_mut_numeric(:,jj);
        for kk = 1:8%obtaining the nt sequence for the current variant
            if(nt_mut_numeric_it(kk) ~= 0)
                SL8_nt_seq_it((3*(kk-1)+1):(3*kk)) = ...
                    nt_translation_210603(nt_mut_numeric_it(kk));
            end
        end
        
        Single_aa_mut_struct_it.nt_mut{jj} = SL8_nt_seq_it;

        nt_mut_ind_it = (WT_SL8_nt_seq ~= SL8_nt_seq_it); %logical array 
        %indicating which nt were mutated (1) or not (0)

        Single_aa_mut_struct_it.nt_mut_path_dist(jj) = sum(nt_mut_ind_it);
        %keyboard
    end
    
    %minimum number of nt mutations:
    Single_aa_mut_struct_it.min_path_dist = min(Single_aa_mut_struct_it.nt_mut_path_dist);
    if (Single_aa_mut_struct_it.min_path_dist < 1)
        disp('ERROR: less than 1 mutation in variant')
    end

    min_path_var_ind_it = find(Single_aa_mut_struct_it.nt_mut_path_dist == ...
        Single_aa_mut_struct_it.min_path_dist);

    Single_aa_mut_struct_it.min_path_count = length(min_path_var_ind_it); 
    %number of nt variants with minimum path distance

    Single_aa_mut_struct_it.min_path_nt_mut = cell(1,Single_aa_mut_struct_it.min_path_count);
    Single_aa_mut_struct_it.min_path_nt_mut_numeric = zeros(8,Single_aa_mut_struct_it.min_path_count);
    Single_aa_mut_struct_it.num_w_min_path = zeros(1, Single_aa_mut_struct_it.min_path_count);

    Single_aa_mut_struct_it.min_path_nt_mut_contained = ...
        zeros(1, num_nt_mut_var_it); %will contain the index of the minimum 
    %non-synonymous mutation contained in this nt variant

    for jj = 1:Single_aa_mut_struct_it.min_path_count
        Single_aa_mut_struct_it.min_path_nt_mut{jj} = ...
            Single_aa_mut_struct_it.nt_mut{min_path_var_ind_it(jj)};

        Single_aa_mut_struct_it.min_path_nt_mut_numeric(:,jj) = ...
            Single_aa_mut_struct_it.nt_mut_numeric(:, min_path_var_ind_it(jj));
        % recording the nt sequences for the variants with the minimum path 
        % distance. recorded as alphabetic sequence and numeric code for 
        % codons

        nt_mut_ind_vec_it = find(WT_SL8_nt_seq ~= ...
            Single_aa_mut_struct_it.min_path_nt_mut{jj}); %logical array 
        %indicating which nt were mutated (1) or not (0)
        
        for kk = 1:num_nt_mut_var_it
            if (strcmp(Single_aa_mut_struct_it.min_path_nt_mut{jj}(nt_mut_ind_vec_it),...
                Single_aa_mut_struct_it.nt_mut{kk}(nt_mut_ind_vec_it)))
                %comparing if the mutations in the current minimum path 
                % variant are present in the kk^th variant

                Single_aa_mut_struct_it.num_w_min_path(jj) = ...
                    Single_aa_mut_struct_it.num_w_min_path(jj) + 1;

                Single_aa_mut_struct_it.min_path_nt_mut_contained(kk) = jj;
                %recording with minimal non-synonymous mutation contained
                %in this variant

            end   
            
        end
    end

    if (ii >1)
        Single_aa_mut_report = [Single_aa_mut_report; Single_aa_mut_struct_it];
    else
        Single_aa_mut_report = Single_aa_mut_struct_it;

    end

    if ((Single_aa_mut_struct_it.stop_codon == 0)&&...
            (Single_aa_mut_struct_it.min_path_dist == 1)) %only recording
        %non-stop mutations and mutations needing only 1 non-synonymous
        %mutation

        if(sum(Single_aa_mut_struct_it.num_w_min_path) ~=...
                Single_aa_mut_struct_it.num_nt_mut_var)
            disp('ERROR: not all nt variants contain a minimal NS mutation');
        end

        for jj = 1:Single_aa_mut_struct_it.min_path_count
            analysis_aa_mut_info_it = analysis_aa_mut_info_init;
            %initializing structure array to hold info for all single nt
            %variants

            analysis_aa_mut_info_it.SL8ID = ...
                Single_aa_mut_struct_it.SL8_ID + (jj-1)/10; %if there are
            %more than 1 minimum path variants, they are indexed by #.1 for
            %the second, #.2 for the 3rd, etc.

            analysis_aa_mut_info_it.SL8aaID = ...
                Single_aa_mut_struct_it.aa_mut_numeric;

            analysis_aa_mut_info_it.orig_SL8ntID_set =...
                Single_aa_mut_struct_it.nt_mut_numeric(:, ...
                Single_aa_mut_struct_it.min_path_nt_mut_contained == jj);
            %recording the nt numeric code for the nt variants that contain
            %the current minimal non-synonymous mutation

            if (exist("analysis_aa_mut_info", "var"))
                analysis_aa_mut_info = [analysis_aa_mut_info;...
                    analysis_aa_mut_info_it'];
            else
                analysis_aa_mut_info = analysis_aa_mut_info_it;
            end
        end
    end

end

save("Single_aa_mutations.mat", "Single_aa_mut_report", "analysis_aa_mut_info");