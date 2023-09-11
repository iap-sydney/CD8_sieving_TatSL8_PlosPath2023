%Steffen Docken
%12-04-21
%Code to record the count of barcodes and SL8 sequences seen in sequencing

%Code edited on 22-4-21 to remove sequences with "Unique" barcodes or SL8 misreads
%from data

function [B, M, W, S, A] = bcode_SL8ID_total_count_210604(Barcode_mut_names_raw_it,...
    Seq_count_raw_it, bcode_ID_vec, SL8_ID_array, T_it, animal_name)

load('Global_params_seq_coding_SL8WTaa.mat', 'SL8_WTaa_seq_alpha');
%loading the sequence of WT SL8 amino acids using alphabetic IDs of amino
%acids

SL8_WTaa_seq_num = zeros(1, 8);
for ii = 1:8
    SL8_WTaa_seq_num(ii) = aa_translation_210603(SL8_WTaa_seq_alpha{ii});
end %recording the sequence of amino acids for WT SL8 using the code 
%defined in Global_params_seq_coding_SL8WTaa.m


num_rows = length(Barcode_mut_names_raw_it); %number of rows of data
num_bcode_it = length(bcode_ID_vec); %number of barcodes in animal
num_SL8_ID_it = length(SL8_ID_array(1,:)); %number of SL8 IDs in animal

Num_seq_runs_it = length(Seq_count_raw_it(1,:)); %number of sequencing runs

S = zeros(1, Num_seq_runs_it);%S_raw is number of sequences for each data 
%time point with barcode and SL8 variant identified
A = zeros(1, Num_seq_runs_it);%A_raw is amplification factor for each data 
%time point (S/T)

B.count = zeros(num_bcode_it,1, Num_seq_runs_it); 
M.count = zeros(1,num_SL8_ID_it, Num_seq_runs_it); 
W.count = zeros(num_bcode_it ,num_SL8_ID_it, Num_seq_runs_it); %initializing arrays to hold barcode,
%SL8 variant, and barcode-SL8 variant combo counts

%% recording counts
for ii = 1:num_rows
    
    bcode_SL8_info1_it = split(Barcode_mut_names_raw_it{ii}, '.');
    
    if (strcmp(bcode_SL8_info1_it{1}(1:6), 'Unique'))
        
        continue %skipping rows with "Unique" barcode
%         bcode_ID = -1; %registering that barcode was not a barcode
%         %previously identified in stock
%         
%         SL8_info1 = strcat(bcode_SL8_info1_it{2:end}); %for 'Unique' barcodes,
%         %SL8 info comes after first '.'
    else
        
        for jj = 3:length(bcode_SL8_info1_it)-1
            bcode_SL8_info1_it{jj} = strcat(bcode_SL8_info1_it{jj},'.');
        end%replacing the '.' in the SL8 code, which separate the amino 
        %acid index and codon sequence
        
        SL8_info1 = strcat(bcode_SL8_info1_it{3:end}); %for known barcodes,
        %SL8 info comes after first '.'
        
        SL8_info2 = split(SL8_info1, ':'); %breaking up var_ID and mutation IDs
        %if they exist.
        
        if (sum(strcmp(SL8_info2{1}, ...
                {'[shortread]', '[badseq]', '[N]', '[indel]'})) > 0)
            continue %not counting sequences with SL8 misread
            
        elseif (strcmp(SL8_info2{1}, 'WT'))
            var_ID = 0; %registering that SL8 is WT
            
        else
            
            if (length(SL8_info2) <2)
                disp(strcat("ERROR: SL8 not captured and not recorded in processed data correctly for row ",...
                    num2str(ii), " of animal ", animal_name));
                %if there is no aa mutation info after the first ':', then the
                %sequence should be one of the options previously checked for.
                %If not, an error is reported
            end
            
            var_ID = str2double(SL8_info2{1});
        end
        
        aa_mut_ID = zeros(8,1);
        nt_mut_ID = zeros(8,1); %initializing column vectors to 
        %indicate amino acid mutations and codon mutations (left as 0 if 
        %fully WT
        
        num_nt_mut_it = length(SL8_info2) - 1;
        
        for jj = 1:num_nt_mut_it %will skip adding any mutation info if 
            %there wasn't any, because SL8 label was just 'WT'
            
            mut_info_it = split(SL8_info2{jj + 1}, '.');
            aa_ind = str2double(mut_info_it{1}); %identifying the codon of
            %which aa was mutated
            
            codon_num = nt_translation_210603(mut_info_it{2}); %getting the
            %numeric code for the mutated nt sequence 
            
            nt_mut_ID(aa_ind) = codon_num; %recording the numeric code 
            %for the mutated nt sequence in the entry for the correct amino
            %acid
            
            aa_num = nt_aa_translation_210604(mut_info_it{2}); %getting the 
            %numeric code for the amino acid coded by the mutated codon
            
            if (aa_num ~= SL8_WTaa_seq_num(aa_ind))
                aa_mut_ID(aa_ind) = aa_num; %recording the numeric code 
            %for the mutated aa in the entry for the correct amino
            %acid. if the nt mutation is just a synonymous mutation, the
            %corresponding aa_mut_ID_it location is left at 0 to indicate
            %the aa is still wild-type.
            end
            
        end
        
        bcode_ID = str2double(bcode_SL8_info1_it{2}); %registering barcode ID
        
        SL8_ID = [var_ID; aa_mut_ID; nt_mut_ID]; %combining var_ID, 
        %aa_mut_ID, and nt_mut_ID into overall SL8 ID
        
    end
    
    bcode_ind = find(bcode_ID_vec == bcode_ID); 
    SL8_ID_ind = find(~sum((SL8_ID_array - SL8_ID).^2));%finding indices corresponding
    %to current barcode and SL8 sequence. by subtracting SL8_ID from
    %SL8_ID_array, the correct column will be all zeros (by squaring the 
    %entries of the resultant array, that insures all non-zero entries are
    %>0, so summing a non-zero column will give a >0 value). find(~sum())
    %then finds the column that sums to zero.
    
    if ((length(bcode_ind) ~= 1)||(length(SL8_ID_ind) ~= 1))
        disp(strcat("ERROR: non-unique or no matching barcode or SL8 ID found for row ",...
                num2str(ii), " of animal ", animal_name));
    end
    
    for jj = 1:Num_seq_runs_it
        B.count(bcode_ind, 1, jj) = B.count(bcode_ind, 1, jj) + ...
            Seq_count_raw_it(ii, jj); %increasing count of this barcode at
        %this time point
        M.count(1, SL8_ID_ind, jj) = M.count(1, SL8_ID_ind, jj) + ...
            Seq_count_raw_it(ii, jj); %increasing count of this SL8 variant at
        %this time point
        
        W.count(bcode_ind, SL8_ID_ind, jj) = W.count(bcode_ind, SL8_ID_ind, jj) + ...
            Seq_count_raw_it(ii, jj); %increasing count of this barcode-SL8 
        %variant at this time point. Adding because the same variant can
        %arise due to different mutations
    end
    
end

%recording S and A for each sequencing run
for ii = 1:Num_seq_runs_it
    S_B = sum(B.count(:,1,ii));
    S_M = sum(M.count(1,:,ii));
    S_W = sum(sum(W.count(:,:,ii)));
    
    if ((S_B ~= S_M)||(S_B ~= S_W))
        disp(strcat("ERROR: Total number of sequences not equal for B, M, and W in animal ", ...
            animal_name, ", column ", num2str(ii)));
    end
    
    S(ii) = S_B;
    A(ii) = S_B/T_it(ii); %recording S and A
end

%% recording LODs for and if barcodes and SL8 variants are above LOD


B.LOD = zeros(size(B.count));
B.LOD_ind = zeros(size(B.count));

M.LOD = zeros(size(M.count));
M.LOD_ind = zeros(size(M.count));

W.LOD = zeros(size(W.count));
W.LOD_ind = zeros(size(W.count));

for ii = 1:Num_seq_runs_it
    
    [B.LOD(:,1,ii), M.LOD(1,:,ii), W.LOD(:,:,ii)] = ...
        bcode_SL8_LOD_210604(B.count(:,1, ii), M.count(1,:,ii), ...
        S(ii), A(ii));
    
end

%keyboard

B.LOD_ind = B.count > B.LOD; %recording if each barcode is 
%above LOD at each time point
M.LOD_ind = M.count > M.LOD; %recording if each SL8 variant
%is above LOD at each time point
W.LOD_ind = W.count > W.LOD; %recording if each SL8 variant on each barcode
%is above LOD at each time point


