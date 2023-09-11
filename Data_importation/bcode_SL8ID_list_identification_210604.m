%Steffen Docken
%09-04-21
%Code to identify unique barcodes and SL8 variants seen in sequencing

%Code edited on 22-4-21 to remove sequences with "Unique" barcodes or SL8 misreads
%from data

function [bcode_ID_it, SL8_ID_it] = bcode_SL8ID_list_identification_210604(Barcode_mut_names_raw_it,...
    animal_name)

load('Global_params_seq_coding_SL8WTaa.mat', 'SL8_WTaa_seq_alpha');
%loading the sequence of WT SL8 amino acids using alphabetic IDs of amino
%acids

SL8_WTaa_seq_num = zeros(1, 8);
for ii = 1:8
    SL8_WTaa_seq_num(ii) = aa_translation_210603(SL8_WTaa_seq_alpha{ii});
end %recording the sequence of amino acids for WT SL8 using the code 
%defined in Global_params_seq_coding_SL8WTaa.m

num_rows = length(Barcode_mut_names_raw_it); %number of rows

bcode_ID_it = [];
var_ID_it = []; 
aa_mut_ID_it = [];
nt_mut_ID_it = [];%initializing arrays to hold barcode ID,
%SL8 variant ID, aa mutation ID, and nt mutation ID (will combine variant,
%aa_mut, and nt_mut all into one at end for "full SL8 ID")


for ii = 1:num_rows
    
    bcode_SL8_info1_it = split(Barcode_mut_names_raw_it{ii}, '.');
    
    if (~strcmp(bcode_SL8_info1_it{1}(1:6), 'Unique')) %only including
        %sequences with non-unique barcodes
    
        for jj = 3:length(bcode_SL8_info1_it)-1
            bcode_SL8_info1_it{jj} = strcat(bcode_SL8_info1_it{jj},'.');
        end%replacing the '.' in the SL8 code, which separate the amino 
        %acid index and codon sequence
        
        SL8_info1 = strcat(bcode_SL8_info1_it{3:end}); %for known barcodes,
        %SL8 info comes after second '.'
        
        SL8_info2 = split(SL8_info1, ':'); %breaking up var_ID and mutation IDs
    %if they exist.
    
        if (sum(strcmp(SL8_info2{1}, ...
                {'[shortread]', '[badseq]', '[N]', '[indel]'})) > 0)
            continue %not counting sequences with SL8 misread

        elseif (strcmp(SL8_info2{1}, 'WT'))
            var_ID_it = [var_ID_it, 0]; %recording that SL8 is WT

        else

            if (length(SL8_info2) <2)
                disp(strcat("ERROR: SL8 not captured and not recorded in processed data correctly for row ",...
                    num2str(ii), " of animal ", animal_name));
                %if there is no aa mutation info after the first ':', then the
                %sequence should be one of the options previously checked for.
                %If not, an error is reported
            end

            var_ID_it = [var_ID_it, str2double(SL8_info2{1})];
        end
        
        aa_mut_ID_temp = zeros(8,1);
        nt_mut_ID_temp = zeros(8,1); %initializing column vectors to 
        %indicate amino acid mutations and codon mutations for current row
        %(left as 0 if fully WT). will add to aa_mut_ID_it and nt_mut_ID_it
        
        num_nt_mut_it = length(SL8_info2) - 1;
        
        for jj = 1:num_nt_mut_it %will skip adding any mutation info if 
            %there wasn't any, because SL8 label was just 'WT'
            
            mut_info_it = split(SL8_info2{jj + 1}, '.');
            aa_ind = str2double(mut_info_it{1}); %identifying the codon of
            %which aa was mutated
            
            codon_num = nt_translation_210603(mut_info_it{2}); %getting the
            %numeric code for the mutated nt sequence 
            
            nt_mut_ID_temp(aa_ind) = codon_num; %recording the numeric code 
            %for the mutated nt sequence in the entry for the correct amino
            %acid
            
            aa_num = nt_aa_translation_210604(mut_info_it{2}); %getting the 
            %numeric code for the amino acid coded by the mutated codon
            
            if (aa_num ~= SL8_WTaa_seq_num(aa_ind))
                aa_mut_ID_temp(aa_ind) = aa_num; %recording the numeric code 
            %for the mutated aa in the entry for the correct amino
            %acid. if the nt mutation is just a synonymous mutation, the
            %corresponding aa_mut_ID_it location is left at 0 to indicate
            %the aa is still wild-type.
            end
            
        end
        
        aa_mut_ID_it = [aa_mut_ID_it, aa_mut_ID_temp];
        nt_mut_ID_it = [nt_mut_ID_it, nt_mut_ID_temp]; %adding column to 
        %indicate amino acid mutations and codon mutations (left as 0 if 
        %fully WT
        
        
        bcode_ID_it = [bcode_ID_it; str2double(bcode_SL8_info1_it{2})];
        %recording barcode ID only if both barcode and SL8 detected
    end
    
    
end

bcode_ID_it = unique(bcode_ID_it);%removing duplicates of bcode_ID

SL8_ID_it = [var_ID_it; ...
    aa_mut_ID_it;...
    nt_mut_ID_it]; %combining the SL8 variant code, the amino acid mutation
%code, and the nt mutation code into one

SL8_ID_transpose = unique(SL8_ID_it', 'rows'); %eliminating repeated SL8_IDs
%need to transpose SL8_ID to use the unique feature 'rows', which
%eliminates duplicated rows instead of duplicated values


SL8_ID_it = SL8_ID_transpose';
