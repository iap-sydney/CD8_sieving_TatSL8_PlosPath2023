%Steffen Docken
%13-4-21
%Code to collect indices for best sequencing runs (most templates) of short
%and SL8 sequencing runs separately

function [final_short_Plasma_ind, final_short_LNMC_ind, final_short_PBMC_ind,...
    final_SL8_Plasma_ind, final_SL8_LNMC_ind, final_SL8_PBMC_ind] = Best_seq_run_indices_210604(dpi_raw,...
    sample_raw, seq_raw, T_raw, S_raw, cutoff_vec)

Template_Plasma_cutoff = cutoff_vec(1);
Template_DNA_cutoff = cutoff_vec(2);
Sequence_Plasma_cutoff = cutoff_vec(3);
Sequence_DNA_cutoff = cutoff_vec(4);

unique_dpi = unique(dpi_raw);
    
best_seq_run_short_Plasma_ind = [];
best_seq_run_short_LNMC_ind = [];
best_seq_run_short_PBMC_ind = [];
best_seq_run_SL8_Plasma_ind = [];
best_seq_run_SL8_LNMC_ind = [];
best_seq_run_SL8_PBMC_ind = []; %these arrays will hold the indices of the best
%sequencing runs for each time point and sample type (for when sequence
%runs were redone) for short and SL8 of each sample type, separately

for jj = 1:length(unique_dpi)
    dpi_ind_it = find(dpi_raw == unique_dpi(jj));

    for kk = 1:3
        switch kk
            case 1
                sample_type_it = 'Plasma';
            case 2
                sample_type_it = 'LNMC';
            case 3
                sample_type_it = 'PBMC';
        end

        for ll = 1:2
            switch ll
                case 1
                    seq_type_it = 'short';
                case 2
                    seq_type_it = 'SL8';
            end
            dpi_sample_ind_it = []; %This array will hold the indices for 
            %the current dpi and sample type combination

            for mm = 1:length(dpi_ind_it)
                if (strcmp(sample_raw{dpi_ind_it(mm)}, sample_type_it)&&...
                        (strcmp(seq_raw{dpi_ind_it(mm)}, seq_type_it)))
                    dpi_sample_ind_it = [dpi_sample_ind_it, dpi_ind_it(mm)];
                    %adding the current index to the list of indices for
                    %the current dpi and sample type combination if it is
                    %the right sample type
                end
            end

            if (~isempty(dpi_sample_ind_it))
                [~, most_templates_ind]= max(T_raw(dpi_sample_ind_it));
                
                if (length(most_templates_ind) > 1) %if multiple runs have 
                    %the same number of templates
                    
                    dpi_sample_ind_it = dpi_sample_ind_it(most_templates_ind); %redefining 
                    %dpi_sample_ind_it to only include those with the max
                    %number of Templates
                    
                    [~, most_templates_ind]= max(S_raw(dpi_sample_ind_it));
                    %getting indices for run with the highest number of
                    %sequences
                end

                %% adding the index of the sequencing run with the most 
                %templates to the list of sequencing runs to be used
                switch ll
                    case 1 %short sequencing
                        switch kk
                            case 1 %Plasma sequencing
                                best_seq_run_short_Plasma_ind = [best_seq_run_short_Plasma_ind, ...
                                    dpi_sample_ind_it(most_templates_ind)]; 
                                
                            case 2 %LNMC sequencing
                                best_seq_run_short_LNMC_ind = [best_seq_run_short_LNMC_ind, ...
                                    dpi_sample_ind_it(most_templates_ind)]; 
                                
                            case 3 %PBMC sequencing
                                best_seq_run_short_PBMC_ind = [best_seq_run_short_PBMC_ind, ...
                                    dpi_sample_ind_it(most_templates_ind)]; 
                        end
                    case 2 %SL8 sequencing
                        switch kk
                            case 1 %Plasma sequencing
                                best_seq_run_SL8_Plasma_ind = [best_seq_run_SL8_Plasma_ind, ...
                                    dpi_sample_ind_it(most_templates_ind)];
                                
                            case 2 %LNMC sequencing
                                best_seq_run_SL8_LNMC_ind = [best_seq_run_SL8_LNMC_ind, ...
                                    dpi_sample_ind_it(most_templates_ind)];
                                
                            case 3 %PBMC sequencing
                                best_seq_run_SL8_PBMC_ind = [best_seq_run_SL8_PBMC_ind, ...
                                    dpi_sample_ind_it(most_templates_ind)];
                        end

                end

            end
        end
        
    end
end

suf_T_Plasma_ind = find(T_raw > Template_Plasma_cutoff);
suf_S_Plasma_ind = find(S_raw > Sequence_Plasma_cutoff);
suf_T_S_Plasma_ind = intersect(suf_T_Plasma_ind, suf_S_Plasma_ind);
%indices of all sequencing runs with Template and Sequence numbers above
%the Plasma cut offs (includes any DNA also above the cut off, but those
%will be eliminated when collecting just best indices)

suf_T_DNA_ind = find(T_raw > Template_DNA_cutoff);
suf_S_DNA_ind = find(S_raw > Sequence_DNA_cutoff);
suf_T_S_DNA_ind = intersect(suf_T_DNA_ind, suf_S_DNA_ind);
%indices of all sequencing runs with Template and Sequence numbers above
%the DNA cut offs (includes any Plasma also above the cut off, but those
%will be eliminated when collecting just best indices)

final_short_Plasma_ind = intersect(suf_T_S_Plasma_ind, best_seq_run_short_Plasma_ind);
final_short_LNMC_ind = intersect(suf_T_S_DNA_ind, best_seq_run_short_LNMC_ind);
final_short_PBMC_ind = intersect(suf_T_S_DNA_ind, best_seq_run_short_PBMC_ind);

final_SL8_Plasma_ind = intersect(suf_T_S_Plasma_ind, best_seq_run_SL8_Plasma_ind);
final_SL8_LNMC_ind = intersect(suf_T_S_DNA_ind, best_seq_run_SL8_LNMC_ind);
final_SL8_PBMC_ind = intersect(suf_T_S_DNA_ind, best_seq_run_SL8_PBMC_ind);