%Steffen Docken
%9-11-21
%Code to fit Poisson model for if variants detected in reactivation are
%predicted by variants before ART

clear
close all

sig_AIC_cutoff = 10; %minimum difference in AIC for significant evidence of
%a better model

N_hypercube = 100;
ln_r_range = [-40; 0];
ln_m_range = [-60; 5];
ln_b_range = [0; 5];
%parameters for latin hypercube sampling of initial parameter estimates

save_paper_figs = 1;
paper_fig_location = 'Figures/';
output_location = 'Poisson_model_variant_reactivation/';


%% bins and vectors for plot of reactivation model fit.
prop_bin_size = 0.25; %log10
prop_bin_min = -4;
prop_bin_max= 0;
prop_bin_edges = logspace(prop_bin_min, prop_bin_max, ...
    (prop_bin_max - prop_bin_min)/prop_bin_size + 1);   

prop_bin_plot = [1e-5, prop_bin_edges(1:end-1)*10^(prop_bin_size/2)];
%center of bars for plotting fraction of variants detected in rebound


%% load data
load UPenn_FTY_barcode_SL8_Mutation_data_210607.mat

Num_animals = length(Data_SL8_Plasma);

%% check correct variants included
load Poisson_model_variant_reactivation/Single_aa_mutations.mat

addpath('Data_importation/')
addpath('General_Functions/')
addpath('Poisson_model_variant_reactivation/');

Data_struct_init.animal_ID = {};
Data_struct_init.ART2_start_SL8_prop = [];%will contain proportion at 
% initiation of ART-2
Data_struct_init.ATI2_SL8_react = [];%will contain if variant detected in 
% reactivation of ATI-2 (1) or not (0)
Data_struct_init.SL8_var_list = []; %will contain the list of SL8 variants
%corresponding to the data held in the other elements of the structure
Data_struct_init.ART2_seq_dpi = zeros(1); %dpi of sequencing at start of 
% second round of ART
Data_struct_init.ATI2_timing = zeros(1);% days post rebound of sequencing 
% after second ATI

SL8_var_interest_react_data(Num_animals) = Data_struct_init;


Model_fit_init.ATI_2_LogL = nan(4, Num_animals); 
Model_fit_init.ATI_2_r = nan(4, Num_animals); 
Model_fit_init.ATI_2_b = nan(4, Num_animals); 
Model_fit_init.ATI_2_m = nan(4, Num_animals);%Loglikelihoods and 
%parameter best fits values for models r_k,j = r_j; r_j + b_j*q_k,j;
% r_j + m_j*f_k; r_j + b_j*q_k,j + m_jf_k for each animal in ATI-2
%reactivation

%definition of structure to contain unrestricted and restricted model 
%definitions for each test of if parameter is important and p-value for in
%each animal
Model_IDs = [1, 0, 0;... %r_j
    1, 0, 1;... %r_j + b_j*q_k,j
    1, 1, 0;...%r_j + m_j*f_k
    1, 1, 1];%r_j + b_j*q_k,j + m_j*f_k


%array to hold results:

Model_AIC_table_names = {'animal', 'ART2_seq_dpi', 'ATI2_rebound_seq_dpi',...
    'AIC_minimal', 'AIC_react', 'AIC_evol', 'AIC_full', 'Inconclusive'};
Model_AIC_table_types = {'cellstr', 'double', 'double', 'double', 'double',...
    'double', 'double', 'double'};

Model_AICs_init = table('Size', [Num_animals, length(Model_AIC_table_names)],...
        'VariableTypes', Model_AIC_table_types,...
        'VariableNames', Model_AIC_table_names);

Model_Lratio_names = {'animal',...
    'react_simple', 'evolve_simple', 'react_full', 'evolve_full'};

Model_Lratio_p_init = table('Size', [Num_animals, length(Model_Lratio_names)],...
        'VariableTypes', Model_AIC_table_types([1, 4:7]),...
        'VariableNames', Model_Lratio_names); %

% getting list of SL8 variants being considered (from 
% Single_aa_mutations.mat)
num_SL8_var_interest = length(analysis_aa_mut_info);
analysis_aa_mut_IDs = zeros(1, num_SL8_var_interest);
for ii = 1:num_SL8_var_interest
    analysis_aa_mut_IDs(ii) = analysis_aa_mut_info(ii).SL8ID;
end

%% collecting day 14 data
mut_prop_d14 = zeros(Num_animals, num_SL8_var_interest);

for ii = 1:Num_animals
    [Data_SL8_Plasma_animal, ~, ~] = ...
        SL8_single_mut_variant_count_221107(Data_SL8_Plasma(ii),...
        Data_SL8_LNMC(ii), Data_SL8_PBMC(ii)); %Data strucs for current
    %animal. will only contain SL8 variants of interest for this analysis

    preART_ind = find(Data_SL8_Plasma_animal.dpi == 14, 1);
    %find index of last time point at or before ART
    if (isempty(preART_ind))
        disp('ERROR: no day 14 data')
        continue %skip if no timepoint found
    end

    [~, preART_SL8_Plasma_animal_frac] = background_subtraction(Data_SL8_Plasma_animal.V.count(:,:,preART_ind),...
        Data_SL8_Plasma_animal.V.LOD(:,:,preART_ind), Data_SL8_Plasma_animal.S(preART_ind));
    % getting VL proportions (after background subtraction) for SL8 
    % variants of interest on day 14

    preART_SL8_Plasma_animal_frac(preART_SL8_Plasma_animal_frac <0) = 0;
    %setting SL8 proportions to 0 for SL8 variants below detection
    
    for kk = 1:num_SL8_var_interest
        SL8_var_ind_it = find(Data_SL8_Plasma_animal.SL8_ID(1,:)== ...
            analysis_aa_mut_IDs(kk));
        if isempty(SL8_var_ind_it)
            disp(strcat('ERROR: Single aa var not found for animal ', ...
                Data_SL8_Plasma_animal.animal_ID))
        end
        mut_prop_d14(ii,kk) = preART_SL8_Plasma_animal_frac(SL8_var_ind_it); 
    end
    %recording the fraction of 
    %the viral load preART that was each SL8 variant (after
    %background subtraction)
        
end

avg_mut_prop_d14 = mean(mut_prop_d14, 1); %averaging day 14 VL
% proportions across all animals

%% collecting animal data

Model_AICs_it = Model_AICs_init; %initializing table for results
Model_deltaAICs_it = Model_AICs_init; %initializing table for results

Model_fit_it = Model_fit_init; %initializing LogL_struct for this 
%iteration

Model_Lratio_p_it = Model_Lratio_p_init; % initializing likelihood ratio 
% p-test table

plot_animal_counter = 1; %initializing counter for plot

for ii = 1:Num_animals
    [Data_SL8_Plasma_animal, ~, ~] = ...
        SL8_single_mut_variant_count_221107(Data_SL8_Plasma(ii),...
        Data_SL8_LNMC(ii), Data_SL8_PBMC(ii)); %Data strucs for current
    %animal. will only contain SL8 variants detected in Plasma above LOD in
    %at least one time point

    Model_AICs_it.animal(ii) = {Data_SL8_Plasma_animal.animal_ID};
    Model_deltaAICs_it.animal(ii) = {Data_SL8_Plasma_animal.animal_ID};
    Model_Lratio_p_it.animal(ii) = {Data_SL8_Plasma_animal.animal_ID};
    
    VL_data_animal = VL_data(ii);

    SL8_var_interest_react_data(ii).animal_ID = Data_SL8_Plasma_animal.animal_ID;
    SL8_var_interest_react_data(ii).SL8_var_list = analysis_aa_mut_IDs;

    if strcmp(VL_data_animal.Group, 'A')
        Group_ART_days = Group_A_ART_days;
    elseif strcmp(VL_data_animal.Group, 'B')
        Group_ART_days = Group_B_ART_days;
    else
        disp('ERROR: group not identified');
    end

    ART2_seq_ind = find(Data_SL8_Plasma_animal.dpi == ...
        Group_ART_days(3), 1);
    %find index of last time point at or before ART

    ATI2_react_ind = find((strcmp(Data_SL8_Plasma_animal.phase,...
        'ATI-2')&(Data_SL8_Plasma_animal.timing < 11)), 1, 'first');
    %finding first sequencing within 10 days of detected reactivation
    %(timing for first day of detected reactivation is listed as 1)

    if ((isempty(ART2_seq_ind))||(isempty(ATI2_react_ind)))
        Model_Lratio_p_it.react_simple(ii) = nan;
        Model_Lratio_p_it.evolve_simple(ii) = nan;
        Model_Lratio_p_it.react_full(ii) = nan;
        Model_Lratio_p_it.evolve_full(ii) = nan;
        continue %skip this animal
    end

    if ((~isempty(VL_data_animal.CD8depl_dpi))&&...
            (Data_SL8_Plasma_animal.dpi(ATI2_react_ind) > VL_data_animal.CD8depl_dpi))
        Model_Lratio_p_it.react_simple(ii) = nan;
        Model_Lratio_p_it.evolve_simple(ii) = nan;
        Model_Lratio_p_it.react_full(ii) = nan;
        Model_Lratio_p_it.evolve_full(ii) = nan;
        continue %skip animal if reactivation sequencing is after CD8 depletion
    end

    for jj = 1:2
        switch jj
            case 1 %ART-2 sequencing
                ind_it = ART2_seq_ind;
            case 2 %ATI-2 rebound seq
                ind_it = ATI2_react_ind;
        end

        [~, Plasma_animal_frac_it] = background_subtraction(Data_SL8_Plasma_animal.V.count(:,:,ind_it),...
            Data_SL8_Plasma_animal.V.LOD(:,:,ind_it), Data_SL8_Plasma_animal.S(ind_it));
        % getting VL proportions (after background subtraction) for SL8 
            % variants of interest on current day

        Plasma_animal_frac_it(Plasma_animal_frac_it < 0) = 0;
        %setting SL8 proportions to 0 for SL8 variants below detection

        switch jj
            case 1 %ART-2 sequencing
                ART2_SL8_Plasma_animal_frac = Plasma_animal_frac_it;
                
                SL8_var_interest_react_data(ii).ART2_seq_dpi = ...
                    Data_SL8_Plasma_animal.dpi(ind_it);    
                
                Model_AICs_it.ART2_seq_dpi(ii) = Data_SL8_Plasma_animal.dpi(ind_it);
                Model_deltaAICs_it.ART2_seq_dpi(ii) = Data_SL8_Plasma_animal.dpi(ind_it);
            case 2 %ATI-2 rebound seq
                ATI2_react_SL8_Plasma_animal_frac = Plasma_animal_frac_it;

                SL8_var_interest_react_data(ii).ATI2_timing = ...
                    Data_SL8_Plasma_animal.timing(ind_it);   

                Model_AICs_it.ATI2_rebound_seq_dpi(ii) = Data_SL8_Plasma_animal.dpi(ind_it);
                Model_deltaAICs_it.ATI2_rebound_seq_dpi(ii) = Data_SL8_Plasma_animal.dpi(ind_it);
        end
        
    end

    for kk = 1:num_SL8_var_interest
        SL8_var_ind_it = find(Data_SL8_Plasma_animal.SL8_ID(1,:)== ...
            analysis_aa_mut_IDs(kk));
        if isempty(SL8_var_ind_it)
            disp(strcat('ERROR: Single aa var not found for animal ', ...
                Data_SL8_Plasma_animal.animal_ID))
        end
        
        SL8_var_interest_react_data(ii).ART2_start_SL8_prop(kk) = ...
            ART2_SL8_Plasma_animal_frac(SL8_var_ind_it); %recording the fraction of 
        %the viral load preART that was each SL8 variant (after
        %background subtraction)
    
        SL8_var_interest_react_data(ii).ATI2_SL8_react(kk) = ...
            ATI2_react_SL8_Plasma_animal_frac(SL8_var_ind_it) > 0;
        %testing if each variant detected in rebound in this animal

    end

    %% fitting Poisson models and testing significance of parameters being 
    % non-zero

    react_det_it = SL8_var_interest_react_data(ii).ATI2_SL8_react;
    preART2_prop = SL8_var_interest_react_data(ii).ART2_start_SL8_prop;
    %indicator of if SL8 variant was detected in rebound and proportion of
    %the viral load at re-initiation of ART consisting of that SL8 variant.

    for Model_ind = 1:4

        model_def_it = Model_IDs(Model_ind,:);

        %LHS:
        rng default %for reproducibility
   
        IC_vec = lhsdesign(N_hypercube, 3); 
        %hypercube sampling of [0,1]^3.

        IC_vec(:,1) = unifinv(IC_vec(:,1), ln_r_range(1), ln_r_range(2));
        %converting the sampling from [0,1] to the "realistic range" for
        %ln(r_j)

        IC_vec(:,2) = unifinv(IC_vec(:,2), ln_m_range(1), ln_m_range(2));
        %converting the sampling from [0,1] to the "realistic range" for
        %ln(m_j)

        IC_vec(:,3) = unifinv(IC_vec(:,3), ln_b_range(1), ln_b_range(2));
        %converting the sampling from [0,1] to the "realistic range" for
        %ln(b_j)

        IC_vec = IC_vec(:,model_def_it == 1); %only parameters used in current model

        neg_logLike_it = @(ln_params)-loglikelihood_Poisson_restricted(exp(ln_params), ...
            model_def_it, react_det_it, preART2_prop, ...
            avg_mut_prop_d14);
        
        IC_logL_array = zeros(N_hypercube, 1);
        IC_ln_pars_array = zeros(N_hypercube, sum(model_def_it));

        parfor ll = 1:N_hypercube

            IC_it = IC_vec(ll,:); %only params for current model

            for kk = 1:10 %recursively optimizing from the last best fit
                [ln_pars_opt_IC, min_neglogLike_IC] = ...
                    fminsearch(neg_logLike_it, IC_it);

                IC_it = ln_pars_opt_IC;
            end

            IC_logL_array(ll) = -min_neglogLike_IC;
            IC_ln_pars_array(ll, :) = ln_pars_opt_IC;
        end

        [~, max_IC_ind] = max(IC_logL_array);
        max_logLike = IC_logL_array(max_IC_ind);
        ln_pars_opt = IC_ln_pars_array(max_IC_ind, :);
                    
        %%reporting results
        Model_fit_it.ATI_2_LogL(Model_ind, ii) = max_logLike;
        Model_fit_it.ATI_2_r(Model_ind, ii) = exp(ln_pars_opt(1));
        ln_pars_opt = ln_pars_opt(2:end); %changing list of opt. paramters so
            %m_j (if fit) will be first parameter. Otherwise, b_j (if
            %fit) is now first parameter
        
        if (model_def_it(2) == 1) %if m_j != 0; Only updating if parameter
            %included in model. Then changing list of opt. paramters so
            %b_j (if fit) will be first parameter
            Model_fit_it.ATI_2_m(Model_ind, ii) = exp(ln_pars_opt(1));
            ln_pars_opt = ln_pars_opt(2:end);
        end
        
        if (model_def_it(3) == 1) %if b_j != 0; Only updating if parameter
            %included in model.
            Model_fit_it.ATI_2_b(Model_ind, ii) = exp(ln_pars_opt(1));
        end
        
        Model_AICs_it{ii, Model_ind + 3} = ...
            AIC_calc_generic(max_logLike, sum(model_def_it), ...
            num_SL8_var_interest, 1); %calculated corrected AIC based
        %on max loglikelihood, number of parameters, number of data
        %points, and '1' indicates to use corrected AIC
        
    end
    
    [minAIC, minAIC_ind] = min(Model_AICs_it{ii, 4:7}); %finding which
    %model is the minimum AIC
    
    delta_AIC_vec = Model_AICs_it{ii, 4:7} - minAIC;
    Model_deltaAICs_it{ii, 4:7} = delta_AIC_vec;
    

    [~, Model_Lratio_p_it.react_simple(ii)] = lratiotest(Model_fit_it.ATI_2_LogL(2, ii),...
        Model_fit_it.ATI_2_LogL(1, ii), 1);
    [~, Model_Lratio_p_it.evolve_simple(ii)] = lratiotest(Model_fit_it.ATI_2_LogL(3, ii),...
        Model_fit_it.ATI_2_LogL(1, ii), 1);
    [~, Model_Lratio_p_it.react_full(ii)] = lratiotest(Model_fit_it.ATI_2_LogL(4, ii),...
        Model_fit_it.ATI_2_LogL(3, ii), 1);
    [~, Model_Lratio_p_it.evolve_full(ii)] = lratiotest(Model_fit_it.ATI_2_LogL(4, ii),...
        Model_fit_it.ATI_2_LogL(2, ii), 1);

    %% plotting reactivation model fit
    prop_bin_frac_det = zeros(1, length(prop_bin_plot));
    %resetting vector that will contain the fraction of variants
    %detected in the rebound.

    for jj = 1:length(prop_bin_plot)
        if (jj == 1)
            bin_var_ind_it = find(preART2_prop < prop_bin_edges(jj));
            %indices of variants below detection threshold
        else
            bin_var_ind_it = find((preART2_prop >= prop_bin_edges(jj-1))&...
                (preART2_prop < prop_bin_edges(jj)));
        end
        num_bin_var_tot_it = length(bin_var_ind_it); %total number of variants
        %whose pre-ART proportion were in the current range
        if(num_bin_var_tot_it >0)
            num_bin_var_det_it = sum(react_det_it(bin_var_ind_it));
        %number of variants whose pre-ART proportion were in the current
        %range that were detected early ATI-2
            prop_bin_frac_det(jj) = num_bin_var_det_it/num_bin_var_tot_it;
        end
    end

    prob_react_vec = react_Poisson_Model_full_211118(prop_bin_plot,...
        Model_fit_it.ATI_2_r(2,ii), 0, ...
        Model_fit_it.ATI_2_b(2,ii), 0);
    %generating expected probability of reactivation vector based on
    %best fit parameters of reservoir model.

    figure(1)
    subplot(4,2,plot_animal_counter)
    bar(log10(prop_bin_plot), prop_bin_frac_det, 'c');
    %set(gca, 'XScale', 'log')
    hold on
    preART_prop_plot = preART2_prop;
    preART_prop_plot(preART_prop_plot==0) = prop_bin_plot(1); %setting 
    %undetected to what will be used to indicate LOD
    plot(log10(preART_prop_plot), react_det_it, '*k',...
        log10(prop_bin_plot), prob_react_vec, '-r');
    %semilogx(preART_prop, react_det_it, '*k',...
    %    prop_bin_edges, prob_react_vec, '-r');
    hold off
    yticks(0:.25:1)
    yticklabels({'0%', '25%', '50%', '75%', '100%'});
    xticks(-5:0)
    xticklabels({'<LOD', '10^{-4}', '10^{-3}', '10^{-2}',...
        '10^{-1}', '1'})
    title(Data_SL8_Plasma_animal.animal_ID, 'Interpreter','none')

    plot_animal_counter = plot_animal_counter + 1;



end

Model_logL_it = Model_AICs_it;
Model_logL_it{1:end,4:7} = Model_fit_it.ATI_2_LogL';

if (save_paper_figs == 1)
    %getting saved export style
    style_name = 'Documents1';
    s=hgexport('readstyle', style_name);
    %applying export style and saving figure
    fig_name = strcat(paper_fig_location,...
        'Single_aa_mut_Poisson_React_Model_fit.eps');
    s.Format = 'eps';
    hgexport(gcf, fig_name, s);
end


Poisson_data_1aa_mut_SL8 = SL8_var_interest_react_data;
Poisson_fit_1aa_mut_SL8 = Model_fit_it;
Poisson_Model_AICs_1aa_mut_SL8 = Model_AICs_it;
Poisson_Model_deltaAICs_1aa_mut_SL8 = Model_deltaAICs_it;
Poisson_Model_logL_1aa_mut_SL8 = Model_logL_it;
Poisson_Model_lr_pval_1aa_mut_SL8 = Model_Lratio_p_it;



save(strcat(output_location, "Poisson_variant_react_model_output_221109.mat"),...
    "Poisson_data_1aa_mut_SL8", "Poisson_fit_1aa_mut_SL8");


writetable(Poisson_Model_AICs_1aa_mut_SL8, strcat(output_location, ...
    "Poisson_Model_AICs_221109.xlsx"), "Sheet","Single_aa_mut_SL8_AIC");

writetable(Poisson_Model_deltaAICs_1aa_mut_SL8, strcat(output_location, ...
    "Poisson_Model_AICs_221109.xlsx"), "Sheet","Single_aa_mut_SL8_deltaAIC");

writetable(Poisson_Model_logL_1aa_mut_SL8, strcat(output_location, ...
    "Poisson_Model_AICs_221109.xlsx"), "Sheet","Single_aa_mut_SL8_logL");

writetable(Poisson_Model_lr_pval_1aa_mut_SL8, strcat(output_location, ...
    "Poisson_Model_AICs_221109.xlsx"), "Sheet","Single_aa_mut_SL8_lr_pval");