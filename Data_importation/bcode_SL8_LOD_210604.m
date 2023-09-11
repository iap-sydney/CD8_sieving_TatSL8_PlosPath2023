%Steffen Docken
%04-06-21
%This code generates the LODs for all barcodes, SL8 sequences, and 
%barcode-SL8 combinations at current time point

function [B_LOD_array, M_LOD_array, W_LOD_array] = bcode_SL8_LOD_210604(bcode_count, ...
    SL8_count, S, A)

in_vitro_SL8_mut_rate = 10^-2.54; %fastest rate at which SL8 mutations 
%arise during in vitro sequencing (from Immonen et al. 2020)

num_bcode_it = length(bcode_count);
num_SL8_it = length(SL8_count);

W_LOD_array = zeros(num_bcode_it, num_SL8_it);


LOD_bcodes = max(2, A); %LOD for barcodes is max(2, Amplification) so 
%that a barcode represents at least 1 template and 2 accounts for if
%amplification was less than 2 (so a single in vitro mutation can't cause
%a new barcode to be registered)

B_LOD_array = LOD_bcodes*ones(size(bcode_count));

LOD_SL8_ID = max(LOD_bcodes, in_vitro_SL8_mut_rate*S);%LOD for SL8 variants
%is max(2, Amplification, 10^-2.54*# sequences) so 
%that a SL8 variant represents at least 1 template, 2 accounts for if
%amplification was less than 2 (so a single in vitro mutation can't cause
%a new SL8 variant to be registered), and from Immonen et al. 2020,
%different SL8 variants arise artificially through mutation during the sequencing
%process at different rates, but the highest rate is such that 10^-2.54 of
%all sequences falsely register with a specific SL8 variant.

M_LOD_array = LOD_SL8_ID*ones(size(SL8_count));

%Documentation is in Supplemental Information S1
recomb_rate_fac = 0.13; %Immonen et al. estimated that if the
%fraction of a barcode that is a specific SL8 variant is less than 13% what
%the overall fraction of the viral load composed of that SL8 variant is,
%then it is merely due to recombination and isn't truely on that barcode

for jj = 1:num_bcode_it
    for kk = 1:num_SL8_it

        W_LOD_array(jj, kk) = max([LOD_bcodes,...
            in_vitro_SL8_mut_rate*bcode_count(jj), ...
            recomb_rate_fac*(SL8_count(kk)/S)*bcode_count(jj)]);
        
    end
end