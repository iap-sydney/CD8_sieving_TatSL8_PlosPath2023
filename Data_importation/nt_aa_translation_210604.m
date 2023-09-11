%Steffen Docken
%4-6-21
%This code translates the codon sequence in T, C, A, G notation into
%the amino acid the codon codes for in numeric notation

%%%%%% numeric code: A = 1, C = 2, D = 3, E = 4, F = 5, G = 6, H = 7, I = 8,
%K = 9, L = 10, M = 11, N = 12, P = 13, Q = 14, R = 15, S = 16, T = 17, 
%V = 18, W = 19, Y = 20, Stop = 21

function aa_num = nt_aa_translation_210604(nt_alpha)

if ~ischar(nt_alpha)
    disp('ERROR: data type of aa_orig not accepted');
end

load('Global_params_seq_coding_SL8WTaa.mat', 'Codon_lookup_table');

[aa_num, ~] = find(strcmp(Codon_lookup_table, nt_alpha));
%get the row of the Codon table that nt sequence is in, which is also the
%numeric code for the amino acid
