%Steffen Docken
%3-6-21
%codes to translate amino acid IDs between letter name and numeric name

%%%%% numeric code: A = 1, C = 2, D = 3, E = 4, F = 5, G = 6, H = 7, I = 8,
%K = 9, L = 10, M = 11, N = 12, P = 13, Q = 14, R = 15, S = 16, T = 17, 
%V = 18, W = 19, Y = 20, Stop = 21

function aa_trans = aa_translation_210603(aa_orig)

aa_alpha_order = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M',...
    'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'Stop'};

if ischar(aa_orig)
    aa_trans = find(strcmp(aa_alpha_order, aa_orig));

elseif isnumeric(aa_orig)
    aa_trans = aa_alpha_order{aa_orig};

else
    disp('ERROR: data type of aa_orig not accepted');
    
end