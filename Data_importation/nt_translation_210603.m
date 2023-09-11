%Steffen Docken
%3-6-21
%This code takes a codon sequence (in either numeric or alpha form) and 
%translates to the other form

%%%% numeric code: T = 0, C = 1, A = 2, G = 3

function codon_seq_trans = nt_translation_210603(codon_seq_orig)

nt_alpha_order = {'T', 'C', 'A', 'G'}; %ordering of nucleotide for numeric 
%coding

if ischar(codon_seq_orig)
    
    a_alpha = codon_seq_orig(1);
    b_alpha = codon_seq_orig(2);
    c_alpha = codon_seq_orig(3); %getting first, second, and third nucleotide of
    %codon (a, b, and c) in alphabetic coding (T, C, A, or G)
    
    a_num = find(strcmp(nt_alpha_order, a_alpha)) - 1;
    b_num = find(strcmp(nt_alpha_order, b_alpha)) - 1;
    c_num = find(strcmp(nt_alpha_order, c_alpha)) - 1; %translating to numeric
    %code for nucleotides (listed at top of code)
    
    codon_seq_trans = a_num*16 + b_num*4 + c_num + 1; %numeric code for full
    %codon (+1 because 0 is reserved for WT)
    
elseif isnumeric(codon_seq_orig)
    
    a_num = floor((codon_seq_orig -1)/ 16);
    b_num = floor(mod(codon_seq_orig-1, 16)/4);
    c_num = mod(codon_seq_orig - 1, 4);%getting first, second, and third nucleotide of
    %codon (a, b, and c) in numeric coding (0, 1, 2, or 3)
    
    a_alpha = nt_alpha_order{a_num + 1};
    b_alpha = nt_alpha_order{b_num + 1};
    c_alpha = nt_alpha_order{c_num + 1};%translating to alphabetic
    %code for nucleotides (listed at top of code)
    
    codon_seq_trans = [a_alpha, b_alpha, c_alpha];
    
else
    
    disp('ERROR: data type of codon_seq_orig not accepted');
    
end
