%Steffen Docken
%07-04-21
%This code generates the look up array for determining which amino acid a
%codon corresponds to.

clear

%will use ismember to check what amino acid codon codes for

A_codons = {'GCT', 'GCC', 'GCA', 'GCG', '', ''};
C_codons = {'TGT', 'TGC', '', '', '', ''};
D_codons = {'GAT', 'GAC', '', '', '', ''};
E_codons = {'GAA', 'GAG', '', '', '', ''};
F_codons = {'TTT', 'TTC', '', '', '', ''};
G_codons = {'GGT', 'GGC', 'GGA', 'GGG', '', ''};
H_codons = {'CAT', 'CAC', '', '', '', ''};
I_codons = {'ATT', 'ATC', 'ATA', '', '', ''};
K_codons = {'AAA', 'AAG', '', '', '', ''};
L_codons = {'TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'};
M_codons = {'ATG', '', '', '', '', ''};
N_codons = {'AAT', 'AAC', '', '', '', ''};
P_codons = {'CCT', 'CCC', 'CCA', 'CCG', '', ''};
Q_codons = {'CAA', 'CAG', '', '', '', ''};
R_codons = {'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'};
S_codons = {'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'};
T_codons = {'ACT', 'ACC', 'ACA', 'ACG', '', ''};
V_codons = {'GTT', 'GTC', 'GTA', 'GTG', '', ''};
W_codons = {'TGG', '', '', '', '', ''};
Y_codons = {'TAT', 'TAC', '', '', '', ''};
Stop_codons = {'TAA', 'TAG', 'TGA', '', '', ''};

Codon_lookup_table = cell(21,6);
for ii = 1:6
    Codon_lookup_table{1,ii} = A_codons{ii};
    Codon_lookup_table{2,ii} = C_codons{ii};
    Codon_lookup_table{3,ii} = D_codons{ii};
    Codon_lookup_table{4,ii} = E_codons{ii};
    Codon_lookup_table{5,ii} = F_codons{ii};
    Codon_lookup_table{6,ii} = G_codons{ii};
    Codon_lookup_table{7,ii} = H_codons{ii};
    Codon_lookup_table{8,ii} = I_codons{ii};
    Codon_lookup_table{9,ii} = K_codons{ii};
    Codon_lookup_table{10,ii} = L_codons{ii};
    Codon_lookup_table{11,ii} = M_codons{ii};
    Codon_lookup_table{12,ii} = N_codons{ii};
    Codon_lookup_table{13,ii} = P_codons{ii};
    Codon_lookup_table{14,ii} = Q_codons{ii};
    Codon_lookup_table{15,ii} = R_codons{ii};
    Codon_lookup_table{16,ii} = S_codons{ii};
    Codon_lookup_table{17,ii} = T_codons{ii};
    Codon_lookup_table{18,ii} = V_codons{ii};
    Codon_lookup_table{19,ii} = W_codons{ii};
    Codon_lookup_table{20,ii} = Y_codons{ii};
    Codon_lookup_table{21,ii} = Stop_codons{ii};
end%will use row index of codon to determine which
%aa it codes for. 
    
SL8_WTaa_seq_alpha = {'S','T','P','E','S','A','N','L'};

save('Global_params_seq_coding_SL8WTaa.mat', 'Codon_lookup_table', 'SL8_WTaa_seq_alpha');