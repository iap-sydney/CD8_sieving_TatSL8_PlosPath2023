%Steffen Docken
%30-4-21
%Background subtraction

function [Back_sub_count, Back_sub_prop] =...
    background_subtraction(seq_count, seq_LOD, total_count)

Back_sub_count = seq_count - seq_LOD;

Back_sub_prop = Back_sub_count/total_count;