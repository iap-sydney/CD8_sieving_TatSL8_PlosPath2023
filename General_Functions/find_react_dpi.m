%Steffen Docken
%20-1-23
%Code to identify dpi of reactivation (react_dpi), index of the 
% corresponding time point in the VL data array (react_ind), and indicator 
% of if reactivation was detected or if the animal was censored at the 
% given time point (0 = censored, 1 = reactivation detected; react_det). 
% Definition for reactivation is the 
% earliest of 1) the first of at least two viral load quantifications above 
% limit of detection or 2) viral load quantification of above 1000 
% copies/ml on a single day. In the case that ART was resumed, CD8s were 
% depleted, or animal was necropsied before one reactivation was detected, 
% animals were right censored at the day of said event.

%input: 
% VL_data: contains the VL data of the animal
% ATI_start: first day without ART
% ATI_end: end of ATI for censoring (start of ART, necropsy, or CD8
% depletion

function [react_dpi, react_ind, react_det] = find_react_dpi(VL_data,...
    ATI_start, ATI_end)

LOD = 60; %LOD = 60 copies/ml

VL_det_2days_ind = ((VL_data.VL(1:end-1) > LOD)&(VL_data.VL(2:end)>LOD));
%indicator of if VL is above detection on each day and the subsequent
%measurement

two_det_ind = find((VL_det_2days_ind)&...
    (VL_data.t(1:end-1) >= ATI_start)&...
    (VL_data.t(1:end-1) <= ATI_end), 1, 'first'); %find the first of at 
% least two viral load quantifications above limit of detection 

VL_g1000_ind = find((VL_data.VL > 1000)&...
    (VL_data.t >= ATI_start)&...
    (VL_data.t <= ATI_end), 1, 'first');
%find first viral load quantification of above 1000 
% copies/ml on a single day.

react_ind = min([two_det_ind, VL_g1000_ind]);

if (isempty(react_ind)) %if reactivation not detected
    react_ind = find(VL_data.t <= ATI_end, 1, 'last');
    react_dpi = VL_data.t(react_ind);
    react_det = 0;
else %if reactivation detected
    react_dpi = VL_data.t(react_ind);
    react_det = 1;
end