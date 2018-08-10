% 
% 07/18/10 :Load the data from Mathematica generated data files 
% 
% m : pump mode
% My: pump mode along Y
% 
% example : load_higherorder_overlap(4,64)
% 
% This code may be used or distributed under terms of MIT License.
% This file is part of the PSA-2D-Mode-Structure-Solver project.
% 
function overlap = load_higherorder_overlap( m, My )
 load(['HO_Overlap_FixScale_Pump_',num2str(m),'_Modes_',num2str(My),'.dat'])
 eval(['overlap = HO_Overlap_FixScale_Pump_',num2str(m),'_Modes_',num2str(My),';'])
 return
end
