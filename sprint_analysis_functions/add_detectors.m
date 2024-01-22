function [tt_new] = add_detectors(tt_1,tt_2)
%ADD_DETECTORS Summary of this function goes here
%   Detailed explanation goes here
if length(tt_1) ~= length(tt_2)
    fprintf('timetags must be equal size')
else
    tt_new = cell(1,length(tt_1));
    for i=1:length(tt_1)
        tt_new{i} = sort([tt_1{i},tt_2{i}]);
    end
    
end

end

