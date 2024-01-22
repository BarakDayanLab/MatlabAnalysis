function [tot_transits,tot_transits_false,SNR,transits_tt,transit_duration,transits_cyc_ind] = find_transits(dt,clicks,min_transit_duration,tt_refl,tt_refl_noatom)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%% clustering from refl data with atoms
TransitsPerCycle = cell(1,length(tt_refl));
for k = 1:length(tt_refl)
    XY = zeros([length(tt_refl{k}),2]); % create an empty vector to input data to clustering function
    XY(:,1) = double(tt_refl{k});
    XY(:,2) = double(tt_refl{k});
    TransitsPerCycle{k} = clusterXYpoints(XY, dt, clicks, 'point', 'merge'); % database of transits per cycle for given condition. 'Merge' is must.
    if ~isempty(TransitsPerCycle{k})
        dur = zeros(1,length(TransitsPerCycle{k}));
        for w=1:length(TransitsPerCycle{k})
            dur(w) = max(TransitsPerCycle{k}{w}(:,1)) - min(TransitsPerCycle{k}{w}(:,1));
        end
        indt= find(dur>min_transit_duration);
        TransitsPerCycle{k} = TransitsPerCycle{k}(indt);
    end
end
countPerCycle = cellfun(@length, TransitsPerCycle); % transit count per cycle
tot_transits = sum(countPerCycle);

%% clustering from refl data without atoms
TransitsPerCycle_false = cell(1,length(tt_refl_noatom));
for k = 1:length(tt_refl_noatom)
    XY2 = zeros([length(tt_refl_noatom{k}),2]); % create an empty vector to input data to clustering function
    XY2(:,1) = tt_refl_noatom{k};
    XY2(:,2) = tt_refl_noatom{k};
    TransitsPerCycle_false{k} = clusterXYpoints(XY2, dt, clicks, 'point', 'merge'); % database of transits per cycle for given condition. 'Merge' is must.
    if ~isempty(TransitsPerCycle_false{k})
        dur = zeros(1,length(TransitsPerCycle_false{k}));
        for w=1:length(TransitsPerCycle_false{k})
            dur(w) = max(TransitsPerCycle_false{k}{w}(:,1)) - min(TransitsPerCycle_false{k}{w}(:,1));
        end
        indt= find(dur>min_transit_duration);
        TransitsPerCycle_false{k} = TransitsPerCycle_false{k}(indt);
    end
end
countPerCycle_false = cellfun(@length, TransitsPerCycle_false); % transit count per cycle
tot_transits_false = sum(countPerCycle_false);


factor = length(tt_refl)/length(tt_refl_noatom);
SNR = (tot_transits - tot_transits_false*factor)./(tot_transits_false*factor);

transits_cycle_ind = find(countPerCycle>0);
    transits_tt = cell(1,sum(countPerCycle));
    transit_duration = zeros(1,sum(countPerCycle));
    x=1;
    transits_cyc_ind = [];
    for q=1:length(transits_cycle_ind)
          ind_cycle = transits_cycle_ind(q);
          for p=1:countPerCycle(ind_cycle)
          transits_tt{x} = sort(TransitsPerCycle{ind_cycle}{p}(:,1));
          transit_duration(x) = max(transits_tt{x}) - min(transits_tt{x});
          x = x+1;
          transits_cyc_ind = [transits_cyc_ind,transits_cycle_ind(q)];
          end
    end

    
end