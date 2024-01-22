function [ind_transits_data_pt,ind_transits_last_det_refl,target_T,target_R,target_B,target_D] = get_qrouter_T(experiment,s1_transit_bin,s2_transit_bin,n_transit_bin,b_transit_bin,d_transit_bin)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
num_det_pulses = sum(ismember(experiment, 'a':'z'));
num_exp_pulses = sum(ismember(experiment, 'A':'Z'));
num_tot_pulses = num_det_pulses+num_exp_pulses;

ind_transits_data_pt = [];
ind_transits_last_det_refl = [];

target_B = 0;
target_D = 0;
target_T = 0;
target_R = 0;


for jj=1:length(n_transit_bin)

rep = (length(n_transit_bin{jj})/num_tot_pulses);

if rep >= 2

ind3 = [(num_det_pulses):(num_tot_pulses)]; 
ind4 = ind3;
if rep >= 3
    for p=1:(rep-2)
        ind4 = [ind4;ind3+p*num_tot_pulses];
    end
end

% find where we got refl on last det pulse
if experiment(num_det_pulses)=='n'
    ff = find((s1_transit_bin{jj}(ind4(:,1)) + s2_transit_bin{jj}(ind4(:,1))) == 1); 
else
    ff = find((n_transit_bin{jj}(ind4(:,1)) + b_transit_bin{jj}(ind4(:,1)) + d_transit_bin{jj}(ind4(:,1)) ) == 1); 
end
   

% find data points
if ~isempty(ff)
    ind_transits_last_det_refl = [ind_transits_last_det_refl,jj];
    for aa = 1:length(ff)
        
        all_det = s1_transit_bin{jj}(ind4(ff(aa),2:end)) + s2_transit_bin{jj}(ind4(ff(aa),2:end)) + n_transit_bin{jj}(ind4(ff(aa),2:end)) + b_transit_bin{jj}(ind4(ff(aa),2:end)) + d_transit_bin{jj}(ind4(ff(aa),2:end));
        S = s1_transit_bin{jj}(ind4(ff(aa),2:end)) + s2_transit_bin{jj}(ind4(ff(aa),2:end));
        N = n_transit_bin{jj}(ind4(ff(aa),2:end)) + b_transit_bin{jj}(ind4(ff(aa),2:end)) + d_transit_bin{jj}(ind4(ff(aa),2:end));

% simple version of experiment - not sending 0c at all because it doesn't
% interact with the atom-cavity anyway, making it a 1-click experiment
% 
%         if (sum(all_det([2,4])) == 1) % one click in slot 2 and 4
%             ind_transits_data_pt = [ind_transits_data_pt,jj];
% 
%                 if sum(N([2,4])) == 1
%                     target_T = target_T + 1;
%                     target_B = target_B + b_transit_bin{jj}(ind4(ff(aa),end));
%                     target_D = target_D + d_transit_bin{jj}(ind4(ff(aa),end));
%                 elseif sum(S([2,4])) == 1
%                     target_R = target_R + 1;  
%                 end         
% 
%             
%         end
%      

% SPRINT version of analysis:
        if (sum(all_det([1])) == 1) % one click in slot 1
            ind_transits_data_pt = [ind_transits_data_pt,jj];

                if sum(N([1])) == 1
                    if experiment(num_det_pulses+1)=='N'
                        target_T = target_T + 1;
                    else
                        target_R = target_R + 1;
                    end
                    target_B = target_B + b_transit_bin{jj}(ind4(ff(aa),2));
                    target_D = target_D + d_transit_bin{jj}(ind4(ff(aa),2));
                elseif sum(S([1])) == 1
                    if experiment(num_det_pulses+1)=='N'
                        target_R = target_R + 1; 
                    else
                        target_T = target_T + 1;
                    end

                end         

            
        end
     

 
    end
end


end

end
end
