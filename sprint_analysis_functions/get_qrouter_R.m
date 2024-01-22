function [ind_transits_data_pt,ind_transits_last_det_refl,control_T,control_R,target_T_cT,target_R_cT,target_T_cR,target_R_cR,B_cR,D_cR,B,D] = get_qrouter_R(experiment,s1_transit_bin,s2_transit_bin,n_transit_bin,b_transit_bin,d_transit_bin)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
num_det_pulses = sum(ismember(experiment, 'a':'z'));
num_exp_pulses = sum(ismember(experiment, 'A':'Z'));
num_tot_pulses = num_det_pulses+num_exp_pulses;

ind_transits_data_pt = [];
ind_transits_last_det_refl = [];
B_cR = 0;
D_cR = 0;
B = 0;
D = 0;
control_T = 0;
control_R = 0;
target_T_cT = 0;
target_R_cT = 0;
target_T_cR = 0;
target_R_cR = 0;

for jj=1:length(n_transit_bin)

rep = (length(n_transit_bin{jj})/num_tot_pulses);

if rep >= 2

ind3 = [(num_det_pulses):(num_tot_pulses)]; 
ind4 = ind3;
if rep >= 3
    for p=1:(rep-1)
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
        
        B = B + b_transit_bin{jj}(ind4(ff(aa),end));
        D = D + d_transit_bin{jj}(ind4(ff(aa),end));
         
        if (sum(all_det([1,3])) == 1) & (sum(all_det([2,4])) == 1)
            ind_transits_data_pt = [ind_transits_data_pt,jj];

            if sum(S([1,3])) == 1 % control came out in the right place
                control_R = control_R + 1;
                if sum(N([2,4])) == 1
                    target_R_cR = target_R_cR + 1;
                    B_cR = B_cR + b_transit_bin{jj}(ind4(ff(aa),end));
                    D_cR = D_cR + d_transit_bin{jj}(ind4(ff(aa),end));
                else
                    target_T_cR = target_T_cR + 1;  
                end
            else
                control_T = control_T + 1;
                if sum(N([2,4])) == 1
                    target_R_cT = target_R_cT + 1;
                else
                    target_T_cT = target_T_cT + 1;  
                end
            end
                    

            


        end
        
      
        
 
    end
end
end

end
end
