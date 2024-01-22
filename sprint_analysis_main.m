%% Load data

total_duration = 8e6;
max_tt = 7.5e6;
min_tt = 0.5e6;%0e6;


router_exp = 'T'; % do not change - means 1click experiment

date = '20240118';

% date2 = '20240109';

load([date,'\time_atoms.mat'])
load([date,'\time_noatoms.mat'])
% load([date,'\seq_coh_check.mat'])


folder_atom = [date,'\with_atoms\'];
folder_noatom = [date,'\no_atoms\'];

t_atoms = time_atoms(9,:) % eventually need to write function that adds different runs
t_noatoms = time_noatoms(9,:)

load([folder_atom,'\tt_',t_atoms,'.mat'])
load([folder_noatom,'\tt_',t_noatoms,'.mat'])

load([folder_atom,'\seq_',t_atoms,'.mat'])

% find experiment sequence
[ss,ii] = sort([peaks_N,peaks_S]);
heights = [heights_N,heights_S];
experiment = repmat('x',1,length(ss));
experiment(heights(ii)>0.44 & ii>length(peaks_N)) = repmat('s',1,sum(heights(ii)>0.44 & ii>length(peaks_N)));
experiment(heights(ii)<0.44 & ii>length(peaks_N)) = repmat('S',1,sum(heights(ii)<0.44 & ii>length(peaks_N)));
experiment(heights(ii)>0.44 & ii<=length(peaks_N)) = repmat('n',1,sum(heights(ii)>0.44 & ii<=length(peaks_N)));
experiment(heights(ii)<0.44 & ii<=length(peaks_N)) = repmat('N',1,sum(heights(ii)<0.44 & ii<=length(peaks_N)));

% experiment = 'snsnsnsnSNSN';   
% patch so I can find the right bins according to tt_SS and tt_NN
% if experiment(8) == 'n'
%     experiment = 'snsnsnsnSNSN';
% else
%     experiment = 'nsnsnsnsNSNS';
% end
num_det_pulses = sum(ismember(experiment, 'a':'z'));
num_exp_pulses = sum(ismember(experiment, 'A':'Z'));
num_tot_pulses = num_det_pulses+num_exp_pulses;


%% Summing detectors to identify transits
% SS = S1 + S2 and NN = N + B + D
tt_SS = add_detectors(tt_S1,tt_S2); 
tt_NN = add_detectors(add_detectors(tt_B,tt_D),tt_N);
tt_SS_noatom = add_detectors(tt_S1_noatom,tt_S2_noatom);
tt_NN_noatom = add_detectors(tt_B_noatom,tt_D_noatom);



%% Search for transits?
% prompt1 = 'Do you want to search for transits? Y/N [Y]: ';
% search_for_transits = input(prompt1,'s');
% if isempty(search_for_transits)
%     search_for_transits = 'Y';
% end
search_for_transits = 'Y';

%% Use previous edges?
prompt = 'Do you want to use edges from workspace? Y/N [Y]: ';
txt = input(prompt,'s');
if isempty(txt)
    txt = 'Y';
end

%% Take timetags up between min_tt and max_tt


for i = 1:length(tt_SS)
    fcutSS = find(tt_SS{i}<max_tt & tt_SS{i}>min_tt);
    tt_SS{i} = tt_SS{i}(fcutSS);
    fcutNN = find(tt_NN{i}<max_tt & tt_NN{i}>min_tt);
    tt_NN{i} = tt_NN{i}(fcutNN);

    fcutN = find(tt_N{i}<max_tt & tt_N{i}>min_tt);
    tt_N{i} = tt_N{i}(fcutN);
    fcutB = find(tt_B{i}<max_tt & tt_B{i}>min_tt);
    tt_B{i} = tt_B{i}(fcutB);
    fcutD = find(tt_D{i}<max_tt & tt_D{i}>min_tt);
    tt_D{i} = tt_D{i}(fcutD);
    fcutS1 = find(tt_S1{i}<max_tt & tt_S1{i}>min_tt);
    tt_S1{i} = tt_S1{i}(fcutS1);
    fcutS2 = find(tt_S2{i}<max_tt & tt_S2{i}>min_tt);
    tt_S2{i} = tt_S2{i}(fcutS2);
end

for i = 1:length(tt_SS_noatom)
    fcutSS = find(tt_SS_noatom{i}<max_tt & tt_SS_noatom{i}>min_tt);
    tt_SS_noatom{i} = tt_SS_noatom{i}(fcutSS);
    fcutNN = find(tt_NN_noatom{i}<max_tt & tt_NN_noatom{i}>min_tt);
    tt_NN_noatom{i} = tt_NN_noatom{i}(fcutNN);

    fcutN = find(tt_N_noatom{i}<max_tt & tt_N_noatom{i}>min_tt);
    tt_N_noatom{i} = tt_N_noatom{i}(fcutN);
    fcutB = find(tt_B_noatom{i}<max_tt & tt_B_noatom{i}>min_tt);
    tt_B_noatom{i} = tt_B_noatom{i}(fcutB);
    fcutD = find(tt_D_noatom{i}<max_tt & tt_D_noatom{i}>min_tt);
    tt_D_noatom{i} = tt_D_noatom{i}(fcutD);
    fcutS1 = find(tt_S1_noatom{i}<max_tt & tt_S1_noatom{i}>min_tt);
    tt_S1_noatom{i} = tt_S1_noatom{i}(fcutS1);
    fcutS2 = find(tt_S2_noatom{i}<max_tt & tt_S2_noatom{i}>min_tt);
    tt_S2_noatom{i} = tt_S2_noatom{i}(fcutS2);
end


%% Find reflection windows in both South and North
edges = 0:1:sequence_duration; 
bin_edges = [edges-0.5,edges(end)+0.5];

ind_trans_n = find(experiment=='n');
ind_trans_s = find(experiment=='s');
ind_refl_n = ind_trans_s;
ind_refl_s = ind_trans_n;

hist_S = zeros(length(tt_SS),length(edges));
hist_N = zeros(length(tt_NN),length(edges));
tt_SS_mod = cell(1,length(tt_SS));
tt_NN_mod = cell(1,length(tt_NN));



% eee = 0:10:20000; %Dorko check
% diff_S = zeros(length(tt_SS),length(eee)-1); %Dorko check
% diff_N = zeros(length(tt_NN),length(eee)-1); %Dorko check



for i = 1:length(tt_SS)
    tt_SS_mod{i} = mod(tt_SS{i},sequence_duration);
%     diff_S(i,:) = histcounts(diff(tt_SS{i}),eee);%Dorko check
    hist_S(i,:) = histcounts(mod(tt_SS{i},sequence_duration),bin_edges);
    tt_NN_mod{i} = mod(tt_NN{i},sequence_duration);
%     diff_N(i,:) = histcounts(diff(tt_NN{i}),eee);%Dorko check
    hist_N(i,:) = histcounts(mod(tt_NN{i},sequence_duration),bin_edges);
end

% figure('name','tt_SS diff per cycle')
% bar(eee(2:end),sum(diff_S)/length(tt_SS))
% axis([min(eee) sequence_duration 0 1.1*max(sum(diff_S)/length(tt_SS))])
% figure('name','tt_NN diff per cycle')
% bar(eee(2:end),sum(diff_N)/length(tt_NN))
% axis([min(eee) sequence_duration 0 1.1*max(sum(diff_N)/length(tt_SS))])



% no atoms
hist_S_noatom = zeros(length(tt_SS_noatom),length(edges));
hist_N_noatom = zeros(length(tt_NN_noatom),length(edges));
tt_SS_mod_noatom = cell(1,length(tt_SS_noatom));
tt_NN_mod_noatom = cell(1,length(tt_NN_noatom));
% diff_S_noatom = zeros(length(tt_SS_noatom),length(eee)-1); %Dorko check
% diff_N_noatom = zeros(length(tt_NN_noatom),length(eee)-1); %Dorko check
for i = 1:length(tt_SS_noatom)
    tt_SS_mod_noatom{i} = mod(tt_SS_noatom{i},sequence_duration);
    tt_NN_mod_noatom{i} = mod(tt_NN_noatom{i},sequence_duration);
    hist_S_noatom(i,:) = histcounts(mod(tt_SS_noatom{i},sequence_duration),bin_edges);
    hist_N_noatom(i,:) = histcounts(mod(tt_NN_noatom{i},sequence_duration),bin_edges);
%     diff_S_noatom(i,:) = histcounts(diff(tt_SS_noatom{i}),eee); %Dorko check
%     diff_N_noatom(i,:) = histcounts(diff(tt_NN_noatom{i}),eee); %Dorko check

end


% figure('name','tt_SS_noatom diff per cycle')
% bar(eee(2:end),sum(diff_S_noatom)/length(tt_SS_noatom))
% axis([min(eee) sequence_duration 0 1.1*max(sum(diff_S_noatom)/length(tt_SS_noatom))])
% figure('name','tt_NN_noatom diff per cycle')
% bar(eee(2:end),sum(diff_N_noatom)/length(tt_NN_noatom))
% axis([min(eee) sequence_duration 0 1.1*max(sum(diff_N_noatom)/length(tt_SS_noatom))])

histogram_S = sum(hist_S);
histogram_N = sum(hist_N);
histogram_S_noatom = sum(hist_S_noatom);
histogram_N_noatom = sum(hist_N_noatom);

%% choose edges
%% CHOOSE EDGES OF REFLECTION PULSES - SOUTH
figure('Name','South')
hold on
bar(edges,histogram_S)
axis([min(edges) max(edges) 0 0.1*max(histogram_S)]) 
    
if txt == 'N'
    %choose points before and after reflection pulse (the small bumps)
    [xs,ys] = ginput(2*sum(experiment=='n')); 
    xs = round(xs);
end
s_edges = xs;
ss_edges = sort([xs;xs]);
ll = zeros(1,2*length(s_edges));
ll(2:4:end) = 1;
ll(3:4:end) = 1;
ll = 1.1*ll*max(histogram_S);
plot(ss_edges,ll,'linewidth',3) 
axis([min(edges) s_edges(3) 0 1.1*max(histogram_S)])  
% axis([min(edges) max(edges) 0 1.1*max(histogram_S)])
if txt == 'N'
    [xS2,yS2] = ginput(1);
end
t_peak_s = round(xS2);


figure('Name','North')
hold on
bar(edges,histogram_N)
axis([min(edges) max(edges) 0 0.1*max(histogram_N)])
if txt == 'N'
    [xn,yn] = ginput(2*sum(experiment=='s'));  
    xn = round(xn); 
end
n_edges =  xn;
nn_edges = sort([xn;xn]);
ll = zeros(1,2*length(n_edges));
ll(2:4:end) = 1;
ll(3:4:end) = 1;
ll = 1.1*ll*max(histogram_N);
plot(nn_edges,ll,'linewidth',3) 
axis([min(edges) s_edges(3) 0 1.1*max(histogram_N)]) 
if txt == 'N'
    [xN2,yN2] = ginput(1);
end
t_peak_n = round(xN2);
axis([min(edges) max(edges) 0 1.1*max(histogram_N)])


close all

% find shift according to pulses peaks
t_shift_S = double(t_peak_s - peaks_S(1));
t_shift_N = double(t_peak_n - peaks_N(1));
s_edges = s_edges-t_shift_S;
n_edges = n_edges-t_shift_N;


% Align timetags according to the shift
for jj=1:length(tt_SS)
    tt_SS{jj} = tt_SS{jj} - t_shift_S;   
    tt_NN{jj} = tt_NN{jj} - t_shift_N;

    tt_S1{jj} = tt_S1{jj} - t_shift_S; 
    tt_S2{jj} = tt_S2{jj} - t_shift_S;
    tt_N{jj} = tt_N{jj} - t_shift_N;
    tt_B{jj} = tt_B{jj} - t_shift_N;
    tt_D{jj} = tt_D{jj} - t_shift_N;


    tt_SS_mod{jj} = tt_SS_mod{jj} - t_shift_S;
    tt_NN_mod{jj} = tt_NN_mod{jj} - t_shift_N;
    hist_S(jj,:) = histcounts(mod(tt_SS{jj},sequence_duration),bin_edges);
    hist_N(jj,:) = histcounts(mod(tt_NN{jj},sequence_duration),bin_edges);
end
histogram_S = sum(hist_S);
histogram_N = sum(hist_N);
for jj=1:length(tt_SS_noatom)
    tt_SS_noatom{jj} = tt_SS_noatom{jj} - t_shift_S;
    tt_NN_noatom{jj} = tt_NN_noatom{jj} - t_shift_N;
    tt_SS_mod_noatom{jj} = tt_SS_mod_noatom{jj} - t_shift_S;
    tt_NN_mod_noatom{jj} = tt_NN_mod_noatom{jj} - t_shift_N;

    tt_S1_noatom{jj} = tt_S1_noatom{jj} - t_shift_S;
    tt_S2_noatom{jj} = tt_S2_noatom{jj} - t_shift_S;
    tt_N_noatom{jj} = tt_N_noatom{jj} - t_shift_N;
    tt_B_noatom{jj} = tt_B_noatom{jj} - t_shift_N;
    tt_D_noatom{jj} = tt_D_noatom{jj} - t_shift_N;

end



%% Get reflection data from South

t_start_S = s_edges(1:2:end);
t_end_S = s_edges(2:2:end);

tt_SS_refl = cell(1,length(tt_SS));
tt_SS_mod_refl = cell(1,length(tt_SS));
hist_refl_S_check = zeros(length(tt_SS),length(edges));

for i = 1:length(tt_SS)
    cond = zeros(1,length(tt_SS_mod{i}));
    for ii=1:(length(s_edges)/2)
        cond = (cond|(tt_SS_mod{i}>t_start_S(ii) & tt_SS_mod{i}<t_end_S(ii))); 
    end
    ff = find(cond>0);
    tt_SS_refl{i} = tt_SS{i}(ff);
    tt_SS_mod_refl{i} = tt_SS_mod{i}(ff);
    hist_refl_S_check(i,:) = histcounts(tt_SS_mod_refl{i},bin_edges);
end

fig_S_seq = figure
hold on
bar(edges,histogram_S)
bar(edges,sum(hist_refl_S_check))
axis([min(edges) max(edges) 0 0.01*max(histogram_S)])
if txt == 'N'
    [xS,yS] = ginput(2*sum(experiment=='N')); 
    S_edges = round(xS);
end
axis([min(edges) max(edges) 0 1.1*max(histogram_S)])

saveas(fig_S_seq, [date,'\',t_atoms,'_',t_noatoms,'\fig_S_seq.fig']);

% no atoms
tt_SS_refl_noatom = cell(1,length(tt_SS_noatom));
tt_SS_mod_refl_noatom = cell(1,length(tt_SS_noatom));

for i = 1:length(tt_SS_noatom)
    cond = zeros(1,length(tt_SS_mod_noatom{i}));
    for ii=1:(length(s_edges)/2)
        cond = (cond|(tt_SS_mod_noatom{i}>t_start_S(ii) & tt_SS_mod_noatom{i}<t_end_S(ii))); 
    end
    ff = find(cond>0);
    tt_SS_refl_noatom{i} = tt_SS_noatom{i}(ff);
    tt_SS_mod_refl_noatom{i} = tt_SS_mod_noatom{i}(ff);
end


%% Get reflection data from North

t_start_N = n_edges(1:2:end);
t_end_N = n_edges(2:2:end);

tt_NN_refl = cell(1,length(tt_NN));
tt_NN_mod_refl = cell(1,length(tt_NN));
hist_refl_N_check = zeros(length(tt_NN),length(edges));

for i = 1:length(tt_SS)
    cond = zeros(1,length(tt_NN_mod{i}));
    for ii=1:(length(n_edges)/2)
        cond = (cond|(tt_NN_mod{i}>t_start_N(ii) & tt_NN_mod{i}<t_end_N(ii))); 
    end
    ff = find(cond>0);
    tt_NN_refl{i} = tt_NN{i}(ff);
    tt_NN_mod_refl{i} = tt_NN_mod{i}(ff);
    hist_refl_N_check(i,:) = histcounts(tt_NN_mod_refl{i},bin_edges);
end

fig_N_seq = figure
hold on
bar(edges,histogram_N)
bar(edges,sum(hist_refl_N_check))
axis([min(edges) max(edges) 0 0.01*max(histogram_N)])
if txt == 'N'
    [xN,yN] = ginput(2*sum(experiment=='S'));
    N_edges = round(xN);
end
axis([min(edges) max(edges) 0 1.1*max(histogram_N)])

saveas(fig_N_seq, [date,'\',t_atoms,'_',t_noatoms,'\fig_N_seq.fig']);

% no atoms
tt_NN_refl_noatom = cell(1,length(tt_NN_noatom));
tt_NN_mod_refl_noatom = cell(1,length(tt_NN_noatom));

for i = 1:length(tt_NN_noatom)
    cond = zeros(1,length(tt_NN_mod_noatom{i}));
    for ii=1:(length(n_edges)/2)
        cond = (cond|(tt_NN_mod_noatom{i}>t_start_N(ii) & tt_NN_mod_noatom{i}<t_end_N(ii))); 
    end
    ff = find(cond>0);
    tt_NN_refl_noatom{i} = tt_NN_noatom{i}(ff);
    tt_NN_mod_refl_noatom{i} = tt_NN_mod_noatom{i}(ff);
end


edges_seq = sort([s_edges;S_edges;n_edges;N_edges]);


%% check B and D port per cycle
tt_B_mod = cell(1,length(tt_B));
tt_D_mod = cell(1,length(tt_B));
B_per_cycle = zeros(1,length(tt_B));
D_per_cycle = zeros(1,length(tt_B));
hist_B = zeros(length(tt_SS),length(edges));
hist_D = zeros(length(tt_SS),length(edges));



for kk=1:length(tt_B)
    tt_B_mod{kk} = mod(tt_B{kk},sequence_duration);
    B_per_cycle(kk) = sum( tt_B_mod{kk}>edges_seq(end-1) & tt_B_mod{kk}<edges_seq(end) );
    hist_B(kk,:) = histcounts(mod(tt_B{kk},sequence_duration),bin_edges);
    tt_D_mod{kk} = mod(tt_D{kk},sequence_duration);
    D_per_cycle(kk) = sum( tt_D_mod{kk}>edges_seq(end-1) & tt_D_mod{kk}<edges_seq(end) );
    hist_D(kk,:) = histcounts(mod(tt_D{kk},sequence_duration),bin_edges);
end

    
histogram_B = sum(hist_B);
histogram_D = sum(hist_D);

figure
bar(edges,histogram_B)
hold on
bar(edges,histogram_D)
axis([edges_seq(end-1) edges_seq(end) 0 10200])


infidBD_per_cycle = D_per_cycle./(B_per_cycle + D_per_cycle);
infidBD = sum(D_per_cycle)/sum(B_per_cycle + D_per_cycle);
infidBD_std = std(infidBD_per_cycle(~isnan(infidBD_per_cycle)));
infidBD_median = median(infidBD_per_cycle(~isnan(infidBD_per_cycle)));
infidBD_totclicks = sum(B_per_cycle + D_per_cycle);
figure; histogram(infidBD_per_cycle);

%% Avg reflections per cycle (nice bump?)
bin_size = 100e3; %[ns]
edges2 = 0:bin_size:total_duration;
bin_edges2 = [edges2 - bin_size/2, edges2(end)+bin_size/2];

hist_refl_S = zeros(length(tt_SS_refl),length(edges2));
for j = 1:length(tt_SS_refl)
    hist_refl_S(j,:) = histcounts(tt_SS_refl{j},bin_edges2);
end
hist_refl_S_sum = sum(hist_refl_S)/length(tt_SS_refl);

hist_refl_N = zeros(length(tt_NN_refl),length(edges2));
for j = 1:length(tt_NN_refl)
    hist_refl_N(j,:) = histcounts(tt_NN_refl{j},bin_edges2);
end
hist_refl_N_sum = sum(hist_refl_N)/length(tt_NN_refl);


refl_per_cycle = sum(hist_refl_N,2) + sum(hist_refl_S,2);
% figure
% bar(edges2,hist_refl_S_sum)
% hold on
% bar(edges2,hist_refl_N_sum)

% no atoms
hist_refl_S_noatom = zeros(length(tt_SS_refl_noatom),length(edges2));
for j = 1:length(tt_SS_refl_noatom)
    hist_refl_S_noatom(j,:) = histcounts(tt_SS_refl_noatom{j},bin_edges2);
end
hist_refl_S_noatom_sum = sum(hist_refl_S_noatom)/length(tt_SS_refl_noatom);

hist_refl_N_noatom = zeros(length(tt_NN_refl_noatom),length(edges2));
for j = 1:length(tt_NN_refl_noatom)
    hist_refl_N_noatom(j,:) = histcounts(tt_NN_refl_noatom{j},bin_edges2);
end
hist_refl_N_noatom_sum = sum(hist_refl_N_noatom)/length(tt_NN_refl_noatom);

% figure
% bar(edges2,hist_refl_S_noatom_sum)
% hold on
% bar(edges2,hist_refl_N_noatom_sum)

figure('Name','Nice Bump?')
subplot(2,1,1);
bar(edges2,hist_refl_S_sum+hist_refl_N_sum)
hold on
bar(edges2,hist_refl_S_noatom_sum+hist_refl_N_noatom_sum,'FaceAlpha',0.4)
subplot(2,1,2);
bar(edges2,(hist_refl_S_sum+hist_refl_N_sum) - (hist_refl_S_noatom_sum+hist_refl_N_noatom_sum))
atom_added_refl_per_cycle = sum((hist_refl_S_sum+hist_refl_N_sum) - (hist_refl_S_noatom_sum+hist_refl_N_noatom_sum));

%% Find transits

% prepare data - number of reflections per pulse for each cycle
edges_det_seq = sort([n_edges;s_edges]);
edges3 = repmat(edges_det_seq,ceil(total_duration/sequence_duration),1);
edge_add = sequence_duration*(ceil((1:(ceil(total_duration/sequence_duration)*length(edges_det_seq)))/length(edges_det_seq))-1);
edges3 = edges3' + edge_add;
S_per_pulse_refl = zeros(length(tt_SS_refl),length(edges3)/2);
N_per_pulse_refl = zeros(length(tt_NN_refl),length(edges3)/2);
S_per_pulse_refl_noatom = zeros(length(tt_SS_refl_noatom),length(edges3)/2);
N_per_pulse_refl_noatom = zeros(length(tt_NN_refl_noatom),length(edges3)/2);


for k=1:length(tt_SS_refl)
    hS = histcounts(tt_SS_refl{k},edges3);
    S_per_pulse_refl(k,:) = hS(1:2:end);
    hN = histcounts(tt_NN_refl{k},edges3);
    N_per_pulse_refl(k,:) = hN(1:2:end);
end

%fix following loop - change index k to index m, change names?

for m=1:length(tt_SS_refl_noatom)
    hS = histcounts(tt_SS_refl_noatom{m},edges3);
    S_per_pulse_refl(k,:) = hS(1:2:end);
    hN = histcounts(tt_NN_refl_noatom{m},edges3);
    N_per_pulse_refl(k,:) = hN(1:2:end);
end


% take all pulses with 2 or more reflections and make it 1 by keeping the first timetag and putting a
% minus sign on the rest of the timetags
[ff_1,ff_2] = find(S_per_pulse_refl>1);
if ~isempty(ff_1)
for jj=1:length(ff_1)
    f2 = find((tt_SS_refl{ff_1(jj)})>edges3(2*ff_2(jj)-1) & (tt_SS_refl{ff_1(jj)})<edges3(2*ff_2(jj)+1));
    tt_SS_refl{ff_1(jj)}(f2(2:end))  = - tt_SS_refl{ff_1(jj)}(f2(2:end));
    f3 = find((tt_SS{ff_1(jj)})>edges3(2*ff_2(jj)-1) & (tt_SS{ff_1(jj)})<edges3(2*ff_2(jj)+1));
    tt_SS{ff_1(jj)}(f3(2:end))  = - tt_SS{ff_1(jj)}(f3(2:end));
end
end

[ff_1,ff_2] = find(N_per_pulse_refl>1);
if ~isempty(ff_1)
for jj=1:length(ff_1)
    f2 = find((tt_NN_refl{ff_1(jj)})>edges3(2*ff_2(jj)-1) & (tt_NN_refl{ff_1(jj)})<edges3(2*ff_2(jj)+1));
    tt_NN_refl{ff_1(jj)}(f2(2:end))  = - tt_NN_refl{ff_1(jj)}(f2(2:end));
    f3 = find((tt_NN{ff_1(jj)})>edges3(2*ff_2(jj)-1) & (tt_NN{ff_1(jj)})<edges3(2*ff_2(jj)+1));
    tt_NN{ff_1(jj)}(f3(2:end))  = - tt_NN{ff_1(jj)}(f3(2:end));
end
end

[ff_1,ff_2] = find(S_per_pulse_refl_noatom>1);
if ~isempty(ff_1)
for jj=1:length(ff_1)
    f2 = find((tt_SS_refl_noatom{ff_1(jj)})>edges3(2*ff_2(jj)-1) & (tt_SS_refl_noatom{ff_1(jj)})<edges3(2*ff_2(jj)+1));
    tt_SS_refl_noatom{ff_1(jj)}(f2(2:end))  = - tt_SS_refl_noatom{ff_1(jj)}(f2(2:end));
    f3 = find((tt_SS_noatom{ff_1(jj)})>edges3(2*ff_2(jj)-1) & (tt_SS_noatom{ff_1(jj)})<edges3(2*ff_2(jj)+1));
    tt_SS_noatom{ff_1(jj)}(f3(2:end))  = - tt_SS_noatom{ff_1(jj)}(f3(2:end));
end
end

[ff_1,ff_2] = find(N_per_pulse_refl_noatom>1);
if ~isempty(ff_1)
for jj=1:length(ff_1)
    f2 = find((tt_NN_refl_noatom{ff_1(jj)})>edges3(2*ff_2(jj)-1) & (tt_NN_refl_noatom{ff_1(jj)})<edges3(2*ff_2(jj)+1));
    tt_NN_refl_noatom{ff_1(jj)}(f2(2:end))  = - tt_NN_refl_noatom{ff_1(jj)}(f2(2:end));
    f3 = find((tt_NN_noatom{ff_1(jj)})>edges3(2*ff_2(jj)-1) & (tt_NN_noatom{ff_1(jj)})<edges3(2*ff_2(jj)+1));
    tt_NN_noatom{ff_1(jj)}(f3(2:end))  = - tt_NN_noatom{ff_1(jj)}(f3(2:end));
end
end


if search_for_transits == 'Y'
%% find clusters
tt_refl = cell(1,length(tt_SS));
for w=1:length(tt_SS)
    v = sort([tt_SS_refl{w},tt_NN_refl{w}]);
    ind = find(v>0,1);
    tt_refl{w} = v(ind:end);
end
tt_refl_noatom = cell(1,length(tt_SS_noatom));
for w=1:length(tt_SS_noatom)
    v2 = sort([tt_SS_refl_noatom{w},tt_NN_refl_noatom{w}]);
    ind2 = find(v2>0,1);
    tt_refl_noatom{w} = v2(ind2:end);
end

min_transit_duration = 350; %[ns]
DT = (300:100:800)';
dt = DT*sqrt(2); %sqrt(2) because the distance is defined in 2d between (t1,t1) and (t2,t2).
clicks = 3:5;
% dt = dt(4);
% clicks = clicks(3);


num_transits = zeros(length(dt),length(clicks));
num_transits_false = zeros(length(dt),length(clicks));
SNR = zeros(length(dt),length(clicks));
transits_tt = cell(length(dt),length(clicks));
transit_duration = cell(length(dt),length(clicks));
transits_cyc_ind = cell(length(dt),length(clicks));

switch router_exp
    case 'R'
        control_T = zeros(length(dt),length(clicks));
        control_R = zeros(length(dt),length(clicks));
        target_T_cT = zeros(length(dt),length(clicks));
        target_R_cT = zeros(length(dt),length(clicks));
        target_T_cR = zeros(length(dt),length(clicks));
        target_R_cR = zeros(length(dt),length(clicks));
        B = zeros(length(dt),length(clicks));
        D = zeros(length(dt),length(clicks));
        B_cR = zeros(length(dt),length(clicks));
        D_cR = zeros(length(dt),length(clicks));
    case 'T'
        target_R = zeros(length(dt),length(clicks));
        target_T = zeros(length(dt),length(clicks));
        target_B = zeros(length(dt),length(clicks));
        target_D = zeros(length(dt),length(clicks));
    case 'Bell'
        fprintf('write something')
end

s1_transit_bin = cell(length(dt),length(clicks));
s2_transit_bin = cell(length(dt),length(clicks));
n_transit_bin = cell(length(dt),length(clicks));
b_transit_bin = cell(length(dt),length(clicks));
d_transit_bin = cell(length(dt),length(clicks));

ind_transits_last_det_refl = cell(length(dt),length(clicks));
ind_transits_data_pt = cell(length(dt),length(clicks));
factor = length(tt_refl)/length(tt_refl_noatom);

for n=1:length(dt)
    for m=1:length(clicks)
    [num_transits(n,m),num_transits_false(n,m),SNR(n,m),transits_tt{n,m},transit_duration{n,m},transits_cyc_ind{n,m}] = find_transits(dt(n),clicks(m),min_transit_duration,tt_refl,tt_refl_noatom);
    [s1_transit_bin{n,m},s2_transit_bin{n,m},n_transit_bin{n,m},b_transit_bin{n,m},d_transit_bin{n,m}] = bin_transits(edges_seq,sequence_duration,transits_tt{n,m},transits_cyc_ind{n,m},tt_S1,tt_S2,tt_N,tt_B,tt_D);
    
    switch router_exp
        case 'R'
            [ind_transits_data_pt{n,m},ind_transits_last_det_refl{n,m},control_T(n,m),control_R(n,m),target_T_cT(n,m),target_R_cT(n,m),target_T_cR(n,m),target_R_cR(n,m),B_cR(n,m),D_cR(n,m),B(n,m),D(n,m)] = get_qrouter_R(experiment,s1_transit_bin{n,m},s2_transit_bin{n,m},n_transit_bin{n,m},b_transit_bin{n,m},d_transit_bin{n,m});
        case 'T'
            [ind_transits_data_pt{n,m},ind_transits_last_det_refl{n,m},target_T(n,m),target_R(n,m),target_B(n,m),target_D(n,m)] = get_qrouter_T(experiment,s1_transit_bin{n,m},s2_transit_bin{n,m},n_transit_bin{n,m},b_transit_bin{n,m},d_transit_bin{n,m});
            % [ind_transits_data_pt{n,m},ind_transits_last_det_refl{n,m},target_T(n,m),target_R(n,m),target_B(n,m),target_D(n,m)] = get_qrouter_T(experiment,s1_transit_bin{n,m},s2_transit_bin{n,m},n_transit_bin{n,m},b_transit_bin{n,m},d_transit_bin{n,m});
        case 'Bell'
%             fprintf('write function')
    end

    end
end


%% one choice of transit condition to plot figures
ind_dt = 3; 
ind_clicks = 1;

edgest = 0:100:2500;
figure 
subplot(2,2,1); histogram(transit_duration{ind_dt,ind_clicks},edgest); title('transit duration histogram')
subplot(2,2,2); histogram(cellfun(@length, transits_tt{ind_dt,ind_clicks})); title('clicks per transit')
first_refl_click = zeros(1,length(transits_tt));
for k=1:length(transits_tt{ind_dt,ind_clicks})
    first_refl_click(k) = transits_tt{ind_dt,ind_clicks}{k}(1);
end
hist_transits = histcounts(first_refl_click,bin_edges2);
subplot(2,2,[3,4]); bar(edges2,hist_transits); title('transits bump')


num_cyc = length(tt_refl);
num_cyc_noatom = length(tt_refl_noatom);



%% Save
[status, msg, msgID] = mkdir([date,'\',t_atoms,'_',t_noatoms]);
save([date,'\',t_atoms,'_',t_noatoms,'\workspace.mat'])

switch router_exp
    case 'R'
        save([date,'\',t_atoms,'_',t_noatoms,'\analysis.mat'],'ind_transits_data_pt','ind_transits_last_det_refl','control_T','control_R','target_T_cT','target_R_cT','target_T_cR','target_R_cR','B','D','B_cR','D_cR')
    case 'T'
        save([date,'\',t_atoms,'_',t_noatoms,'\analysis.mat'],'ind_transits_data_pt','ind_transits_last_det_refl','target_T','target_R','target_B','target_D')
    case 'Bell'
        fprintf('write function')
end

% 
% for n=1:length(dt)
%     for m=1:length(clicks)
%         for k=1:length(s1_transit_bin{n,m})
%             S_transit_bin{n,m}{k} = s1_transit_bin{n,m}{k} + s2_transit_bin{n,m}{k};
%             N_transit_bin{n,m}{k} = n_transit_bin{n,m}{k} + b_transit_bin{n,m}{k} + d_transit_bin{n,m}{k};
%         end
%         if ~isempty(S_transit_bin{n,m})
%         S_transit_bin_sum{n,m} = cellfun(@sum,S_transit_bin{n,m});
%         end
%         if ~isempty(S_transit_bin{n,m})
%         N_transit_bin_sum{n,m} = cellfun(@sum,N_transit_bin{n,m});
%         end
%     end
% end
% diff_vec = [];
% for j = 1:length(S_transit_bin{1,1})
% 
% diff_vec = [diff_vec,diff(find(S_transit_bin{1,1}{j}==1))]
% end
% hh = histogram(diff_vec);
% figure
% plot(hh.Values/sum(hh.Values))


% S_transit_bin_sumsum = cellfun(@sum,S_transit_bin_sum);
% N_transit_bin_sumsum = cellfun(@sum,N_transit_bin_sum);
% 
% ratio = S_transit_bin_sumsum./N_transit_bin_sumsum;
% target_T_tot = target_T_1 + target_T_2 + target_T_3 + target_T_4 + target_T_5 + target_T_6;
% target_R_tot = target_R_1 + target_R_2 + target_R_3 + target_R_4 + target_R_5 + target_R_6;
% T_fid = target_T_tot./(target_R_tot + target_T_tot)
% T_totpts = (target_R_tot + target_T_tot)


end