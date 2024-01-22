function [s1_transit_bin,s2_transit_bin,n_transit_bin,b_transit_bin,d_transit_bin] = bin_transits(edges_seq,sequence_duration,transits_tt,transits_cyc_ind,tt_S1,tt_S2,tt_N,tt_B,tt_D)

%BIN_TRANSITS Summary of this function goes here
%   Detailed explanation goes here
s1_transit_bin = cell(1,length(transits_tt));
s2_transit_bin = cell(1,length(transits_tt));
n_transit_bin = cell(1,length(transits_tt));
b_transit_bin = cell(1,length(transits_tt));
d_transit_bin = cell(1,length(transits_tt));

for jj=1:length(transits_tt)

    fS1 = find(tt_S1{transits_cyc_ind(jj)}>=transits_tt{jj}(1) & tt_S1{transits_cyc_ind(jj)}<=transits_tt{jj}(end) );
    sss1  = tt_S1{transits_cyc_ind(jj)}(fS1);

    fS2 = find(tt_S2{transits_cyc_ind(jj)}>=transits_tt{jj}(1) & tt_S2{transits_cyc_ind(jj)}<=transits_tt{jj}(end) );
    sss2  = tt_S2{transits_cyc_ind(jj)}(fS2);

    fN = find(tt_N{transits_cyc_ind(jj)}>=transits_tt{jj}(1) & tt_N{transits_cyc_ind(jj)}<=transits_tt{jj}(end) );
    nnn = tt_N{transits_cyc_ind(jj)}(fN);

    fB = find(tt_B{transits_cyc_ind(jj)}>=transits_tt{jj}(1) & tt_B{transits_cyc_ind(jj)}<=transits_tt{jj}(end) );
    bbb = tt_B{transits_cyc_ind(jj)}(fB);

    fD = find(tt_D{transits_cyc_ind(jj)}>=transits_tt{jj}(1) & tt_D{transits_cyc_ind(jj)}<=transits_tt{jj}(end) );
    ddd = tt_D{transits_cyc_ind(jj)}(fD);


    start_t_transit = min([floor(double(sss1)/sequence_duration),floor(double(sss2)/sequence_duration),floor(double(nnn)/sequence_duration),floor(double(bbb)/sequence_duration),floor(double(ddd)/sequence_duration)])*sequence_duration;
%     start_t_transit = min([ceil(sss1/sequence_duration),ceil(sss2/sequence_duration),ceil(nnn/sequence_duration),ceil(bbb/sequence_duration),ceil(ddd/sequence_duration)])*sequence_duration;
    sss1 = sss1 - start_t_transit;
    sss2 = sss2 - start_t_transit;
    nnn = nnn - start_t_transit;
    bbb = bbb - start_t_transit;
    ddd = ddd - start_t_transit;
    

    rep = ceil(double(max([sss1,sss2,nnn,bbb,ddd]))/sequence_duration);
    edges_transit = repmat(edges_seq',1,rep);
    edge_add = sequence_duration*(ceil((1:(rep*length(edges_seq)))/length(edges_seq))-1);
    edges_transit = edges_transit + edge_add;


    hs1 = histcounts(sss1,edges_transit);
    s1_transit_bin{jj} = hs1(1:2:end);
    
    hs2 = histcounts(sss2,edges_transit);
    s2_transit_bin{jj} = hs2(1:2:end);
    
    hn = histcounts(nnn,edges_transit);
    n_transit_bin{jj} = hn(1:2:end);
    
    hb = histcounts(bbb,edges_transit);
    b_transit_bin{jj} = hb(1:2:end);
    
    hd = histcounts(ddd,edges_transit);
    d_transit_bin{jj} = hd(1:2:end);


end

end

