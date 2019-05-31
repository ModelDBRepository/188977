function [ripples,spcount,recruit,filtLFP] = CountRipples(T,lfp,tsp_E,tsp_I,NE,NI,inpseq)
tx = 0:0.001:T; %[=]s
fs = 1000; %Hz sampling of LFP

order = 2;
lowFreq = 100;
hiFreq = 300;
[b,a] = butter(order, [lowFreq hiFreq]/(fs/2), 'bandpass');
f3lfp = filtfilt(b,a,lfp); % for duration
mm = mean(f3lfp(tx<=1));
smm = 5*std(f3lfp(tx<=1)); % define duration thresholds
disp(smm)

lowFreq = 100; %60;
hiFreq = 300; %350;
[b2,a2] = butter(order, [lowFreq hiFreq]/(fs/2), 'bandpass');
flfp = filtfilt(b2,a2,lfp); % for frequency


% find ripples from inputs
L = length(inpseq.on);
ripples.n = L;
for idr = 1:L
    xon = inpseq.on(idr)/1000-.02; %[=]s
    xoff =xon+1.5*(inpseq.length)/1000; %[=]s
    k1 = find(tx>=xon,1);   
    k2 = find(tx>=xoff,1);
    if isempty(k2)
        k2 = length(tx);
    end
    t = tx(k1:k2);
    lfpr = f3lfp(k1:k2);
    
    k1 = find(lfpr>(smm+mm),1);
    k2 = find(lfpr>(smm+mm),1,'last');    
    
    if k2>k1+3
        lfpr = lfpr(k1:k2);
        t = t(k1:k2);
    end
    if isempty(lfpr)
        disp('error');
    else
        ripples.time(idr) = t(1);
        ripples.length(idr) = length(t);
    end
        
        % find ripple frequency and duration
     kon = find(tx>= ripples.time(idr),1);
     if kon+150>length(flfp)     
       [pxx,fx] = pwelch(flfp(kon:end),20,2,400,fs);    
        lfpr = flfp((length(flfp)-150):end);
     else
        [pxx,fx] = pwelch(flfp(kon:kon+150),20,2,400,fs);
        lfpr = flfp(kon:kon+150);
    end
    figure(98)
    plot(fx,pxx);
    hold on;
    
    [~,loc] = max(pxx);
    if ~isempty(loc)
        rf = fx(loc);
    else
        rf = NaN;
    end
     ripples.freq(idr) = rf;
     ripples.amps(idr) = max(lfpr);

end
% spikes per ripple


tEall = tsp_E.times; %[=]s
tIall = tsp_I.times;
for idr = 1:L
    xon = ripples.time(idr);%[=]s
    xoff = xon+ripples.length(idr)/1000;
    
    k1 = find(tEall>=xon,1);
    k2 = find(tEall>=xoff,1);
    nEr = tsp_E.celln(k1:k2);
    b = hist(nEr,1:NE);
    spcount.E(idr,:) = b;
    recruit.E(idr) = length(b(b>0))/NE;
    
    k1 = find(tIall>=xon,1);
    k2 = find(tIall>=xoff,1);
    nIr = tsp_I.celln(k1:k2);
    b = hist(nIr,1:NI);
    spcount.I(idr,:) = b;
    recruit.I(idr) = length(b(b>0))/NI;
end
filtLFP = f3lfp;
return
    



