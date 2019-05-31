function [conn,vbar,veg,lfp,tspE,tspI,Isynbar,inp,seqs]= NetworkRunSeqt(p,input,NE,NI,T,opt)
%NE, NI,NA, Edc, Idc, Adc, JmpE ,JmpI ,JmpA,gEE,gII,gAA,gEI,gEA,gIE,gIA,gAE,gAI,novar,nonoise, storecurrs)

CE = p.CE;
glE = p.glE;
ElE = p.ElE;
aE = p.aE ;
bE = p.bE;
slpE = p.slpE;
twE = p.twE;
VtE = p.VtE;
VrE = p.VrE;

CI = p.CI;
glI = p.glI;
ElI = p.ElI;
aI = p.aI;
bI = p.bI;
slpI = p.slpI;
twI = p.twI;
VtI = p.VtI;
VrI = p.VrI;

VrevE = p.VrevE;
VrevI = p.VrevI;

% NE excitatory neurons NI inhibitory neurons
% T [=] s duration of sim

%% pre select the sequence
NEseq = [];
nL = 10;
nu = 0;
if opt.seqassign
    while nu<nL
        NEseq = randi(80,[1,nL]) ;
        NEu = unique(NEseq);
        nu = length(NEu);
    end
    NEseq = NEseq+(0:80:(NE-80));
end
seqs = sort(NEseq); %I am progressing CA3 input in time with CA1 index
% so cells with increasing indexes will have progressively delayed inputs



% opt options: nonoise novar noiseprc

if opt.nonoise;
    nfrc = 0;
else
    nfrc = opt.noiseprc/100;
end
gnoiseE = p.gnoiseE*nfrc;
gnoiseI = p.gnoiseI*nfrc;

Edc = p.Edc;
Idc = p.Idc;

Istd = p.DCstdI*Idc;
Estd = p.DCstdE*Edc;
w = random('norm',0,1,[1,NE]);
Edc_dist = w*Estd+Edc;

w = random('norm',0,1,[1,NI]);
Idc_dist = w*Istd+Idc;


if ~isempty(NEseq)
    w = random('norm',0,1,[1,length(NEseq)]);
    w = sort(w,'descend');
    mm = max(Edc_dist)+p.dcbias*Estd;
    Edc_dist(NEseq) = mm*(1+w*0.02);
end
inp.Edc = Edc_dist';
inp.Idc = Idc_dist';

%% inputs

MX = zeros(100,NE); %matrix to sum only some incoming inputs of CA3
for idcn = seqs%1:NE
    kkS = ceil((idcn/NE)*100);
    idX = kkS;
    %   ll = max(1,kkS-1);
    %   lr = min(100,kkS+1);
    %   idX = ll:lr;
    MX(idX, idcn) = 1;
end
%% synapses
display('wire ntwk - ')
tic
if opt.novar
    p.gvarEE = 0;
    p.gvarII = 0;
    p.gvarEI = 0;
    p.gvarIE = 0;
end

W = random('norm',0,1,[NE,NE]);   %mean 0 var 1
mn = p.gmaxEE/NE;
vr = p.gvarEE*mn;
GEE = sqrt(vr)*W+mn; %var(aX) = a^2var(X)
GEE(GEE<0) = 0;
conn.EtoE = GEE;

W = random('norm',0,1,[NI,NI]);  %mean 0 var 1
mn = p.gmaxII/NI;
vr = p.gvarII*mn;
GII = sqrt(vr)*W+mn;
GII(GII<0) = 0;
conn.ItoI = GII;

W = random('norm',0,1,[NE,NI]);
mn = p.gmaxEI/NE;
vr = p.gvarEI*mn;
GEI = sqrt(vr)*W+mn;
GEI(GEI<0) = 0;
conn.EtoI = GEI;

W = random('norm',0,1,[NI,NE]);
mn = p.gmaxIE/NI;
vr = p.gvarIE*mn;
GIE = sqrt(vr)*W+mn;
GIE(GIE<0) = 0;
conn.ItoE = GIE;

toc

%% initialize sim

% allocate simulation output(GIE'*sIE)'*(vE-VrevI)
tspE.times = [];
tspE.celln = [];
tspI.times = [];
tspI.celln = [];
vbar.E = zeros(size(0:.001:T));
vbar.I = zeros(size(0:.001:T));
veg.E = zeros(size(0:.001:T));
veg.I = zeros(size(0:.001:T));
lfp = zeros(size(0:.001:T));

% time
dt = 0.001; % [=] ms integration step
t = 0:dt:1000; % one sec time axis

%peak values of biexps signals
pvsE = exp((1/(1/p.tauEd-1/p.tauEr)*log(p.tauEd/p.tauEr))/p.tauEr)...
    - exp((1/(1/p.tauEd-1/p.tauEr)*log(p.tauEd/p.tauEr))/p.tauEd);
fdE = exp(-dt/p.tauEd); %factor of decay
frE = exp(-dt/p.tauEr); %factor of rise

pvsI = exp((1/(1/p.tauId-1/p.tauIr)*log(p.tauId/p.tauIr))/p.tauIr)...
    - exp((1/(1/p.tauId-1/p.tauIr)*log(p.tauId/p.tauIr))/p.tauId);
fdI = exp(-dt/p.tauId); %factor of decay
frI = exp(-dt/p.tauIr); %factor of rise

pvsIE = exp((1/(1/p.tauIEd-1/p.tauIEr)*log(p.tauIEd/p.tauIEr))/p.tauIEr)...
    - exp((1/(1/p.tauIEd-1/p.tauIEr)*log(p.tauIEd/p.tauIEr))/p.tauIEd);
fdIE = exp(-dt/p.tauIEd); %factor of decay
frIE = exp(-dt/p.tauIEr); %factor of rise

pvsEI = exp((1/(1/p.tauEId-1/p.tauEIr)*log(p.tauEId/p.tauEIr))/p.tauEIr)...
    - exp((1/(1/p.tauEId-1/p.tauEIr)*log(p.tauEId/p.tauEIr))/p.tauEId);
fdEI = exp(-dt/p.tauEId); %factor of decay
frEI = exp(-dt/p.tauEIr); %factor of rise


Eeg = randi(NE,1) ;
veg.ne = Eeg;
Ieg = randi(NI,1) ;
veg.ni = Ieg;

if opt.storecurrs
    Isynbar.ItoE = zeros(NE,length(0:.001:T));
    Isynbar.EtoE = zeros(NE,length(0:.001:T));
else
    Isynbar = [];
end

Einptrace = zeros(NE,length(0:0.001:1));
Iinptrace = zeros(NI,length(0:0.001:1));

jmpE = p.jmpE;
jmpI = p.jmpI;


starts = input.on; %ripple start times

for secNum = 1:T
    tic
    display('second number: ');
    display(secNum);
    
    % make background noises
    display('making noises');
    tic
    Enoise = NaN(NE,length(t));
    for idE = 1:NE
        noise = NoiseGenerator(1,dt);
        Enoise(idE,:) = (noise);
    end
    
    Enoise = gnoiseE*Enoise+repmat(Edc_dist',1,length(Enoise));
    
    Inoise = NaN(NI, length(t));
    for idI = 1:NI
        noise = NoiseGenerator(1,dt);
        Inoise(idI,:) = (noise);
    end
    
    Inoise = gnoiseI*Inoise+repmat(Idc_dist',1,length(Inoise));
    
    clear noise
    
    toc
    display('integrating ODE');
    
    % initialize 1st instant of simulation
    vE = rand([NE,1])*(70+VrE)-70;
    wE = aE*(vE-ElE);
    sE = zeros(NE,1); % outgoing synapses
    sEI = zeros(NE,1);
    
    vI = rand([NI,1])*(70+VrI)-70;
    wI = aI*(vI-ElI);
    sI = zeros(NI,1);
    sIE = zeros(NI,1);
    
    
    if secNum==1
        erE = zeros(NE,1);
        edE = zeros(NE,1);
        erEI = zeros(NE,1);
        edEI = zeros(NE,1);
        erI = zeros(NI,1);
        edI = zeros(NI,1);
        erIE = zeros(NI,1);
        edIE = zeros(NI,1);
    else
        vE = vEend;
        wE = wEend;
        sE = sEend;
        sEI = sEIend;
        vI = vIend;
        wI = wIend;
        sI = sIend;
        sIE = sIEend;
        
        % there can be more than one input per second
        tmin = (secNum-1)*1000-100;
        k1 = find(starts>=tmin,1);
        tmax = secNum*1000+20;
        k2 = find(starts>=tmax,1);
        if isempty(k2)
            k2 = length(starts);
        end
        stsec = starts(k1:k2-1);
        
        bmps = zeros(100,length(t));
        bt = zeros(1,length(t));
        ebt = zeros(1,length(t));
        
        if ~isempty(stsec) %if there is no input in that second, skip
            stsec = stsec-(secNum-1)*1000;
            for idr = 1:length(stsec)
                %inside each ripple event
                rplton = stsec(idr);
                rpltoff = rplton+input.length;
                L = input.length-input.slp-2;
                L0 = 0;%input.slp/input.length;
                L1 = 1;%-L0;
                xL = (L0:1/99:L1)*L;
                tbins = xL;
                tons = rplton+input.slp+2+tbins; % start the Ecells bumps after the I cells are inhibiting already
                toffs = tons+input.length/99;
                for idkk = 1:length(tons)
                    bmps(idkk,:) = bmps(idkk,:)+1./(1+exp(-(t-tons(idkk))/1.5))*1./(1+exp((t-toffs(idkk))/1.5));
                end
                bt = bt+1./(1+exp(-(t-rplton)/input.slp))*1./(1+exp((t-rpltoff)/input.slp));
                ebt = ebt+1./(1+exp(-(t-rplton)/(input.slp/2)))*1./(1+exp((t-rpltoff)/(input.slp/2)));
                
                
            end
        end
        AEX = 5*(MX'*bmps)+repmat(ebt,[NE 1]);
        AIX = repmat(bt,[NI 1]);
        Escale = max(ebt);
        if Escale>0
            ff = jmpE/Escale;
        else
            ff = 0;
        end
        Enoise = Enoise+ff*AEX;
        
        Iscale =max(bt);
        if Iscale>0
            gg = jmpI/Iscale;
        else
            gg = 0;
        end
        Inoise = Inoise+gg*AIX;
        
        Einptrace = [Einptrace ff*AEX(:,1000:1000:end)];
        Iinptrace = [Iinptrace gg*AIX(:,1000:1000:end)];
    end
    
    for idt = 1:length(t)
        if and(mod(idt,1000)==1, idt>1)
            ids = (idt-1)/1000+1+1000*(secNum-1);
            vE(vE>0) = 50;
            vI(vI>0) = 50;
            vbar.E(ids) = sum(vE)/NE;
            vbar.I(ids) = sum(vI)/NI;
            veg.E(ids) = vE(Eeg);
            veg.I(ids) = vI(Ieg);
            
            lfp(ids) =(((GEE'*sE)'*(vE-VrevE)+(GIE'*sIE)'*(vE-VrevI)))/NE;
            if opt.storecurrs
                Isynbar.EtoE(:,ids) = (GEE'*sE).*(vE-VrevE);
                Isynbar.ItoE(:,ids) = (GIE'*sIE).*(vE-VrevI);
            end
            
        end
        % E cells
        idX = find(vE>=0);
        wEst = wE(idX);
        
        IsEtoE = (GEE'*sE).*(vE-VrevE);
        IsItoE = (GIE'*sIE).*(vE-VrevI);
        Iapp = Enoise(:,idt);
        Iion = -glE*(vE-ElE)...
            +glE*slpE*exp((vE-VtE)/slpE)-wE;
        Ielec = 0;
        Isyn = IsEtoE+IsItoE;
        dvdt = (Iapp+Iion-Isyn-Ielec)/CE;
        dwdt = (aE*(vE-ElE)-wE)/twE;
        vE = vE+dvdt*dt;
        wE = wE+dwdt*dt;
        %syn gates evolution
        edE = fdE*edE;
        erE = frE*erE;
        sE = erE-edE;
        edEI = fdEI*edEI;
        erEI = frEI*erEI;
        sEI = erEI-edEI;
        if ~isempty(idX)
            %update dynamic vars
            vE(idX) = VrE;
            wE(idX) = wEst+bE;
            % update syn gates
            edE(idX) = edE(idX)+1/pvsE;
            erE(idX) = erE(idX)+1/pvsE;
            edEI(idX) = edEI(idX)+1/pvsEI;
            erEI(idX) = erEI(idX)+1/pvsEI;
            for idn = 1:length(idX)
                n = idX(idn);
                tspE.times(end+1) =  t(idt)/1000+secNum-1;
                tspE.celln(end+1) = n;
            end
        end
        % I cells
        
        idX = find(vI>=0);
        wIst = wI(idX);
        
        IsEtoI = (GEI'*sEI).*(vI-VrevE);
        IsItoI = (GII'*sI).*(vI-VrevI);
        Iapp = Inoise(:,idt);
        Iion = -glI*(vI-ElI)...
            +glI*slpI*exp((vI-VtI)/slpI)-wI;
        Ielec = 0;
        Isyn = IsEtoI+IsItoI;
        dvdt = (Iapp+Iion-Isyn-Ielec)/CI;
        dwdt = (aI*(vI-ElI)-wI)/twI;
        vI = vI+dvdt*dt;
        wI = wI+dwdt*dt;
        %syn gates evolution
        edI = fdI*edI;
        erI = frI*erI;
        sI = erI-edI;
        edIE = fdIE*edIE;
        erIE = frIE*erIE;
        sIE = erIE-edIE;
        if ~isempty(idX)
            %update dynamic vars
            vI(idX) = VrI;
            wI(idX) = wIst+bI;
            % update syn gates
            edI(idX) = edI(idX)+1/pvsI;
            erI(idX) = erI(idX)+1/pvsI;
            edIE(idX) = edIE(idX)+1/pvsIE;
            erIE(idX) = erIE(idX)+1/pvsIE;
            for idn = 1:length(idX)
                n = idX(idn);
                tspI.times(end+1) =  t(idt)/1000+secNum-1;
                tspI.celln(end+1) = n;
            end
        end
    end
    vEend = vE;
    wEend = wE;
    sEend = sE;
    sEIend = sEI;
    vIend = vI;
    wIend = wI;
    sIend = sI;
    sIEend = sIE;
    % store the states at spike times
    toc
end
inp.Etrace = Einptrace;
inp.Itrace = Iinptrace;


return



