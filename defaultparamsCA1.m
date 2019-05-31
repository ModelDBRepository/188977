function [p] = defaultparamsCA1()

%params
p.CE = 200;
p.glE = 10;
p.ElE = -58;
p.aE = 2;
p.bE = 100;
p.slpE = 2;
p.twE = 120;
p.VtE = -50;
p.VrE = -46;

p.CI = 200;
p.glI = 10;
p.ElI = -70;
p.aI = 2;
p.bI = 10;
p.slpI = 2;
p.twI = 30;
p.VtI = -50;
p.VrI = -58;


p.gnoiseE = 80;
p.gnoiseI = 90;

% EtoE 
p.tauEr = 0.5;
p.tauEd = 3.5; 
% EtoB
p.tauEIr = 0.9;
p.tauEId = 3; %EPSP on inrn faster than on principal cells
% BtoB 
p.tauIr = 0.3;
p.tauId = 2;  % Bartos et al 2002
% BtoP is slower than BtoB
p.tauIEr = 0.3;
p.tauIEd = 3.5;

p.gmaxII = 3/(80/100);
p.gmaxEI = 2/(30/100);
p.gmaxIE = 5/(60/100);
p.gmaxEE = p.gmaxEI/10;%1on1 syn strength /percent synch efficacy

p.VrevE = 0; % mV reversal of excitatory syn
p.VrevI = -80;% mV rev of inhibitory syn

p.gvarEE = 0.01;
p.gvarII = 0.01;

p.gvarEI = 0.01;
p.gvarIE = 0.01;

p.DCstdI = 0.1;
p.DCstdE = 0.1;

p.Edc = 40;
p.jmpE = 210;
p.Idc = 180;
p.jmpI = 700;

p.seqsize = 10;
p.dcbias = 2;
return