%{
Function to get MOXE (Chang) rbc/bar hemoglobin correction scaling factor

Inputs:
    b: signal normalization constant
    delta: barrier thickness in m
    d: total septal wall thickness in meters
    t_x: capillary transit time in seconds
    T_exch: principal exchange time constant (s) per eq 5 of MOXE
    nterms: number of terms to include in series summations
    myTime: time point in rbc2bar CSSR recovery
    hgbRef: reference hemoglobin concentration (g/dL)
    hgbIn: input hgb to be scaled

Outputs:
    rbc2barFactor: hgb correction factor based on rbc/bar value at hgbRef

Aryil Bechtel 2021
%}

function [rbc2barFactor,rbcFactor,barFactor] = correctHgbMOXE(b,delta,d,t_x,T_exch,nterms,myTime,hgbRef,hgbIn)

lambda_RBC = 0.19; % solubility in RBCs -  Ladefoged J, Physiol Med Biol 1967;12:353?358.
lambda_p = 0.09; % solubility in plasma - Ladefoged J, Physiol Med Biol 1967;12:353?358.

%calculate reference eta and input eta from hgbRef and hgbIn
hgbToHct = 2.953/100; %conversion factor of hgb g/dl to hct
hctRef = hgbRef*hgbToHct;
etaRef = lambda_RBC*hctRef/(lambda_RBC*hctRef + lambda_p*(1-hctRef));
hctIn = hgbIn*hgbToHct;
etaIn = lambda_RBC*hctIn/(lambda_RBC*hctIn + lambda_p*(1-hctIn));

syms n t x % define symbolic variables
S_RBC1(t) = (1-2*delta/d)*t/t_x - (8*T_exch/(pi^2*t_x))*...
    symsum(cos((2*n+1)*pi*delta/d).*(1-exp(-t*(2*n+1)^2/T_exch))./(2*n+1)^4,n,[0 nterms]);
S_RBC2(t) = (1-t/t_x)*((1-2*delta/d)-(8/pi^2)*...
    symsum(cos((2*n+1)*pi*delta/d).*exp(-t*(2*n+1)^2/T_exch)./(2*n+1)^2,n,[0 nterms]));

%use etaRef and etaIn to calc two RBC signals
S_RBC_ref = @(t)b*etaRef*(2*S_RBC1(t) + S_RBC2(t)); % return S_RBC as a function handle dependent only on t
S_RBC_in = @(t)b*etaIn*(2*S_RBC1(t) + S_RBC2(t));

S_bar1(t) = 2*delta/d - (8/pi^2)*...
    symsum((1-cos((2*n+1)*pi*delta/d)).*exp(-t*(2*n+1)^2/T_exch)./(2*n+1)^2,n,[0 nterms]);
S_bar2(t) = (1-2*delta/d)*t/t_x - (8*T_exch/(pi^2*t_x))*...
    symsum(cos((2*n+1)*pi*delta/d).*(1-exp(-t*(2*n+1)^2/T_exch))./(2*n+1)^4,n,[0 nterms]);
S_bar3(t) = S_RBC2(t);

%use etaRef and etaIn to calc two bar signals
S_bar_ref = @(t)b*(S_bar1(t) + (1-etaRef)*(2*S_bar2(t) + S_bar3(t))); % return S_bar as a function handle dependent only on t
S_bar_in = @(t)b*(S_bar1(t) + (1-etaIn)*(2*S_bar2(t) + S_bar3(t)));

rbc2barRef = S_RBC_ref(myTime)/S_bar_ref(myTime);
rbc2barIn = S_RBC_in(myTime)/S_bar_in(myTime);

rbc2barFactor = double(rbc2barIn/rbc2barRef);
rbcFactor = double(S_RBC_in(myTime)/S_RBC_ref(myTime));
barFactor = double(S_bar_in(myTime)/S_bar_ref(myTime));

end