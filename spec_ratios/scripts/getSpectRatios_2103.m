function [rbc2gas,bar2gas,linewidth,gas,disRes,gasRes] = getSpectRatios_2103(twix_obj)

%Aryil Bechtel 2021

%{
takes twix_obj and returns:
        -rbc2gas and bar2gas ratios
        -dedicated gas signal spectral linewidth (hz) and area from spectroscopy
        -resnorm of bi-exp fit to dis and gas FIDs

written for calibration sequences

2103: detrends dissolved and gas FIDs using bi-exp fit to FID peaks
%}

theFID = squeeze(double(twix_obj.image())); %get all FIDs
FIDlength = numel(theFID(:,1)); %get number of pts in one FID

% Define number of calibration, number to skip, number of gas, target flip
nCal = 20; % number of flip angle calibration frames following, normally 20
nDis = numel(theFID(1,:)) - (nCal + 2); %end of dis signal
TR = twix_obj.hdr.Config.TR/1e6; %time between rf pulses (in sec)
nAvg = round(1/TR); % Dissolved fids to average. Set to match Elly's 1-sec static. Set to >600 to use all available
%nAvg = nDis - 1;
nSkip = nDis - nAvg - 1;
%nSkip = 0;
Voigt = 1; % Barrier lorentzian for 0, Voigt if 1. 
           % Requires Elly's NMR_fit_v, NMR_mix_v, and NMR_TimeFit_v for Voigt=1
           
% parse out the various fids - dissolved, gas, and calibration
disFID = theFID(:,2:nDis); %all dissolved data, skip first FID
xdataDis = 1:numel(disFID(1:end)); 
gasFIDall = theFID(:,(end-nCal+2):end); %all gas data, skip first FID
xdataGas = 1:numel(gasFIDall(1:end));

%get values and indices of FID peaks
[disPeaks,disX] = max(abs(disFID),[],[1]);
disPeaksX = FIDlength*(0:(size(disFID,2)-1))+disX;
[gasPeaks,gasX] = max(abs(gasFIDall),[],[1]);
gasPeaksX = FIDlength*(0:(size(gasFIDall,2)-1))+gasX;

%fit bi-exponential to FID peaks
fitfunct = @(coefs,xdata)coefs(1)*exp(-xdata*coefs(2))+coefs(3)*exp(-xdata*coefs(4));
myVal = max(abs(disFID(1:end)))/exp(1); %estimate tau1 for both fits
[minVal,I] = min(abs(abs(disFID(1:end)) - myVal));
disTau1 = 1/I;
myVal = max(abs(gasFIDall(1:end)))/exp(1);
[minVal,I] = min(abs(abs(gasFIDall(1:end)) - myVal));
gasTau1 = 1/I;

[disFitParams, disRes] = biExpFIDfit(disPeaks,disPeaksX,disTau1,-1e-50);
[gasFitParams, gasRes] = biExpFIDfit(gasPeaks,gasPeaksX,gasTau1,-1e-50);


%detrend FIDs
disFitArray = disFitParams(1)*exp(-xdataDis*disFitParams(2)) ... 
                +disFitParams(3)*exp(-xdataDis*disFitParams(4));
%disFIDdetrend = abs(disFID(1:end)) ./ disFitArray;
disFIDdetrend = abs(disFID(1:end)) - disFitArray;

gasFitArray = gasFitParams(1)*exp(-xdataGas*gasFitParams(2)) ... 
                +gasFitParams(3)*exp(-xdataGas*gasFitParams(4));
%gasFIDdetrend = abs(gasFIDall(1:end)) ./ gasFitArray;
gasFIDdetrend = abs(gasFIDall(1:end)) - gasFitArray;

%separate detrended data back into indiv. FIDs
disLength = length(disFIDdetrend)/FIDlength;
gasLength = length(gasFIDdetrend)/FIDlength;
disData = reshape(disFIDdetrend,[FIDlength disLength]);
gasDataAll = reshape(gasFIDdetrend,[FIDlength gasLength]);

% extract favorite variables from the header
dwell_time = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1,1};  % appears most robust and generic
dwell_time=dwell_time*1E-9; % convert to seconds for calculations
 
% parse out the various fids - dissolved, gas, and calibration
%gasData = gasDataAll(:,end-nCal+2); %skip a few detrend gas FID then pick one
gasData = theFID(:,end-nCal+1);
if (nSkip+1+nAvg) > nDis
    fprintf('\n Requested averages exceeds available; using %0.0f\n',nDis);
    disData1 = theFID(:,nSkip+1:nDis); % all available dissolved data with skipping
else
    disData1 = theFID(:,nSkip+1:(nSkip+1+nAvg)); % requested averages of dissolved
end %if
disData1_avg = mean(disData1,2);  % average of dissolved after skipping
t = double((0:(length(disData)-1))*dwell_time');

%fit gas spectrum
gasfitObj = NMR_TimeFit(gasData,t,1e-4,-84,30,0,0,10000);
gasfitObj.fitTimeDomainSignal();
    %figure('Name','Gas Phase Analysis')
    %gasfitObj.plotTimeAndSpectralFit;
    %xlim([0 0.01]);
gasfitObj.describe();

%fit disolved spectrum
if Voigt == 1
    %fprintf('Using Voigt Lineshape to fit Barrier\n');
    disfitObj = NMR_TimeFit_v(disData1_avg,t,[1 1 1],[0 -700  -7400],[250 200 30],[0 200 0],[0 0 0],0,length(t)); % first widths lorenzian, 2nd are gauss
else
    %fprintf('Using Lorentzian Lineshape to fit Barrier\n');
    disfitObj = NMR_TimeFit(disData1_avg,t,[1 1 1],[0 -700  -7400],[240 240 40],[0 0 0],0,length(t));
end % if
disfitObj = disfitObj.fitTimeDomainSignal();
%figure('Name','Dissolved Phase Analysis')
%disfitObj.plotTimeAndSpectralFit;
 

rbc2gas = disfitObj.area(1)/gasfitObj.area;
bar2gas = disfitObj.area(2)/gasfitObj.area;
linewidth = gasfitObj.fwhm;
gas = gasfitObj.area;



function [fitparams, resnorm] = biExpFIDfit(FID,xdata,tau1,tau2)

%input: real-valued, 1D series of FIDs, guess for both decay constants
%(tau1 and tau2)

%written for dissolved and dedicated gas cali sequence

fitfunct = @(coefs,xdata)coefs(1)*exp(-xdata*coefs(2))+coefs(3)*exp(-xdata*coefs(4));
guess(1) = mean(FID);
guess(2) = tau1;     
guess(3) = mean(FID);
guess(4) = tau2;

ydata = FID;
 
fitoptions = optimoptions('lsqcurvefit','Display','off');
[fitparams, resnorm]  = lsqcurvefit(fitfunct,guess,xdata,ydata,[],[],fitoptions);

end

end
