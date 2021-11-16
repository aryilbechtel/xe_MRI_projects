function [refVolt,resnorm] = getRefVolt(twix_obj)

%Aryil Bechtel 2021
%takes twix_obj and returns the true reference voltage for Xe MRI

% Define number of calibration, number to skip, number of gas, target flip
FlipTarget = 20; % target flip angle from calibration
nDis = 200;
nSkip = 100; % number of dissolved fids to skip to ensure steady-state
nAvg = 50; % Dissolved fids to average. Set to 50 to match Elly's 1-sec static. Set to >600 to use all available
nCal = 20; % number of flip angle calibration frames following, normally 20Voigt = 1; % Barrier lorentzian for 0, Voigt if 1. Requires Elly's NMR_fit_v, NMR_mix_v, and NMR_TimeFit_v for Voigt=1
Voigt = 1;

% extract favorite variables from the header
dwell_time = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1,1};  % appears most robust and generic

if isfield(twix_obj.hdr.Phoenix, 'sWiPMemBlock')
    % Duke twix file
    if isfield(twix_obj.hdr.Phoenix.sWiPMemBlock,'adFree')
        VRef = twix_obj.hdr.Phoenix.sWiPMemBlock.adFree{4};
        %rf_amp3=1; %dummy override for when using hard pulse
    elseif isfield(twix_obj.hdr.Phoenix.sWiPMemBlock,'alFree')
        % Duke twix file using UVA sequence
        VRef = twix_obj.hdr.Phoenix.sWiPMemBlock.alFree{1};
    else
        disp('WARNING: twix file type not supported, cannot determine reference voltage')
    end 
elseif isfield(twix_obj.hdr.Phoenix, 'sWipMemBlock')
    % UVA twix file
    VRef = twix_obj.hdr.Phoenix.sWipMemBlock.alFree{1};
else
    disp('WARNING: twix file type not supported, cannot determine reference voltage')
end 

dwell_time=dwell_time*1E-9; % convert to seconds for calculations

% parse out the various fids - dissolved, gas, and calibration
theFID = squeeze(double(twix_obj.image()));
disData = theFID(:,1:nDis); % all dissolved data
gasData = theFID(:,end-nCal+1);
if nAvg > nDis
    fprintf('\n Requested averages exceeds available; using %0.0f\n',nDis);
    disData1 = theFID(:,nSkip+1:nDis); % all available dissolved data with skipping
else
    disData1 = theFID(:,nSkip+1:(nSkip+1+nAvg)); % requested averages of dissolved
end %if
disData1_avg = mean(disData1,2);  % average of dissolved after skipping
calData = theFID(:,end-nCal+1:end); % the data left for flip angle calculations.
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
 
% Find Amplitudes faster by using max amplitudes from each FID in calibration
flipCalAmps = max(abs(calData));

% calculate flip angle
fitfunct = @(coefs,xdata)coefs(1)*cos(coefs(2)).^(xdata-1);%+coefs(3);   % cos theta decay
guess(1)=max(flipCalAmps);
guess(2)=20*pi/180;       % just guess 10 degrees
%guess(3)=0;

xdata=1:length(flipCalAmps);
ydata = flipCalAmps;

fitoptions = optimoptions('lsqcurvefit','Display','off');
[fitparams,resnorm]  = lsqcurvefit(fitfunct,guess,xdata,ydata,[],[],fitoptions);
flip_angle=abs(fitparams(2)*180/pi);

% provide reference amplitude and warnings
VRefScaleFactor=FlipTarget/flip_angle; % How much should Ref Voltage Be Scaled?
refVolt = VRef*VRefScaleFactor;


end %end getRefVolt fn