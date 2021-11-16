%Aryil Bechtel 2021
%{
open a series of patient folders 
-> Xe_Dynamic_Cali (129Xe MRI) .dat files
-> extract variables:
    -rbc/gas and bar/gas ratios
    -linewidth of dedicated gas spectral peak
    -area under gas spectral fit

if more than one cali sequence for a subject, use last one in directory
%}

clear;close all;

mainDir = 'Duke_18_22/Xenon_lab/project_spec_ratios/data/additional_subjects/'; %directory with patient folders
%mainDir = '~/OneDrive - Duke University/Raw Data/2019 Q4 Raw Files/';
folders = dir(mainDir);
directoryNames = {folders([folders.isdir]).name};
directoryNames = directoryNames(~ismember(directoryNames,{'.','..'}));

%declare output arrays
directoryNamesCropped = directoryNames;
rbc2gas=zeros(1,size(directoryNames,2));
bar2gas=zeros(1,size(directoryNames,2));
linewidth=zeros(1,size(directoryNames,2));
gas=zeros(1,size(directoryNames,2));
disRes=zeros(1,size(directoryNames,2));
gasRes=zeros(1,size(directoryNames,2));


%loop through folders in mainDir
for i=1:size(directoryNames,2)
    
    %check if directory has the files we want
    contents = dir(strcat(mainDir,char(directoryNames(i))));
    checkCali = 0;
    for k = 1:size(contents,1)
        if contains(contents(k).name(),'Xe_Dynamic_Cali')
            checkCali=checkCali+1;
        end
    end
    
    if  checkCali>0
        %get file names
            %input of '1' to files().name selects first occurence of
            %sequence in directory
            
            %inputof 'end' to files().name selects last occurance of
            %sequence in directory
        files = dir([mainDir,char(directoryNames(i)),'/','*','Xe_Dynamic_Cali','*']);
        XefileName=char(strcat(mainDir,directoryNames(i),'/',files(end).name));

        %create twix object
        Xetwix_obj = mapVBVD(XefileName);

        % dealing with fields that may or may not be in header for all sequences
        %****************Xe twix object ************************
        if(~isfield(Xetwix_obj.hdr.Meas, 'flTransRefAmpl'))
              header.imagetwix_obj.hdr.Meas.flTransRefAmpl = 0;
        end

        if(~isfield(Xetwix_obj.hdr.MeasYaps.sTXSPEC.aRFPULSE{1,1}, 'flAmplitude'))
              Xetwix_obj.hdr.MeasYaps.sTXSPEC.aRFPULSE{1,1}.flAmplitude = 0;
        end

        if(~isfield(Xetwix_obj.hdr.Meas, 'GradDelayTime'))
              Xetwix_obj.hdr.Meas.GradDelayTime = 0;
        end

        if(~isfield(Xetwix_obj.hdr.Config, 'NCha'))
              Xetwix_obj.hdr.Config.NCha = 0;
        end

        if(~isfield(Xetwix_obj.hdr.Config, 'nCha'))  % if nCha also isn't a field, then just one channel
              Xetwix_obj.hdr.Config.nCha = 1;
        end


        if(~isfield(Xetwix_obj.hdr.Config, 'NSlc'))
              Xetwix_obj.hdr.Config.NSlc = 0;
        end

        if(~isfield(Xetwix_obj.hdr.Config, 'NoOfFourierPartitions'))  % if reading 3D data
              Xetwix_obj.hdr.Config.NoOfFourierPartitions = 1;
        end

        if(~isfield(Xetwix_obj.hdr.Meas, 'lGain'))
              Xetwix_obj.hdr.Meas.lGain = 0;
        end

        %get variables
        [rbc2gas(i),bar2gas(i),linewidth(i),gas(i),disRes(i),gasRes(i)] = getSpectRatios_2103(Xetwix_obj);
        
        fprintf('rbc2gas: ');
        fprintf('%4.4f',rbc2gas(i));
        fprintf(newline);

        fprintf('bar2gas: ');
        fprintf('%4.4f',bar2gas(i));
        fprintf(newline);

        fprintf('dedicated gas linewidth: ');
        fprintf('%4.4f',linewidth(i));
        fprintf(newline);
        
        fprintf('dissolved resnorm: ');
        fprintf('%4.9f',disRes(i));
        fprintf(newline);
        
        fprintf('gas resnorm: ');
        fprintf('%4.9f',gasRes(i));
        fprintf(newline);
    else
        directoryNamesCropped=setdiff(directoryNamesCropped,directoryNames(i));
    end
    
end

bar2gas(bar2gas==0)=[];
rbc2gas(rbc2gas==0)=[];
linewidth(linewidth==0)=[];
gas(gas==0)=[];



