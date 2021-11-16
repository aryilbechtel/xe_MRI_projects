%Aryil Bechtel 2021
%{
open a series of patient folders 
-> read BHUTE (1H MRI) and Xe_Dynamic_Cali (129Xe MRI) .dat files
-> extract header variables and 1H ref voltage from BHUTE
-> extract 129Xe ref voltage from Xe_Dynamic

if more than one of each file type for a subject: use first one in
directory
%}

clear;close all;

mainDir = '~/OneDrive - Duke University/Raw Data/2020 Q3-4 Raw Files/'; %directory with patient folders
folders = dir(mainDir);
directoryNames = {folders([folders.isdir]).name};
directoryNames = directoryNames(~ismember(directoryNames,{'.','..'}));

%declare output arrays
directoryNamesCropped = directoryNames;
weight=zeros(1,size(directoryNames,2));
Href_amp=zeros(1,size(directoryNames,2));
true_refVolt=zeros(1,size(directoryNames,2));
resnorm=zeros(1,size(directoryNames,2));

%loop through folders in mainDir
for i=1:size(directoryNames,2)
    
    %check if directory has the files we want
    contents = dir(strcat(mainDir,char(directoryNames(i))));
    checkBHUTE=0;
    checkCali = 0;
    for k = 1:size(contents,1)
        if contains(contents(k).name(),'BHUTE')
            checkBHUTE=checkBHUTE+1;
        end
    end
    for k = 1:size(contents,1)
        if contains(contents(k).name(),'Xe_Dynamic_Cali')
            checkCali=checkCali+1;
        end
    end
    
    if and(checkBHUTE>0,checkCali>0)
        %get file names
            %input of '1' to files().name selects first occurence of
            %sequence in directory
            
            %inputof 'end' to files().name selects last occurance of
            %sequence in directory
        files = dir([mainDir,char(directoryNames(i)),'/','*','BHUTE','*']);
        HfileName = char(strcat(mainDir,directoryNames(i),'/',files(1).name));
        files = dir([mainDir,char(directoryNames(i)),'/','*','Xe_Dynamic_Cali','*']);
        XefileName=char(strcat(mainDir,directoryNames(i),'/',files(1).name));

        %create twix objects
        Htwix_obj = mapVBVD(HfileName);
        Xetwix_obj = mapVBVD(XefileName);


        % dealing with fields that may or may not be in header for all sequences

        %****************H twix object ************************
        if(~isfield(Htwix_obj.hdr.Meas, 'flTransRefAmpl'))
          header.imagetwix_obj.hdr.Meas.flTransRefAmpl = 0;
        end

        if(~isfield(Htwix_obj.hdr.MeasYaps.sTXSPEC.aRFPULSE{1,1}, 'flAmplitude'))
          Htwix_obj.hdr.MeasYaps.sTXSPEC.aRFPULSE{1,1}.flAmplitude = 0;
        end

        if(~isfield(Htwix_obj.hdr.Meas, 'GradDelayTime'))
          Htwix_obj.hdr.Meas.GradDelayTime = 0;
        end

        if(~isfield(Htwix_obj.hdr.Config, 'NCha'))
          Htwix_obj.hdr.Config.NCha = 0;
        end

        if(~isfield(Htwix_obj.hdr.Config, 'nCha'))  % if nCha also isn't a field, then just one channel
          Htwix_obj.hdr.Config.nCha = 1;
        end


        if(~isfield(Htwix_obj.hdr.Config, 'NSlc'))
          Htwix_obj.hdr.Config.NSlc = 0;
        end

        if(~isfield(Htwix_obj.hdr.Config, 'NoOfFourierPartitions'))  % if reading 3D data
          Htwix_obj.hdr.Config.NoOfFourierPartitions = 1;
        end

        if(~isfield(Htwix_obj.hdr.Meas, 'lGain'))
          Htwix_obj.hdr.Meas.lGain = 0;
        end

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

        %extract header variables
        weight(i) = Htwix_obj.hdr.Dicom.flUsedPatientWeight;
        Href_amp(i) = Htwix_obj.hdr.Spice.TransmitterReferenceAmplitude;
        [true_refVolt(i),resnorm(i)]=getRefVolt(Xetwix_obj);

        fprintf('weight: ');
        fprintf('%4.4f',weight(i));
        fprintf(newline);

        fprintf('1H ref amp: ');
        fprintf('%4.4f',Href_amp(i));
        fprintf(newline);

        fprintf('Xe ref folvt: ');
        fprintf('%4.4f',true_refVolt(i));
        fprintf(newline);
    else
        directoryNamesCropped=setdiff(directoryNamesCropped,directoryNames(i));
    end
    
end

weight(weight==0)=[];
Href_amp(Href_amp==0)=[];
true_refVolt(true_refVolt==0)=[];
resnorm(resnorm==0)=[];



