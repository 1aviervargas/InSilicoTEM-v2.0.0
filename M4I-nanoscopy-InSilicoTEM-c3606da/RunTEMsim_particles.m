%%  RunTEMsim % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%                                                                         %
%     Script to simulate micrographs under different conditions:          %
%       - with/without phase plate                                        %
%       - defocus range                                                   %
%       - motion blur (radiation damage)                                  %
%       - correction factor of the microscope                             %
%       - electron dose                                                   %
%       - size of the micrograph (pixel number)                           %
%       - pixel size of the micrograph                                    %
%       - distance between particles                                      %  
%                                                                         %
%     To generate the micrographs the script calls the function TEMsim,   %
%       which is an adaptation to the code found in Vulovic, 2013         %
%                                                                         %
%     The main changes to the original code include:                      %
%       -Addition of PHASE PLATE shift                                    %
%       -Optimization of memory usage for parallel processing of          %
%         micrographs                                                     %
%                                                                         %
%    IMPORTANT                                                            %
%       The simulator is based on the DIPimage Matlab Toolbox.            %
%       A working version of DIPimage must be installed on this computer  %
%       in order for the simulation to function properly.                 %
%       Download links can be found at the following address:             %
%           http://www.diplib.org/                                        %
%                                                                         %
%                                                                         %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  

clc;
close all;
clear;
addpath('./extra');

time = tic;

% Parameters
mg = 1;               % Number of micrographs to generate
pp = 1;                 % Phase plate (0 = no; 1 = yes)

%phase_shifts = linspace(0,pi,10);%pi/2*rand(1,10);
%phase_shifts = rand(1,50)*(pi);
%phase_shifts = [pi/2, pi/2, pi/2];
phase_shifts = pi/2;

df_range= [0 2500];    % Defocus range [nm]
mb_series = [0.0];      % Motion blur series (If multiple MB, enter them as a vector)
cf_series = [14];       % Correction factor  (If multiple CF, enter them as a vector)

dose      = 10000000*ones(1,2);
%dose      = [100];       % Electron dose to the specimen [e-/A2]
%pix = 1000;             % Number of pixels 
pix = 800 %500

pixsize = 0.5;          % Pixel size [A]
mindist = 120/pixsize;  %80/PIXSIZE  % Minimum distance between particles divided by pixsize, 
                            % Depends on type of protein (apo ferrtin ~150)
partNum = 1;
cropFraction = 0.5;     % Fraction of the simulated image side to keep (0 < value <= 1)
                            
%pixsize = 1.1;          % Pixel size [A]
%mindist = 150/pixsize;   % Minimum distance between particles divided by pixsize, 
%dir1 =  './Micrographs'; % Select folder where to save micrographs
%dir1 = '/media/jvargas/DATA/jvargas/micSimul/WPP';
dir1 = '/media/jvargas/DATA/jvargas/micSimul/WOPP';
stackName = 'kk.mrcs';
starName = 'kk.star';
stackPath = [dir1 filesep stackName];
starPath = [dir1 filesep starName];
relionMagnification = 100000; % Nominal magnification for metadata export


tot = mg*length(mb_series)*length(cf_series)*length(dose);
count  = 0;
rng('shuffle');
if cropFraction <= 0 || cropFraction > 1
    error('cropFraction must be within (0,1].');
end
finalPix = max(2,floor(pix*cropFraction));
if finalPix > pix
    finalPix = pix;
end
cropOffset = floor((pix - finalPix)/2);
cropStart = cropOffset + 1;
cropEnd = cropStart + finalPix - 1;


disp(...
   [char("###################### Starting new Simulation ######################") newline...
    char("       Micrographs: " + tot) newline...
    char("       Size:        " + pix + "x" + pix) newline...
    char("       Particles:   " + "Randomized") newline...
    char("       Defocus:     [" + df_range(1) + " - " + df_range(2)+"]nm") newline...
    char("       Motion blur: " + mb_series) newline...
    char("       Corr factor: " + cf_series) newline...
    char("       Dose:        " + dose(1)) 'e/A^2' newline...
    '                                 ' char(datetime(now,'ConvertFrom','datenum'))])


% Read PDB File
tic
 disp([newline 'Loading PDB file...'])
 %rpdb = pdbread('./PDBs/1dpx.pdb'); % Select PDB file to read
 rpdb = pdbread('./PDBs/2GTL.pdb'); % Select PDB file to read
 toc
disp(' ')               
                

imageNames = cell(tot,1);
microNames = cell(tot,1);
coordX = zeros(tot,1);
coordY = zeros(tot,1);
originX = zeros(tot,1);
originY = zeros(tot,1);
angleRot = zeros(tot,1);
angleTilt = zeros(tot,1);
anglePsi = zeros(tot,1);
defocusU = zeros(tot,1);
defocusV = zeros(tot,1);
defocusAngle = zeros(tot,1);
voltageKV = zeros(tot,1);
csMM = zeros(tot,1);
phaseShiftDeg = zeros(tot,1);
accDose = zeros(tot,1);
pixelSizeA = zeros(tot,1);
detectorPixUM = zeros(tot,1);
magnifications = zeros(tot,1);
stack = zeros(finalPix,finalPix,tot);

for micro = 1:mg
    for cf = 1 : length(cf_series)
        for d = 1:length(dose)
            for mb = 1 : length (mb_series)    
                try
                    count = count + 1; 
                  
                    ou = tic;
                    
                    % Randomize particles positions for each micrograph
                    disp('Generating random particle positions...')
                    tic
                    %JV
                    %circles = Randomposition(pix, mindist);
                    %particles = length(circles);
                    circles = [0 0 mindist/2]; 
                    particles = 1;
                    %JV
                    toc
                    
                    % Randomize defocus for each micrograph
                    defocus = (df_range(2)-df_range(1))*rand(1,1)+df_range(1);

                    disp([newline newline newline])                  
                    disp( "-----------------------   Progress:   "+count+" / "+tot+" Micrographs   -----------------------")
                    disp(' ')
                    disp(char( "Particles:    "+particles) )
                    disp([char( "Defocus:      "+defocus) 'nm' newline])
                    disp(' ')

                    % Run simulation to generate micrograph
                    [out, simMeta] = TEMsim(defocus,mb_series(mb),cf_series(cf),1,rpdb,particles,dose(d),circles,pix,pixsize,pp,phase_shifts,mindist);
                    %ImageOut = out.series;                    
                    %ImageOut = out.series;
                    delete('./Raw/Particles/*.raw')
                    
                    
                    % Write micrograph to .MRC file
                    %JV
                    %s1 = "MicrographNr"+ (micro) +"_size"+pix+"x"+pix+"_pixsize"+pixsize*100+"_partnr"+particles+"_dose"+dose(d)+"_cf"+cf_series(cf)+"_mb"+mb_series(cf)+"_df"+round(defocus)+"_PP"+pp+".mrcs";
                    %WriteMRC(double(out),pixsize,[dir1 filesep char(s1)]);
                    %JV
                    
                    croppedOut = out(cropStart:cropEnd,cropStart:cropEnd);
                    stack(:,:,count) = croppedOut;
                    imageNames{count} = sprintf('%d@%s',count,stackName);
                    microNames{count} = sprintf('SimMicrograph_%05d',count);
                    accDose(count) = simMeta.dose_e_per_A2;
                    defocusAng = simMeta.defocus_nm*10;
                    defocusU(count) = defocusAng;
                    defocusV(count) = defocusAng;
                    defocusAngle(count) = simMeta.astigmatismAngle_deg;
                    voltageKV(count) = simMeta.voltage_kV;
                    csMM(count) = simMeta.cs_mm;
                    pixelSizeA(count) = simMeta.pixelSize_A;
                    phaseShiftDeg(count) = simMeta.phaseShift_rad(1)*180/pi;
                    magnifications(count) = relionMagnification;
                    detectorPixUM(count) = simMeta.pixelSize_A*relionMagnification/10000;

                    pos = simMeta.positions;
                    if isnumeric(pos) && numel(pos) >= 2
                        if size(pos,1) > 1
                            pos = pos(1,:);
                        end
                        shiftX = pos(1);
                        shiftY = pos(2);
                        centerCoord = (finalPix+1)/2;
                        coordX(count) = shiftX + centerCoord;
                        coordY(count) = shiftY + centerCoord;
                        originX(count) = -shiftX;
                        originY(count) = -shiftY;
                        angleRot(count) = pos(4);
                        angleTilt(count) = pos(5);
                        anglePsi(count) = pos(6);
                    else
                        centerCoord = (finalPix+1)/2;
                        coordX(count) = centerCoord;
                        coordY(count) = centerCoord;
                        originX(count) = 0;
                        originY(count) = 0;
                        angleRot(count) = 0;
                        angleTilt(count) = 0;
                        anglePsi(count) = 0;
                    end
                    disp(' ')
                    disp('Successful particle')
                    

                    
                catch ex
                      
                    % In case the simulator fails to produce a micrograph,
                    % it will skip it and display a warning
                      warning(['Failed Micrograph' newline ex.identifier  ...
                                    newline ex.message      ...
                                    newline ex.stack.name])
                      disp(' ')
                      
                end
                
            end
        end
    end
end

%JV
WriteMRC(double(stack),pixsize,stackPath);
starStruct = struct();
starStruct.rlnImageName = imageNames;
starStruct.rlnMicrographName = microNames;
starStruct.rlnCoordinateX = coordX;
starStruct.rlnCoordinateY = coordY;
starStruct.rlnOriginX = originX;
starStruct.rlnOriginY = originY;
starStruct.rlnAngleRot = angleRot;
starStruct.rlnAngleTilt = angleTilt;
starStruct.rlnAnglePsi = anglePsi;
starStruct.rlnDefocusU = defocusU;
starStruct.rlnDefocusV = defocusV;
starStruct.rlnDefocusAngle = defocusAngle;
starStruct.rlnVoltage = voltageKV;
starStruct.rlnSphericalAberration = csMM;
starStruct.rlnPhaseShift = phaseShiftDeg;
starStruct.rlnDetectorPixelSize = detectorPixUM;
starStruct.rlnMagnification = magnifications;
starStruct.rlnCtfImagePixelSize = pixelSizeA;
starStruct.rlnAccumulatedDose = accDose;
WriteStarFileStruct(starStruct,'particles',starPath);
%JV
disp(' ')
disp(...
   [char("###################### End of Simulation ######################") newline newline...
    '                     ' char(datetime(now,'ConvertFrom','datenum')) newline])
toc(time)
