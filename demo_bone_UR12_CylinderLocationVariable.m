%clc
clear
close all

%% File requirement
%addpath('C:\Program Files\MATLAB\R2023b\toolbox\mcx\mcxlab');
addpath('C:\Users\Mohammad\OneDrive - University of Rochester\PhD_Lab\Raman-Monte-Carlo-Extreme-main');
addpath('C:\Program Files\MATLAB\R2023b\toolbox\mcx\mcxtools');
addpath('C:\Program Files\MATLAB\R2023b\toolbox\mcx\iso2mesh');

%% Geometry

%%create cfg.shapes

% tissue labels:0-ambient air,1-skin,2-Diaphysis, 3-Epiphysis, 4-Cortical

unitinmm = 0.22;
cfg.unitinmm = unitinmm;
% define the dimensions in mm
gridsizeInmm = [40, 25, 25]; 

cylinder2_C0Inmm = [15, 12, 17.5];
cylinder2_C1Inmm = [40, 12, 17.5];
cylinder2_RInmm  = 5;
cylinder3_C0Inmm = [0,  12, 17.5];
cylinder3_C1Inmm = [15, 12, 17.5];
cylinder3_RInmm  = 5;
cylinder4_C0Inmm = [0,  12, 17.5];
cylinder4_C1Inmm = [40, 12, 17.5];
cylinder4_RInmm  = 1;

% convert from mm to grid

gridSize = gridsizeInmm / unitinmm;
gridSize = floor(gridSize);
zLayers = [1, gridSize(3), 1];
cylinder2_C0 = (cylinder2_C0Inmm/unitinmm);
cylinder2_C1 = (cylinder2_C1Inmm/unitinmm);
cylinder2_R = (cylinder2_RInmm/unitinmm);
cylinder3_C0 = (cylinder3_C0Inmm/unitinmm);
cylinder3_C1 = (cylinder3_C1Inmm/unitinmm);
cylinder3_R = (cylinder3_RInmm/unitinmm);
cylinder4_C0 = (cylinder4_C0Inmm/unitinmm);
cylinder4_C1 = (cylinder4_C1Inmm/unitinmm);
cylinder4_R = (cylinder4_RInmm/unitinmm);


template = ['{"Shapes":[{"Grid":{"Tag": 0,"Size": GRID_SIZE_PLACEHOLDER}},' ...
              '{"ZLayers":[Z_LAYERS_PLACEHOLDER]},' ...
              '{"Cylinder":{"Tag":2,"C0":CYLINDER2_C0_PLACEHOLDER,"C1":CYLINDER2_C1_PLACEHOLDER,"R":CYLINDER2_R_PLACEHOLDER}},' ...
              '{"Cylinder":{"Tag":3,"C0":CYLINDER3_C0_PLACEHOLDER,"C1":CYLINDER3_C1_PLACEHOLDER,"R":CYLINDER3_R_PLACEHOLDER}},' ...
              '{"Cylinder":{"Tag":4,"C0":CYLINDER4_C0_PLACEHOLDER,"C1":CYLINDER4_C1_PLACEHOLDER,"R":CYLINDER4_R_PLACEHOLDER}}]}'];

template = strrep(template, 'GRID_SIZE_PLACEHOLDER', regexprep(mat2str(gridSize),'\s*',','));
template = strrep(template, 'Z_LAYERS_PLACEHOLDER', regexprep(mat2str(zLayers),'\s*',','));

template = strrep(template, 'CYLINDER2_C0_PLACEHOLDER', regexprep(mat2str(cylinder2_C0),'\s*',','));
template = strrep(template, 'CYLINDER2_C1_PLACEHOLDER', regexprep(mat2str(cylinder2_C1),'\s*',','));
template = strrep(template, 'CYLINDER2_R_PLACEHOLDER', regexprep(mat2str(cylinder2_R),'\s*',','));

template = strrep(template, 'CYLINDER3_C0_PLACEHOLDER', regexprep(mat2str(cylinder3_C0),'\s*',','));
template = strrep(template, 'CYLINDER3_C1_PLACEHOLDER', regexprep(mat2str(cylinder3_C1),'\s*',','));
template = strrep(template, 'CYLINDER3_R_PLACEHOLDER', regexprep(mat2str(cylinder3_R),'\s*',','));

template = strrep(template, 'CYLINDER4_C0_PLACEHOLDER', regexprep(mat2str(cylinder4_C0),'\s*',','));
template = strrep(template, 'CYLINDER4_C1_PLACEHOLDER', regexprep(mat2str(cylinder4_C1),'\s*',','));
template = strrep(template, 'CYLINDER4_R_PLACEHOLDER', regexprep(mat2str(cylinder4_R),'\s*',','));


cfg.shapes = template;

%mcxpreview(cfg)
% Define a new variable to change the location of the cylinder
Z_CInmm = [19.5 17.5 15.5];
%Z_CInmm = 12;
Z_C = Z_CInmm/unitinmm;

OutputData = struct();

for k = 1:1:length(Z_C)


shapesdata = jsondecode(cfg.shapes);
shapesdata.Shapes{3}.Cylinder.C0 = [cylinder2_C0(1),cylinder2_C0(2),Z_C(k)];
shapesdata.Shapes{3}.Cylinder.C1 = [cylinder2_C1(1),cylinder2_C1(2),Z_C(k)];
shapesdata.Shapes{4}.Cylinder.C0 = [cylinder3_C0(1),cylinder3_C0(2),Z_C(k)];
shapesdata.Shapes{4}.Cylinder.C1 = [cylinder3_C1(1),cylinder3_C1(2),Z_C(k)];
shapesdata.Shapes{5}.Cylinder.C0 = [cylinder4_C0(1),cylinder4_C0(2),Z_C(k)];
shapesdata.Shapes{5}.Cylinder.C1 = [cylinder4_C1(1),cylinder4_C1(2),Z_C(k)];
%shapesdata.Shapes{3}.Cylinder.R = 15;
cfg.shapes = jsonencode(shapesdata);

%%export from json file
%Convert JSON data to the cfg structure for mcxlab using json2mcx
% addpath('/home/shossei3/Desktop/UR1/iso2mesh');
% cfg = json2mcx('skinvessel_bone2.json');

%neccessary if want to use json2mcx
cfg.autopilot=1;

%% prepare cfg for MCX simulation
%clear cfg
%%secify number of photons 
cfg.nphoton=1e6;
%%specify excitation wavelength
lambda_exc = 785;

% tissue labels:0-ambient air,1-skin,2 bone
% first is vaccumm as default

%parameters from demo_head example 1- skin 2- diaphysis
%cfg.prop=[0 0 1 1;0.9172 3.8127 0.8900 1.3700;1.9062 1.6579 0.8900 1.3700]; %[mua, mus, g, n];

%Parameters from Sandell review paper 2011; 1- skin 2- bone diaphysis
mua_bone = mean(0.07:0.01:0.09)/10;         % devided by 10 to convert from 1/cm to 1/mm
mus_prime_bone = mean(11.9:0.1:14.1)/10;   % 1/mm
g = 0.89;
mus_bone = mus_prime_bone./(1-g);

mua_skin = mean(0.16:0.01:0.23)/10;
mus_prime_skin = mean(6.80:0.1:9.84)/10;
mus_skin = mus_prime_skin./(1-g);


%[mua, mus, g, n];  
cfg.prop=[0 0 1 1;
          mua_skin mus_skin 0.89 1.37;
          mua_bone mus_bone 0.89 1.37;
          mua_bone mus_bone 0.89 1.37;
          mua_skin mus_skin 0.89 1.37];

% 100% bone
% cfg.prop=[0 0 1 1;
%           mua_bone mus_bone 0.89 1.37;
%           mua_bone mus_bone 0.89 1.37;
%           mua_bone mus_bone 0.89 1.37;
%           mua_bone mus_bone 0.89 1.37];

% 100% skin
% cfg.prop=[0 0 1 1;
%           mua_skin mus_skin 0.89 1.37;
%           mua_skin mus_skin 0.89 1.37;
%           mua_skin mus_skin 0.89 1.37;
%           mua_skin mus_skin 0.89 1.37];


% light source
%cfg.srcpos=[size(cfg.vol,1)/2,size(cfg.vol,2)/2,size(cfg.vol,3)]; %inward-pointing source
srcposInmm = [20,12.5,25];
cfg.srcpos=srcposInmm/unitinmm;
cfg.srcdir=[0,0,-1]; 
cfg.srctype = 'disk';
srcdiameterInmm = 1;
srcdiameter = srcdiameterInmm/unitinmm;
srcradious = srcdiameter/2;
cfg.srcparam1 = [srcradious,0,0,0];

% time windows
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=5e-9;

% other simulation parameters
cfg.isspecular=0;
cfg.isreflect=1;
cfg.autopilot=1;
cfg.gpuid=1;

%% run MCX simulation for excitation
%addpath('/home/shossei3/Desktop/rmcx/2020/b2/src/mcxlab');
%addpath('/software/rmcx/2023.08.02/b1/src/mcx/mcxlab');
% addpath('/home/shossei3/Desktop/UR2/mcx/mcxlab');
% addpath('/home/shossei3/Desktop/UR1/mcx/mcxtools'); % should be investigat
% addpath('/home/shossei3/Desktop/UR1/iso2mesh');

cfg.detpos = [00 00 00 0];
cfg = rmfield(cfg,'detpos'); 

%cfg.issaveexit = 1;
%mcxpreview(cfg)
%cfg.vol = uint8(zeros((gridSize)));
[flux_exc,~,vol]=mcxlab(cfg);
cfg.vol =uint8(vol.data);


%%detector position and size
%cfg.detpos = [cfg.srcpos 10]; %%make detector a 10 mm diameter (???This is radious not diamter, also unitimm shoul be considered here) incident on launch position
%cfg.detpos = [30 30 60 10]; 

% cfg.detpos = [00 30 60 2;
%               10 30 60 2;
%               20 30 60 2;
%               30 30 60 2;
%               40 30 60 2;
%               50 30 60 2;
%               60 30 60 2]; 

% cfg.detpos = [080 50 100 2;
%               090 50 100 2;
%               100 50 100 2]; 

%cfg.detpos = [80 50 0 40]; 

%location of 3 detectors: 1- 0mm-1.2 mm diameter, 2- (68-80)*1/4=-3mm, 1.8mm
% diameter, 3- (56-80)*1/4=6mm, 2.4 mm diameter


% real location and real diameters 
% [Location: -0mm-->(80-80)*1/4=-0mm, 300um diamter-->300um*4=1.2mm diameter-->0.6unitmm radious]
% [Location: -3mm-->(68-80)*1/4=-3mm, 450um diamter-->450um*4=1.8mm diameter-->0.9unitmm radious]
% [Location: -6mm-->(56-80)*1/4=-6mm, 600um diamter-->600um*4=2.4mm diameter-->1.2unitmm radious]


% detZposInmm = 25+(unitinmm);
% detposInmm = [20 12.5 detZposInmm 0.150;
%               17 12.5 detZposInmm 0.225;
%               14 12.5 detZposInmm 0.300];

cfg.detpos = [20/unitinmm 12.5/unitinmm  floor((25+unitinmm)/unitinmm) 0.150/unitinmm;
              17/unitinmm 12.5/unitinmm  floor((25+unitinmm)/unitinmm) 0.225/unitinmm; 
              14/unitinmm 12.5/unitinmm  floor((25+unitinmm)/unitinmm) 0.300/unitinmm;];

%cfg.detpos = (detposInmm/unitinmm);
% cfg.detpos = [080 50 100 1.1;
%               068 50 100 1.65;
%               056 50 100 2.20]; 

% cfg.detpos = [76 48 96 0.5769;
%               65.38 48.769  96.1538 0.8654;
%               53.8462 48.0769 96.1538 1.1538];

% diameters are greater than real to prevent zero detected photons
% cfg.detpos = [080 50 100 10;
%               070 50 100 10;
%               060 50 100 10]; 

%cfg.detpos = [30 30 60 2];


%mcxpreview(cfg)
%% show representative image of excitation flux
%close all

%mcxplotvol(log10(flux_exc.data));
%mcxplotvol((vol.data));
%figure(3);
%mcxpreview(cfg)
%imagesc(rot90(squeeze(log10(flux_exc.data(80,:,:)))))
%axis equal

%% begin setup for emmission launches
cfg.srctype='weighed';
cfg.srcpos=[1 1 1];
%cfg.srcparam1=[size(cfg.vol,1) size(cfg.vol,2) size(cfg.vol,3)];
cfg.srcparam1=gridSize;

%%calculate how much of the fluence is in each label

Skin = cfg.vol == 1;
CDiaphysis = cfg.vol == 2;
CEpiphysis = cfg.vol == 3;
MCanal = cfg.vol == 4;

%%also get total flux through the system in terms of labels 
Total = sum(flux_exc.data(:));

nphotonsSkin = Skin.*flux_exc.data;
nphotonsSkin = round(double(cfg.nphoton*sum(nphotonsSkin(:))/Total));

nphotonsCDiaphysis = CDiaphysis.*flux_exc.data;
nphotonsCDiaphysis = round(double(cfg.nphoton*sum(nphotonsCDiaphysis(:))/Total));

nphotonsCEpiphysis = CEpiphysis.*flux_exc.data;
nphotonsCEpiphysis = round(double(cfg.nphoton*sum(nphotonsCEpiphysis(:))/Total));

nphotonsMCanal = MCanal.*flux_exc.data;
nphotonsMCanal = round(double(cfg.nphoton*sum(nphotonsMCanal(:))/Total));

%nphotonsAir = Air.*flux_exc.data;
%nphotonsAir = round(double(cfg.nphoton*sum(nphotonsAir(:))/Total));

%%time-gate for emission launches 
cfg.tstart=0;
cfg.tend=1e-9; 
cfg.tstep=1e-9;


%%set srcpattern to results from flux_exc
cfg.srcpattern = double(flux_exc.data);


load ThreeSpectra.mat
% plot(wavenum,AnitaSpec2)
% xlabel('wavenum')
% ylabel('Peak intensities')
% legend('Diaphysis','epiphysis','Skin');

% wavenum = 920:2:1020;
% wavenum = wavenum';
RamanShift = wavenum;


%smh: single wave simulation to calculate photon numbers for each detector
%RamanShift = RamanShift(1:2);

% to reduce time of simulation;
%RamanShift = RamanShift(1:end);

em_wave = (-(RamanShift/1e7 - 1/785)).^-1;

%addpath('/home/shossei3/Desktop/UR2/Raman-Monte-Carlo-Extreme-master');
parfor_progress(length(RamanShift));
for ind = 1:length(em_wave)
    cfglocal = cfg;
    
    lambda_em = em_wave(ind);
    % tissue labels:0-ambient air,1-scalp,2-skull,3-csf,4-gray matter,5-white matter,6-air cavities
%     cfglocal.prop=[0,0,1,1;
%         mua_lambda(lambda_em,.75,.1,.1,.05,0,.1,.01) mus_lambda(lambda_exc,46.1,.421) 0.89 1.37;
%         mua_lambda(lambda_em,.9,.2,.1,.05,0.01,.1,.0) mus_lambda(lambda_em,22.9,.716) 0.89 1.37;
%         0.001,0.001,0.89,1.37;
%         mua_lambda(lambda_em,.9,.5,.1,0,0,.5,.05) mus_lambda(lambda_em,24.2,1.611) 0.89 1.37;
%         mua_lambda(lambda_em,.92,.6,.1,0,0,.4,.05) mus_lambda(lambda_em,22,1.4) 0.84 1.37;
%         0,0,1,1];
    %cfg.prop=[0 0 1 1;0.5 1 0 1.37;0.2 10 0.9 1.37];
    %cfg.prop=[0 0 1 1;0.9172 3.8127 0.8900 1.3700;1.9062 1.6579 0.8900 1.3700];

    %%launch emmission bounded to launch in just the scalp, with
    %%appropriate number of photons as determined from the excitation
    %%matrix
    if(nphotonsSkin>0)
        cfglocal.nphoton = nphotonsSkin;
        cfglocal.srcpattern = cfg.srcpattern.*Skin;
        %cfglocal.srcpattern = cfg.srcpattern;
        %addpath('/home/shossei3/Desktop/rmcx/2020/b2/src/mcxlab');
        %mcxplotvol((vol.data));
        cfglocal.issaveexit = 1;
        [~,det]=mcxlab(cfglocal);
        det2 = det;
        det.data = det.data(1:6,:);
        AA = size(det.data);


        detIDs = 1:size(cfg.detpos,1);
        DetPhotonosIds = det.data(1,:);
        DetPhotonosNSc = det.data(2,:);
        DetPhotonosPM1 = det.data(3,:);
        DetPhotonosPM2 = det.data(4,:);
        DetPhotonosPM3 = det.data(5,:);
        DetPhotonosPM4 = det.data(6,:);

        Det1Indices = DetPhotonosIds == 1;
        Det1Data = [DetPhotonosIds(Det1Indices);
                        DetPhotonosNSc(Det1Indices);
                        DetPhotonosPM1(Det1Indices);
                        DetPhotonosPM2(Det1Indices);
                        DetPhotonosPM3(Det1Indices);
                        DetPhotonosPM4(Det1Indices)];

        Det2Indices = DetPhotonosIds == 2;
        Det2Data = [DetPhotonosIds(Det2Indices);
                        DetPhotonosNSc(Det2Indices);
                        DetPhotonosPM1(Det2Indices);
                        DetPhotonosPM2(Det2Indices);
                        DetPhotonosPM3(Det2Indices);
                        DetPhotonosPM4(Det2Indices)];

        Det3Indices = DetPhotonosIds == 3;
        Det3Data = [DetPhotonosIds(Det3Indices);
                        DetPhotonosNSc(Det3Indices);
                        DetPhotonosPM1(Det3Indices);
                        DetPhotonosPM2(Det3Indices);
                        DetPhotonosPM3(Det3Indices);
                        DetPhotonosPM4(Det3Indices)];

        %%calculate reflectance from the each detected photons 
        weight1 = length(Det1Data)./cfg.nphoton;
        %
        RSkin1(ind) = sum(mean((weight1*exp(-cfglocal.prop(2:end,1).*Det1Data(3:end,:))),2));
        %%multiply reflectance by raman feature for the sckin
        RamanSkin1(ind) = interp1(RamanShift,AnitaSpec2(1:size(RamanShift,1),3),RamanShift(ind))*RSkin1(ind);

        %%calculate reflectance from the each detected photons 
        weight2 = length(Det2Data)./cfg.nphoton;
        %
        RSkin2(ind) = sum(mean((weight2*exp(-cfglocal.prop(2:end,1).*Det2Data(3:end,:))),2));
        %%multiply reflectance by raman feature for the sckin
        RamanSkin2(ind) = interp1(RamanShift,AnitaSpec2(1:size(RamanShift,1),3),RamanShift(ind))*RSkin2(ind);

        %%calculate reflectance from the each detected photons 
        weight3 = length(Det3Data)./cfg.nphoton;
        %
        RSkin3(ind) = sum(mean((weight3*exp(-cfglocal.prop(2:end,1).*Det3Data(3:end,:))),2));
        %%multiply reflectance by raman feature for the sckin
        RamanSkin3(ind) = interp1(RamanShift,AnitaSpec2(1:size(RamanShift,1),3),RamanShift(ind))*RSkin3(ind);


        %Total
        %%calculate reflectance from the detected photons 
        weight = length(det.data)./cfg.nphoton;
        %
        RSkin(ind) = sum(mean((weight*exp(-cfglocal.prop(2:end,1).*det.data(3:end,:))),2));
        %%multiply reflectance by raman feature for the scalp
        RamanSkin(ind) = interp1(RamanShift,AnitaSpec2(1:size(RamanShift,1),3),RamanShift(ind))*RSkin(ind);

        % at least 2 values shoule be selected for interp1
        %RamanSkin(k) = interp1(RamanShift,AnitaSpec2(1:2,3),RamanShift(ind))*RSkin(ind);

    end
    if(nphotonsCDiaphysis)
        cfglocal.nphoton = nphotonsCDiaphysis;
        cfglocal.srcpattern = cfg.srcpattern.*CDiaphysis;
        [~,det]=mcxlab(cfglocal);
        det3 = det;
        det.data = det.data(1:6,:);

        detIDs = 1:size(cfg.detpos,1);
        DetPhotonosIds = det.data(1,:);
        DetPhotonosNSc = det.data(2,:);
        DetPhotonosPM1 = det.data(3,:);
        DetPhotonosPM2 = det.data(4,:);
        DetPhotonosPM3 = det.data(5,:);
        DetPhotonosPM4 = det.data(6,:);

        Det1Indices = DetPhotonosIds == 1;
        Det1Data = [DetPhotonosIds(Det1Indices);
                        DetPhotonosNSc(Det1Indices);
                        DetPhotonosPM1(Det1Indices);
                        DetPhotonosPM2(Det1Indices);
                        DetPhotonosPM3(Det1Indices);
                        DetPhotonosPM4(Det1Indices)];

        Det2Indices = DetPhotonosIds == 2;
        Det2Data = [DetPhotonosIds(Det2Indices);
                        DetPhotonosNSc(Det2Indices);
                        DetPhotonosPM1(Det2Indices);
                        DetPhotonosPM2(Det2Indices);
                        DetPhotonosPM3(Det2Indices);
                        DetPhotonosPM4(Det2Indices)];

        Det3Indices = DetPhotonosIds == 3;
        Det3Data = [DetPhotonosIds(Det3Indices);
                        DetPhotonosNSc(Det3Indices);
                        DetPhotonosPM1(Det3Indices);
                        DetPhotonosPM2(Det3Indices);
                        DetPhotonosPM3(Det3Indices);
                        DetPhotonosPM4(Det3Indices)];

        %%calculate reflectance from the each detected photons 
        weight1 = length(Det1Data)./cfg.nphoton;
        %
        if (weight1 ~=0)
        RCDia1(ind) = sum(mean((weight1*exp(-cfglocal.prop(2:end,1).*Det1Data(3:end,:))),2));
        %%multiply reflectance by raman feature for the sckin
        RamanCDia1(ind) = interp1(RamanShift,AnitaSpec2(1:size(RamanShift,1),1),RamanShift(ind))*RCDia1(ind);
        end

        if (weight1 ==0)
            RCDia1(ind) = 0;
            RamanCDia1(ind) = 0;
        end

        %%calculate reflectance from the each detected photons 
        weight2 = length(Det2Data)./cfg.nphoton;
        %
        if (weight2 ~=0)
        RCDia2(ind) = sum(mean((weight2*exp(-cfglocal.prop(2:end,1).*Det2Data(3:end,:))),2));
        %%multiply reflectance by raman feature for the sckin
        RamanCDia2(ind) = interp1(RamanShift,AnitaSpec2(1:size(RamanShift,1),1),RamanShift(ind))*RCDia2(ind);
        end

        if (weight2 ==0)
            RCDia2(ind) = 0;
            RamanCDia2(ind) = 0;
        end

        %%calculate reflectance from the each detected photons 
        weight3 = length(Det3Data)./cfg.nphoton;
        %
        if (weight3 ~=0)
        RCDia3(ind) = sum(mean((weight3*exp(-cfglocal.prop(2:end,1).*Det3Data(3:end,:))),2));
        %%multiply reflectance by raman feature for the sckin
        RamanCDia3(ind) = interp1(RamanShift,AnitaSpec2(1:size(RamanShift,1),1),RamanShift(ind))*RCDia3(ind);
        end

         if (weight3 ==0)
            RCDia3(ind) = 0;
            RamanCDia3(ind) = 0;
        end       

        %Totoal
        %%calculate reflectance from the detected photons 
        weight = length(det.data)./cfg.nphoton;
        RCDiaphysis(ind) = sum(mean((weight*exp(-cfglocal.prop(2:end,1).*det.data(3:end,:))),2));
        RamanCDiaphysis(ind) = interp1(RamanShift,AnitaSpec2(1:size(RamanShift,1),1),RamanShift(ind))*RCDiaphysis(ind);

        % at least 2 values shoule be selected for interp1
        %RamanDiaphysis(k) = interp1(RamanShift,AnitaSpec2(1:2,1),RamanShift(ind))*RDiaphysis(ind);
    end
    if(nphotonsCEpiphysis)
        cfglocal.nphoton = nphotonsCEpiphysis;
        cfglocal.srcpattern = cfg.srcpattern.*CEpiphysis;
        [~,det]=mcxlab(cfglocal);
        det4 = det.data;
        det.data = det.data(1:6,:);

        detIDs = 1:size(cfg.detpos,1);
        DetPhotonosIds = det.data(1,:);
        DetPhotonosNSc = det.data(2,:);
        DetPhotonosPM1 = det.data(3,:);
        DetPhotonosPM2 = det.data(4,:);
        DetPhotonosPM3 = det.data(5,:);
        DetPhotonosPM4 = det.data(6,:);

        Det1Indices = DetPhotonosIds == 1;
        Det1Data = [DetPhotonosIds(Det1Indices);
                        DetPhotonosNSc(Det1Indices);
                        DetPhotonosPM1(Det1Indices);
                        DetPhotonosPM2(Det1Indices);
                        DetPhotonosPM3(Det1Indices);
                        DetPhotonosPM4(Det1Indices)];

        Det2Indices = DetPhotonosIds == 2;
        Det2Data = [DetPhotonosIds(Det2Indices);
                        DetPhotonosNSc(Det2Indices);
                        DetPhotonosPM1(Det2Indices);
                        DetPhotonosPM2(Det2Indices);
                        DetPhotonosPM3(Det2Indices);
                        DetPhotonosPM4(Det2Indices)];

        Det3Indices = DetPhotonosIds == 3;
        Det3Data = [DetPhotonosIds(Det3Indices);
                        DetPhotonosNSc(Det3Indices);
                        DetPhotonosPM1(Det3Indices);
                        DetPhotonosPM2(Det3Indices);
                        DetPhotonosPM3(Det3Indices);
                        DetPhotonosPM4(Det3Indices)];


        %%calculate reflectance from the each detected photons 
        weight1 = length(Det1Data)./cfg.nphoton;
        %
        if (weight1 ~=0)
        RCEpi1(ind) = sum(mean((weight1*exp(-cfglocal.prop(2:end,1).*Det1Data(3:end,:))),2));
        %%multiply reflectance by raman feature for the sckin
        RamanCEpi1(ind) = interp1(RamanShift,AnitaSpec2(1:size(RamanShift,1),2),RamanShift(ind))*RCEpi1(ind);
        end

        if (weight1 ==0)
            RCEpi1(ind) = 0;
            RamanCEpi1(ind) = 0;
        end

        %%calculate reflectance from the each detected photons 
        weight2 = length(Det2Data)./cfg.nphoton;
        %
        if (weight2 ==0)
            RCEpi2(ind) = 0;
            RamanCEpi2(ind) = 0;
        end
        if (weight2 ~=0)
        RCEpi2(ind) = sum(mean((weight2*exp(-cfglocal.prop(2:end,1).*Det2Data(3:end,:))),2));
        %%multiply reflectance by raman feature for the sckin
        RamanCEpi2(ind) = interp1(RamanShift,AnitaSpec2(1:size(RamanShift,1),2),RamanShift(ind))*RCEpi2(ind);
        end

        %%calculate reflectance from the each detected photons 
        weight3 = length(Det3Data)./cfg.nphoton;
        %
        if (weight3 ==0)
            RCEpi3(ind) = 0;
            RamanCEpi3(ind) = 0;
        end
        if (weight3 ~=0)
        RCEpi3(ind) = sum(mean((weight3*exp(-cfglocal.prop(2:end,1).*Det3Data(3:end,:))),2));
        %%multiply reflectance by raman feature for the sckin
        RamanCEpi3(ind) = interp1(RamanShift,AnitaSpec2(1:size(RamanShift,1),2),RamanShift(ind))*RCEpi3(ind);
        end

        %Totoal
        %%calculate reflectance from the detected photons 
        weight = length(det.data)./cfg.nphoton;
        RCEpiphysis(ind) = sum(mean((weight*exp(-cfglocal.prop(2:end,1).*det.data(3:end,:))),2));
        RamanCEpiphysis(ind) = interp1(RamanShift,AnitaSpec2(1:size(RamanShift,1),2),RamanShift(ind))*RCEpiphysis(ind);

        % at least 2 values shoule be selected for interp1
        %RamanDiaphysis(k) = interp1(RamanShift,AnitaSpec2(1:2,1),RamanShift(ind))*RDiaphysis(ind);
    end
    if(nphotonsMCanal)
        cfglocal.nphoton = nphotonsMCanal;
        cfglocal.srcpattern = cfg.srcpattern.*MCanal;
        [~,det]=mcxlab(cfglocal);
        det5 = det.data;
        det.data = det.data(1:6,:);

        detIDs = 1:size(cfg.detpos,1);
        DetPhotonosIds = det.data(1,:);
        DetPhotonosNSc = det.data(2,:);
        DetPhotonosPM1 = det.data(3,:);
        DetPhotonosPM2 = det.data(4,:);
        DetPhotonosPM3 = det.data(5,:);
        DetPhotonosPM4 = det.data(6,:);

        Det1Indices = DetPhotonosIds == 1;
        Det1Data = [DetPhotonosIds(Det1Indices);
                        DetPhotonosNSc(Det1Indices);
                        DetPhotonosPM1(Det1Indices);
                        DetPhotonosPM2(Det1Indices);
                        DetPhotonosPM3(Det1Indices);
                        DetPhotonosPM4(Det1Indices)];

        Det2Indices = DetPhotonosIds == 2;
        Det2Data = [DetPhotonosIds(Det2Indices);
                        DetPhotonosNSc(Det2Indices);
                        DetPhotonosPM1(Det2Indices);
                        DetPhotonosPM2(Det2Indices);
                        DetPhotonosPM3(Det2Indices);
                        DetPhotonosPM4(Det2Indices)];

        Det3Indices = DetPhotonosIds == 3;
        Det3Data = [DetPhotonosIds(Det3Indices);
                        DetPhotonosNSc(Det3Indices);
                        DetPhotonosPM1(Det3Indices);
                        DetPhotonosPM2(Det3Indices);
                        DetPhotonosPM3(Det3Indices);
                        DetPhotonosPM4(Det3Indices)];


        medfactor = 0.05;


        %%calculate reflectance from the each detected photons 
        weight1 = length(Det1Data)./cfg.nphoton;
        %
        RMCanal1(ind) = sum(mean((weight1*exp(-cfglocal.prop(2:end,1).*Det1Data(3:end,:))),2));
        %%multiply reflectance by raman feature for the sckin
        RamanMCanal1(ind) = interp1(RamanShift,medfactor*AnitaSpec2(1:size(RamanShift,1),3),RamanShift(ind))*RMCanal1(ind);

        %%calculate reflectance from the each detected photons 
        weight2 = length(Det2Data)./cfg.nphoton;
        %
        RMCanal2(ind) = sum(mean((weight2*exp(-cfglocal.prop(2:end,1).*Det2Data(3:end,:))),2));
        %%multiply reflectance by raman feature for the sckin
        RamanMCanal2(ind) = interp1(RamanShift,medfactor*AnitaSpec2(1:size(RamanShift,1),3),RamanShift(ind))*RMCanal2(ind);

        %%calculate reflectance from the each detected photons 
        weight3 = length(Det3Data)./cfg.nphoton;
        %
        RMCanal3(ind) = sum(mean((weight3*exp(-cfglocal.prop(2:end,1).*Det3Data(3:end,:))),2));
        %%multiply reflectance by raman feature for the sckin
        RamanMCanal3(ind) = interp1(RamanShift,medfactor*AnitaSpec2(1:size(RamanShift,1),3),RamanShift(ind))*RMCanal3(ind);


        %Totoal
        %%calculate reflectance from the detected photons 
        weight = length(det.data)./cfg.nphoton;
        RMCanal(ind) = sum(mean((weight*exp(-cfglocal.prop(2:end,1).*det.data(3:end,:))),2));
        RamanMCanal(ind) = interp1(RamanShift,medfactor*AnitaSpec2(1:size(RamanShift,1),3),RamanShift(ind))*RMCanal(ind);

        % at least 2 values shoule be selected for interp1
        %RamanDiaphysis(k) = interp1(RamanShift,AnitaSpec2(1:2,1),RamanShift(ind))*RDiaphysis(ind);
    end

    parfor_progress;
    
    %calculating ration of Bone to Skin
    %R(k) = RamanDiaphysis(k)/RamanSkin(k);
    %R2 = photonsCounts1_bone./photonsCounts1_skin;
end


RamanTotalDet1 = RamanSkin1+RamanCDia1+RamanCEpi1+RamanMCanal1;
RamanTotalDet2 = RamanSkin2+RamanCDia2+RamanCEpi2+RamanMCanal2;
RamanTotalDet3 = RamanSkin3+RamanCDia3+RamanCEpi3+RamanMCanal3;

RamanAllDets = RamanSkin+RamanCDiaphysis+RamanCEpiphysis+RamanMCanal;
MTotal = RamanAllDets;


S = AnitaSpec2(1:size(RamanShift,1),3);
B_CD = AnitaSpec2(1:size(RamanShift,1),1);
B_CE = AnitaSpec2(1:size(RamanShift,1),2);
B_MC = AnitaSpec2(1:size(RamanShift,1),3);

M1 = RamanTotalDet1;
coefficients1 = [S,B_CD,B_CE,B_MC]\M1';
a1S = coefficients1(1);
a1B_CD = coefficients1(2);
a1B_CE = coefficients1(3);
a1B_MC = coefficients1(4);


M2 = RamanTotalDet2;
coefficients2 = [S,B_CD,B_CE,B_MC]\M2';
a2S = coefficients2(1);
a2B_CD = coefficients2(2);
a2B_CE = coefficients2(3);
a2B_MC = coefficients2(4);

M3 = RamanTotalDet3;
coefficients3 = [S,B_CD,B_CE,B_MC]\M3';
a3S = coefficients3(1);
a3B_CD = coefficients3(2);
a3B_CE = coefficients3(3);
a3B_MC = coefficients3(4);

M1Prime = a1S*S+a1B_CD*B_CD+a1B_CE*B_CE+a1B_MC*B_MC;
M2Prime = a2S*S+a2B_CD*B_CD+a2B_CE*B_CE+a2B_MC*B_MC;
M3Prime = a3S*S+a3B_CD*B_CD+a3B_CE*B_CE+a3B_MC*B_MC;

OutputData(k).RamanShift = RamanShift;
OutputData(k).M1 = M1;
OutputData(k).M2 = M2;
OutputData(k).M3 = M3;
OutputData(k).M1Prime = M1Prime;
OutputData(k).M2Prime = M2Prime;
OutputData(k).M3Prime = M3Prime;
OutputData(k).MTotal = MTotal;


end

% Save the entire structure array to a single file
save('OutputData1mmMcanal.mat', 'OutputData');



% figure 
% plot (Z_C,R);
% xlabel('Cylinder location')
% ylabel('Ratio = Raman(PO4)/Raman(Protein)')







% figure
% plot(RamanShift,M1);
% hold on;
% plot(RamanShift,M1Prime);
% hold off;
% 
% figure
% plot(RamanShift,M2);
% hold on;
% plot(RamanShift,M2Prime);
% hold off;
% 
% figure
% plot(RamanShift,M3);
% hold on;
% plot(RamanShift,M3Prime);
% hold off;


% figure
% plot(RamanShift,RamanSkin)
% hold on
% plot(RamanShift,RamanDiaphysis)

%plot(RamanShift,RamanGrey)
%plot(RamanShift,RamanGrey)
%plot(RamanShift,RamanWhite)
%plot(RamanShift,RamanAir)
% xlabel('Raman Shift')
% ylabel('Intensity (a.u.)')

% fname = '/home/shossei3/Desktop/UR1';
% saveas(gcf,fullfile(fname,'fig10'),'png');