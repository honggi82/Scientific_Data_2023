% It is sample MATLAB codes for the manuscript entitled 
% "A magnetoencephalography dataset during three-dimensional reaching movements for brain-computer interfaces". 
% A "topoplot" function and related sub-functions in the EEGLAB toolbox are needed to plot a topography. 
% The EEGLAB software is available at https://sccn.ucsd.edu/eeglab30. 
% For time-frequency analysis or FTF analysis, our custom MATLAB codes are required. 
% It can be found at https://github.com/honggi82/Scientific_Data_2023. 
% To use these toolboxes, the folders containing the toolboxes should be added to the "Set Path" of MATLAB. 
% The "Out of memory" error could occur when MEG data is processed because the data size is quite large. 
% To this end, we recommend using a computer with a large amount of RAM and using 64-bit MATLAB and a 64-bit operating system (OS). 
% Copyright (c) 2023 Hong Gi Yeom

clear all; close all; clc;

% adding the data folder to setpath
addpath(genpath('C:\SciData')) % The location 'C:\SciData' should be changed
addpath(genpath('E:\MEG센터 자료\matlab\down\eeglab2021.0')) 

%====================== 1. Parameter setting ======================
sf = 600.615;          % sampling frequency
ch_n = 306;            % channel number (102 for mag, 204 for grad, 306 for all)
wnd_size=[-1 2];       % time duration for analysis
baseline=[-1 0];       % time duration for baseline
f_scale=1;             % frequency resolution for time-frequency analysis
freq_band=[0.1 100];   % frequency band for time-frequency analysis
normal=1;              % normalization or not
scr_sz=get(0,'ScreenSize'); % screen size

% channel location
switch ch_n            
    case 102
        [position]=neuromag_layout_mag;
    case 204
        [position]=neuromag_layout_grad;
    case 306
        [position]=neuromag_layout_all;
end

% iteration for all subject and session
for sub=1:9       % subject number
    for ses=1:2   % session number

        %====================== 2. Loading data ======================
        % file name
        switch ch_n            
            case 102
                file_n=['Sub_',num2str(sub),'_ses_',num2str(ses), '_epoched_mag.mat'];
            case 204
                file_n=['Sub_',num2str(sub),'_ses_',num2str(ses), '_epoched_grad.mat'];
            case 306
                file_n=['Sub_',num2str(sub),'_ses_',num2str(ses), '_epoched_all.mat'];
        end
        
        load(file_n); % loading data
        eve_n = length(epoched_data); % number of events
        
        %====================== 3. Reshaping the data ======================
        for i=1:eve_n
            data(:,:,:,i)=epoched_data{i}(1:ch_n,:,:);
            acc{i}=epoched_data{i}(317:319,:,:);
        end
        clear epoched_data
        sz = size(data); % data size

        %====================== 4. Plotting position ====================== 
        [velocity, position]=acc2pos(acc,sf);

        m_st = round(1.3*sf); % movement start
        m_end = round(2.3*sf); % movement end
        figure(); set(gcf,'Color','w'); hold on;
        for j=1:size(position{1},3);
            plot3(position{1}(1,m_st:m_end,j),position{1}(2,m_st:m_end,j),position{1}(3,m_st:m_end,j),'k'); 
            plot3(position{2}(1,m_st:m_end,j),position{2}(2,m_st:m_end,j),position{2}(3,m_st:m_end,j),'r'); 
            plot3(position{3}(1,m_st:m_end,j),position{3}(2,m_st:m_end,j),position{3}(3,m_st:m_end,j),'g'); 
            plot3(position{4}(1,m_st:m_end,j),position{4}(2,m_st:m_end,j),position{4}(3,m_st:m_end,j),'b');     
        end        

        %====================== 5. Calculation of time-freq power spectra ====================== 
        TF=zeros(sz(1), ceil(freq_band(2)-freq_band(1)/f_scale), sz(2), sz(3), sz(4));
        for ch=1:sz(1)
            for j=1:sz(3)
                for i=1:sz(4)
                    TF(ch,:,:,j,i) = timefreq_anal(data(ch,:,j,i), sf, wnd_size, baseline,f_scale,freq_band, normal);
                end
            end
            fprintf(['Calculation of time-freq anal on channel ', num2str(ch), ' is finished\n']);
        end
        clear data
        
        m_TF(:,:,:,(sub-1)*2+ses)=mean(mean(TF,5),4); % averaging by trials and event types
        %====================== 6. FTF analysis ================== 
        % We used the following code to prevent "Out of memory", 
        % although you can use "FTF_anal" function directly.

        temp1=mean(TF, 4);
        temp2=var(TF, 0, 4);
        clear TF

        ni = sz(3);     % # of trials of i-th group
        K = sz(4);      % # of group
        N = K*ni;       % # of total trials
        B = zeros(sz(1), ceil(freq_band(2)-freq_band(1)/f_scale), sz(2), sz(3), sz(4));
        W = zeros(sz(1), ceil(freq_band(2)-freq_band(1)/f_scale), sz(2), sz(3), sz(4));
        B = squeeze(var(temp1, 0, 5))/(K-1);  % Between variance
        W = squeeze(mean(temp2,5))/(N-K);     % Within variance
        F=B./W;                               % F-value

        file_n=['TimeFreq_Ch_',num2str(ch_n),'_sub_',num2str(sub),'_ses_',num2str(ses),'.mat'];
        save(file_n, 'F', 'm_TF');            % saving F-value
        clear B W F

    end
end

%====================== 7. Plotting time-freq power spectra ====================== 
AVG_TF=mean(m_TF,4);                                     % averaging by subjects

t=linspace(wnd_size(1),wnd_size(2),size(AVG_TF,3));      % time points
fr=linspace(freq_band(1), freq_band(2),size(AVG_TF,2));  % frequency points

% for whole channel plot
figure('Position',[0 0 scr_sz(3) scr_sz(4)]);
for ch=1:ch_n 
    subplot('Position', position(ch,:));
    pcolor(t,fr,squeeze(AVG_TF(ch,:,:))); 
    shading 'interp'; 
    colormap(jet);
    caxis([-0.5 0.5]);      
    x1=[0 0]; y1=freq_band;
    line(x1,y1,'Color','red', 'LineWidth', 1);
    set(gca,'xtick',[],'ytick',[]);
    ylim([0 60])
end
set(gcf,'Color','w')

% for 43-45 channel plot
figure('Position',[0 0 scr_sz(3)/4 scr_sz(4)/4]);
for ch=43:45
    switch ch
        case 43
            subplot('Position', [0.1 0.1 0.4 0.4]);
        case 44
            subplot('Position', [0.1 0.55 0.4 0.4]);
        case 45
            subplot('Position', [0.55 0.325 0.4 0.4]);
    end
    pcolor(t,fr,squeeze(AVG_TF(ch,:,:))); 
    shading 'interp'; 
    colormap(jet);
    caxis([-0.5 0.5]);      
    x1=[0 0]; y1=freq_band;
    line(x1,y1,'Color','red', 'LineWidth', 1);
    if ch==44||ch==45; set(gca,'xtick',[],'ytick',[]); 
    else; xticks([-1 0 1 2]); yticks([0 10 20 30 40 50 60]); fontsize(gca,24,"pixels"); end
    if ch==45; colorbar; fontsize(gca,24,"pixels"); end
    ylim([0 60])
end
set(gcf,'Color','w')

%====================== 8. Plotting topography of time-frequency spectra ====================== 
freq = [14 30]; % [1 8], [9 13]
time = [-0.1 1];
interval = 0.1;

topo_ch=2; % 1: Mag, 2: Grad
switch topo_ch
    case 1
        MEG=AVG_TF(3:3:306,:,:);
    case 2
        grad(:,:,:,1)=AVG_TF(1:3:306,:,:);
        grad(:,:,:,2)=AVG_TF(2:3:306,:,:);
        MEG=mean(grad,4); clear grad
end

switch freq(2)
    case 8
        c_axis=[-0.8 0.8];
    case 13
        c_axis=[-0.7 0.7];
    otherwise
        c_axis=[-0.5 0.5];
end

AVG_freq=squeeze(mean(MEG(:,freq,:),2));

figure('Position',[0 0 scr_sz(3) scr_sz(4)]); colormap(jet);
j=1;
for t_point=round((time(1)-wnd_size(1))*sf):round(interval*sf):round((time(2)-wnd_size(1))*sf)
    subplot(1,ceil((time(2)-time(1))/interval)+2, j);
    topoplot(AVG_freq(:,t_point),'sensor_102mag.locs','style','map','electrodes','off', 'maplimits', c_axis, 'whitebk','on', 'plotrad',0.548, 'headrad',0.5, 'shading','interp');
    title([num2str(round((t_point/sf)+wnd_size(1),2)), 's'])
    j=j+1;
end

%====================== 9. plotting F-value on 43-45 channels ====================== 
for sub=1:9
    for ses=1:2
        file_n=['TimeFreq_Ch_',num2str(ch_n),'_sub_',num2str(sub),'_ses_',num2str(ses),'.mat'];
        load(file_n);
        A_F(:,:,:,(sub-1)*2+ses)=F;
    end
end
AVG_F=mean(A_F,4);

figure('Position',[0 0 scr_sz(3)/4 scr_sz(4)/4]);
for ch=43:45
    switch ch
        case 43
            subplot('Position', [0.1 0.1 0.4 0.4]);
        case 44
            subplot('Position', [0.1 0.55 0.4 0.4]);
        case 45
            subplot('Position', [0.55 0.325 0.4 0.4]);
    end
    pcolor(t,fr,squeeze(AVG_F(ch,:,:))); 
    shading 'interp'; 
    colormap(parula);
    caxis([0 3]);      
    x1=[0 0]; y1=freq_band;
    line(x1,y1,'Color','red', 'LineWidth', 1);
    if ch==44||ch==45; set(gca,'xtick',[],'ytick',[]); 
    else; xticks([-1 0 1 2]); yticks([0 10 20 30 40 50 60]); fontsize(gca,24,"pixels"); end
    if ch==45; colorbar; fontsize(gca,24,"pixels"); end
    ylim([0 60])
end
set(gcf,'Color','w')