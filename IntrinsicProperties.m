%% Clean up
clear
close all
clc

imgshow = 0;
addpath(genpath('J:\Data\Matlab'));             
%% Load file
filepath = uigetdir;
cd(filepath);
        allDat = dir('*.dat');
        allData = {};
        for aa = 1:size(allDat,1)
            inputFile = HEKA_Importer(allDat(aa).name);
            dataPart = inputFile.RecTable.dataRaw;
            for dd = 1:size(dataPart,1)           
                allData{end+1,1} = dataPart{dd};
            end
        end
        
filepath

%% Define the parameters         (Check the inputFile.RecTable)
% exp = 13;        % define the session to be extracted
SR = inputFile.RecTable.SR;     % sampling rate
gain_sig = 10^-3;
gain_stim = 10^-12;
% ExperimentName: # slice
% Experiment: # cell
%% Process the data and Plot the data
list = [];
for exp = 1:size(inputFile.RecTable,1)
    k = size(allData{exp},2);
    for n = 1:k
        if inputFile.RecTable.Rec(exp) > 9
            a = char(strcat(inputFile.RecTable.ExperimentName(exp),'_',num2str(inputFile.RecTable.Experiment(exp)),'_',num2str(inputFile.RecTable.Rec(exp)),'_',inputFile.RecTable.Stimulus(exp),'_',num2str(n))); 
        else
            a = char(strcat(inputFile.RecTable.ExperimentName(exp),'_',num2str(inputFile.RecTable.Experiment(exp)),'_0',num2str(inputFile.RecTable.Rec(exp)),'_',inputFile.RecTable.Stimulus(exp),'_',num2str(n))); 
        end
   
        if isempty(inputFile.RecTable.ChUnit{exp}{1,n}) ~= 1   
            if inputFile.RecTable.ChUnit{exp}{1,n}(1,1) == 'V'
                allData{exp}{1,n} = allData{exp}{1,n}./gain_sig;
   
                if imgshow
                    figure('Name',a);
                    plot((1:length(inputFile.RecTable.dataRaw{exp,1}{1,n}))/SR(exp,1),allData{exp}{1,n})  % Recording signals  (mV/s)
                end
   
            else % == 'A'
                allData{exp}{1,n} = allData{exp}{1,n}./gain_stim;
   
                if imgshow
                    figure('Name',a);
                    plot((1:length(inputFile.RecTable.dataRaw{exp,1}{1,n}))/SR(exp,1),allData{exp}{1,n})  % Stimulus protocol (pA/s)
                end
            end 
        else
            if imgshow
                figure('Name',a);
                plot((1:length(inputFile.RecTable.dataRaw{exp,1}{1,n}))/SR(exp,1),allData{exp}{1,n})
            end
        end
        
        if n == 1
            list{end+1} = a;
        end
      %  T = table (allData{exp,1}{1,n}); 
      %  if ~contains(a,'*') == 1
      %     writetable(T, [a,'.txt']);  %% Save file as txt
      %  else
      %      writetable(T, [erase(a,"*"),'.txt']);  %% Save file as txt
      %  end
    end          
end

%% define parameters 2
baseline = 0.2;    % YC: 0.2 s
stim_dur = 1;    % YC: 1 s
threshold = 50; % AP threshold: 50 V/s; https://www.nature.com/articles/nn.2203   ; https://www.jneurosci.org/content/26/6/1854
hump_time = 0.3; % 0.3 s; https://www.jneurosci.org/content/jneuro/39/1/125.full.pdf
AHP_dur = 0.2; % time window search for AHP minimum, second
AP_dur = 0.005; % AP duration, second

A = exist('sweeps','var');
if A
    clear sweeps stimulus
end

[indx,tf] = listdlg('ListString',list,'ListSize',[200,600],'PromptString','Select data for intrinsic properties analysis');
% sample_point = size(allData{indx(1)}{1,1},1);  
list_selected = list(indx)';
for i = 1: size(list_selected,1)
    list_selected_cellnum(i) = inputFile.RecTable.Experiment(indx(i)) ;
end
list_selected_cellnum_uniq = unique(list_selected_cellnum)';


%% Analysis
for b = 1:size(indx,2)
    sweeps{1,b} = allData{indx(b)}{1,1}; %% voltage trace, mV
    sweeps{2,b} = allData{indx(b)}{1,2}; %% stimulus (current injection), pA
    sweeps{3,b} = SR(indx(b))*diff(sweeps{1,b})/1000;  %% differential of voltage trace, V/s,
    count = 1;
    for bb = 1:size(sweeps{2,b},2)
        %% input resistance
        base = mean(sweeps{2,b}(1:baseline*SR(indx(b)),bb));
        stim_mean = mean(sweeps{2,b}(baseline*SR(indx(b)):stim_dur*SR(indx(b)),bb));
        stimulus(bb,b) = (stim_mean - base);    %% nonfilled blanks will be zeros
        
        base = median(sweeps{1,b}(1:baseline*SR(indx(b)),bb));
        stim_mean = median(sweeps{1,b}(baseline*SR(indx(b)):stim_dur*SR(indx(b)),bb));
        response(bb,b) = (stim_mean - base);    %% nonfilled blanks will be zeros 
        %% First spike latency
        MinPeakDistance = 2; % 2 ms
        if max(sweeps{1,b}(baseline*SR(indx(b)):stim_dur*SR(indx(b)),bb)) > 0
            [pks,locs] = findpeaks(sweeps{1,b}(baseline*SR(indx(b)):stim_dur*SR(indx(b)),bb),'MinPeakHeight',0,'MinPeakDistance',MinPeakDistance*SR(indx(b))/1000);
        else
            pks = nan;
            locs = nan;
        end               
        
        %% Results gathering
        if ~isnan(pks)
            SpikeNumber(bb,b) = size(pks,1);
            FSL(bb,b) = (1000*locs(1)/SR(indx(b)));    % millisecond           
            FS_number(bb,b) = nan;
            while count == 1
                FS_number(bb,b) = SpikeNumber(bb,b);
                count = count -1;
            end
                        
            threshold_idx = find(sweeps{3,b}(baseline*SR(indx(b)):stim_dur*SR(indx(b)),bb) > threshold, 1); %% AP threshold, mV
            AP_threshold(bb,b) = sweeps{1,b}(baseline*SR(indx(b))+threshold_idx,bb); %% AP threshold, mV
            AP_amp(bb,b) = pks(1) - AP_threshold(bb,b); %% AP amplitude, mV
            Hump_ES(bb,b) = AP_threshold(bb,b)-sweeps{1,b}(hump_time*SR(indx(b)),bb); %% Depo Hump amp, mV
            Hump_LS(bb,b) = max(sweeps{1,b}(baseline*SR(indx(b)):hump_time*SR(indx(b)),bb))-sweeps{1,b}(hump_time*SR(indx(b)),bb); %% Depo Hump amp, mV
            
            if SpikeNumber(bb,b) > 1
                threshold_sec_idx = find(sweeps{3,b}((baseline+AP_dur)*SR(indx(b))+threshold_idx:end,bb) > threshold, 1); %% second AP threshold, mV
                [AHP_minimum(bb,b),AHP_min_idk(bb,b)] = max(-sweeps{1,b}(baseline*SR(indx(b))+threshold_idx:(baseline+AP_dur)*SR(indx(b))+threshold_idx+threshold_sec_idx,bb)); %% AHP peak, mV
            elseif SpikeNumber(bb,b) == 1
                [AHP_minimum(bb,b),AHP_min_idk(bb,b)] = max(-sweeps{1,b}(baseline*SR(indx(b))+threshold_idx:baseline*SR(indx(b))+threshold_idx+AHP_dur*SR(indx(b)),bb)); %% AHP peak, mV
            end
            AHP_minimum(bb,b) = -AHP_minimum(bb,b); %% AHP peak, mV
            AHP_latency(bb,b) = 1000*AHP_min_idk(bb,b)/SR(indx(b)); %% AHP latency, millisecond
            AHP_amp(bb,b) = AP_threshold(bb,b)-AHP_minimum(bb,b); %% AHP amplitude, mV
            AHP_half_amp(bb,b) = AHP_amp(bb,b)/2; %% AHP half amplitude, mV
            
            ii = baseline*SR(indx(b))+threshold_idx+AHP_min_idk(bb,b)+2;
            kk = baseline*SR(indx(b))+threshold_idx+AHP_min_idk(bb,b)+2;
            while sweeps{1,b}(ii,bb)-AHP_minimum(bb,b)< AHP_half_amp(bb,b)
                ii = ii-1;
            end
            
            while sweeps{1,b}(kk,bb)-AHP_minimum(bb,b)< AHP_half_amp(bb,b)
                kk = kk+1;
                if kk > size(sweeps{1,b},1)
                    break
                end
            end
            if kk > size(sweeps{1,b},1)
                AHP_half_width(bb,b) = nan;   %  AHP_half_width not available
            else
                AHP_half_width(bb,b) = 1000*(kk-ii)/SR(indx(b)); %% AHP half width, ms
            end
            
        else
            SpikeNumber(bb,b) = nan;  % distinguish from auto filled zeros
            FSL(bb,b) = nan;
            FS_number(bb,b) = nan;
            AP_amp(bb,b) = nan;
            AP_threshold(bb,b) = nan;
            Hump_ES(bb,b) = nan;
            Hump_LS(bb,b) = nan;
            AHP_minimum(bb,b) = nan;
            AHP_latency(bb,b) = nan;
            AHP_amp(bb,b) = nan;
            AHP_half_amp(bb,b) = nan;
            AHP_half_width(bb,b) = nan;
        end
    end
    stimulus(stimulus == 0) = nan;   
    response(response == 0) = nan;
    
    i = size(sweeps{2,b},2);
    k = 1;
    while stimulus(i,b) > 31    %  30 to -30 pA
        i = i-1;
        if i == 0
            break
        end
    end
    
    while stimulus(k,b) < -31   %  30 to -30 pA
        k = k+1;
    end 
    
    if i ~= 0 && k ~= 0 && i-k > 3
        x = stimulus(k:i,b);
        y = response(k:i,b);
        mdl = fitlm(x,y);
        slope = (mdl.Fitted(end)-mdl.Fitted(1))/(stimulus(i,b)-stimulus(k,b));    % input resistance, giga ohm
        LinearRegression(1,b) = slope;
        LinearRegression(2,b) = mdl.Rsquared.Ordinary;
    else
        LinearRegression(1,b) = nan;
        LinearRegression(2,b) = nan;
    end
end

% figure % input resistance
% plot(stimulus,response,'x')
%% Sag ratio and membrane time constant
sag_window = 0.3;  % sag detection from 0 to 0.3 s after stim
smooth_point = 5;  % total number of sample points used for mean to determine the peak 
sag_ss_window = 0.3; % time window to measure steady state voltage at last 0.3 s of current injection
fit_window = 0.1;   % time window for exponential fitting, second
sag_ratio = nan(1,size(indx,2));


for b = 1:size(indx,2)
    for bb = 1:size(sweeps{1,b},2)
        if round(stimulus(bb,b)) == -50 && -200  % -50 or -200 pA current inection
            base = median(sweeps{1,b}(1:baseline*SR(indx(b)),bb));
            [sag_amp,sag_amp_idk] = max(-sweeps{1,b}(baseline*SR(indx(b)):(baseline+sag_window)*SR(indx(b)),bb));
            sag_steadystate = median(sweeps{1,b}((baseline+stim_dur-sag_ss_window)*SR(indx(b)):(baseline+stim_dur)*SR(indx(b)),bb));
            sag_steadystate_all(bb,b) = sag_steadystate-base;
            
            smooth_bin = (smooth_point-1)/2;
            sag_amp_smooth = 0;          
            while smooth_bin > 0
                sag_amp_smooth = sag_amp_smooth + sweeps{1,b}(baseline*SR(indx(b))+sag_amp_idk+smooth_bin,bb) + sweeps{1,b}(baseline*SR(indx(b))+sag_amp_idk-smooth_bin,bb);              
                smooth_bin = smooth_bin-1;
            end
            sag_amp_smooth = sag_amp_smooth + sweeps{1,b}(baseline*SR(indx(b))+sag_amp_idk);
            sag_amp_smooth = sag_amp_smooth/smooth_point;
            sag_amp_smooth_all(bb,b) = (sag_amp_smooth - base);    %% nonfilled blanks will be zeros
            sag_ratio(1,b) = sag_steadystate_all(bb,b)/sag_amp_smooth_all(bb,b);  %% sag ratio, median voltage from steady state devided by smoothed peak voltage of sag period
            
            %% rebound spike check
            [pks_rebound,locs_rebound] = findpeaks(sweeps{1,b}((baseline+stim_dur)*SR(indx(b)):(baseline+stim_dur+fit_window)*SR(indx(b)),bb),'MinPeakHeight',0,'MinPeakDistance',MinPeakDistance*SR(indx(b))/1000);
            
            %% exponential fitting for time constant measurement           
            y = (sweeps{1,b}((baseline+stim_dur)*SR(indx(b)):(baseline+stim_dur+fit_window)*SR(indx(b)),bb))';
            x = 1000*(1:1:size(y,2))/SR(indx(b));  % from 1.2 to 1.3 s, unit: ms
            
            % Define the exponential fit model
            expFitModel = fittype('AA * exp(-x/tau) + c', 'independent', 'x', 'dependent', 'y');

            % Set initial parameter values for the fit
            initialParams = [-10, 1, 1];

            % Fit the data to the model
            if isempty(pks_rebound)
                [fitResult,gof] = fit(x', y', expFitModel, 'StartPoint', initialParams);

            % Display the coefficients
%             disp('Exponential Fit Coefficients:');
%             disp(['AA: ', num2str(fitResult.AA)]);
%             disp(['tau: ', num2str(fitResult.tau)]); % ms
%             disp(['c: ', num2str(fitResult.c)]);
%             disp(['R-square: ', num2str(gof.rsquare)]);
            
            % Plot the original data and the fitted curve
                Plot = 0;
                if Plot           
                    list_name = (list(indx(b)));
                    figure('Name',list_name{:});
                    plot(x, y, 'DisplayName', char({'Original Data'}));
                    legend('Location', 'Best');
                    hold on;
                    plot(fitResult, 'r-')  %, 'DisplayName', 'Exponential Fit');
                    title('Exponential Fitting');
                    xlabel('time (ms)');
                    ylabel('Vm (mV)');          
                    text(60, -70, sprintf('f(x) = %.1f\\cdote^{-x/%.3f}%+.1f \n R^2 = %.4f', fitResult.AA,fitResult.tau,fitResult.c,gof.rsquare));
                    legend('Location', 'Best');
                    grid on;
                    hold off;
                end
                Expo_fit(1,b) = fitResult.tau;  % fitted tau, ms
                Expo_fit(2,b) = gof.rsquare;    % Rsquare
            else
                Expo_fit(1,b) = nan;    % fitted tau, ms
                Expo_fit(2,b) = nan;    % Rsquare
            end
        elseif b == size(indx,2) && bb == size(sweeps{1,b},2) && round(stimulus(bb,b)) ~= (-50 && -200)
            Expo_fit(1,b) = nan;    % fitted tau, ms
            Expo_fit(2,b) = nan;    % Rsquare
        else              
            sag_amp_smooth_all(bb,b) = nan;
            sag_steadystate_all(bb,b) = nan;           
        end
    end
end

%% Finalize results

sz = [size(list_selected_cellnum_uniq,1) 17];
varTypes = {'double','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double'};
varNames = {'Date' 'Cell#' 'Rs (MOhm)' 'Rs_comp (%)' 'Input resistance (GOhm)' 'Rin Rsquare' 'Membrane time constant (ms)' 'fitted tau Rsquare' 'Sag ratio' 'First spike number' 'First spike latency (ms)' 'Depo hump amp ES (mV)' 'Depo hump amp LS (mV)' 'AP amp (mV)' 'AP threshold (mV)' 'First AHP latency (ms)' 'AHP half-width (ms)'};
T = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);

for i = 1:size(list_selected_cellnum_uniq,1)
    
    T{i,2} = list_selected_cellnum_uniq(i);  % Cell#
    
    cell_idx = find(list_selected_cellnum == list_selected_cellnum_uniq(i));
    T{i,3} = inputFile.RecTable.Rs_uncomp{indx(cell_idx(1))}{1,1}(1)/1000000; % Rs-uncompensated value, MOhm
    T{i,4} = inputFile.RecTable.RsFractionComp{indx(cell_idx(1))}{1,1}(1)*100; % Rs-compensation fraction, in percentage
    
    Rsquare_max_idx = find(LinearRegression(2,:) == max(LinearRegression(2,cell_idx)));
    if ~isempty(Rsquare_max_idx)
        T{i,5} = LinearRegression(1,Rsquare_max_idx); % Input resistance (GOhm)
        T{i,6} = LinearRegression(2,Rsquare_max_idx); % Linear regression Rsquare for Rin
    else
        T{i,5} = nan;
        T{i,6} = nan;
    end
    
    
    if max(Expo_fit(2,cell_idx)) ~= 0
        Rsquare_max_idx = find(Expo_fit(2,:) == max(Expo_fit(2,cell_idx)));
    else
        Rsquare_max_idx = [];
    end
    if ~isempty(Rsquare_max_idx)
        T{i,7} = Expo_fit(1,Rsquare_max_idx); % Membrane time constant (ms)
        T{i,8} = Expo_fit(2,Rsquare_max_idx); % fitted tau Rsquare
    else
        T{i,7} = nan;
        T{i,8} = nan;
    end
    
    temp = sag_ratio(cell_idx);
    temp_idx = find(temp > 0);
    T{i,9} = mean(temp(temp_idx),'omit'); % mean Sag ratio
    
    cat_temp = reshape(FS_number(:,cell_idx),[size(cell_idx,2)*size(FS_number(:,cell_idx),1),1]);   
    cat_temp_idx = find(cat_temp > 0);
    T{i,10} = mean(cat_temp(cat_temp_idx));  % mean First Spike number       
 
    cat_temp = reshape(FSL(:,cell_idx),[size(cell_idx,2)*size(FSL(:,cell_idx),1),1]);
    T{i,11} = mean(cat_temp(cat_temp_idx),'omitnan'); % mean First Spike Latency, ms
    
    cat_temp = reshape(Hump_ES(:,cell_idx),[size(cell_idx,2)*size(Hump_ES(:,cell_idx),1),1]);
    T{i,12} = mean(cat_temp(cat_temp_idx),'omit'); % mean Depo hump amp ES, mV
    
    cat_temp = reshape(Hump_LS(:,cell_idx),[size(cell_idx,2)*size(Hump_LS(:,cell_idx),1),1]);
    T{i,13} = mean(cat_temp(cat_temp_idx),'omit'); % mean Depo hump amp LS, mV
    
    cat_temp = reshape(AP_amp(:,cell_idx),[size(cell_idx,2)*size(AP_amp(:,cell_idx),1),1]);
    T{i,14} = mean(cat_temp(cat_temp_idx),'omit'); % mean first AP amp, mV
    
    cat_temp = reshape(AP_threshold(:,cell_idx),[size(cell_idx,2)*size(AP_threshold(:,cell_idx),1),1]);
    T{i,15} = mean(cat_temp(cat_temp_idx),'omit'); % mean first AP threshold, mV
    
    cat_temp = reshape(AHP_latency(:,cell_idx),[size(cell_idx,2)*size(AHP_latency(:,cell_idx),1),1]);
    T{i,16} = mean(cat_temp(cat_temp_idx),'omit'); % mean first AHP latency, ms
    
    cat_temp = reshape(AHP_half_width(:,cell_idx),[size(cell_idx,2)*size(AHP_half_width(:,cell_idx),1),1]);
    T{i,17} = mean(cat_temp(cat_temp_idx),'omit'); % mean first AHP half-width, ms
end
%% Output file
cd(filepath);
str = 'J:\Data\mPFC-ACC\mPFC-ACC\';
T{:,1} = str2double(erase(filepath,str));
writetable(T,strcat(erase(filepath,str),'.xlsx'))


%%
disp('job done');
sound(sin(1:3000));


% % figure
% % plot(sweeps{3,b}(baseline*SR(indx(1)):stim_dur*SR(indx(b)),bb),'x')
% % figure
% % plot(sweeps{1,b}(baseline*SR(indx(b)):stim_dur*SR(indx(b)),bb),'x')
% % 
% % figure
% % plot(sweeps{1,b}(:,bb),'x')
% % hold on
% % plot(ii,sweeps{1,b}(ii,bb),'rx')
% % plot(kk,sweeps{1,b}(kk,bb),'r+')
% 
% 
% function sse = sseval(x,tdata,ydata)
%             A = x(1);
%             lambda = x(2);
%             sse = sum((ydata - A*exp(-lambda*tdata)).^2);
% end
% 
% x=[10 12.5 15 17.5 20 22.5 25 27.5 30 32.5 35 37.5 40 42.5 45 47.5 50];
% y=[62.1 77.3 92.5 104 112.9 121.9 125 129.4 134 138.2 142.3 143.2 144.6 147.2 147.8 149.1 150.9];
% F = @(c,x) c(1).*exp(c(2).*x)+c(3);                                     % Objective Function
% B = fminsearch(@(c) norm(y - F(c,x)), [-200; -1; 100])                  % Estimate Parameters
% figure
% plot(x, y, 'pg')
% hold on
% plot(x, F(B,x), '-r')
% hold off
% grid
% xlabel('x')
% ylabel('f(x)')
% text(27, 105, sprintf('f(x) = %.1f\\cdote^{%.3f\\cdotx}%+.1f', B))
% 
% 
% %%
% 
% x=[10 12.5 15 17.5 20 22.5 25 27.5 30 32.5 35 37.5 40 42.5 45 47.5 50];
% y=[62.1 77.3 92.5 104 112.9 121.9 125 129.4 134 138.2 142.3 143.2 144.6 147.2 147.8 149.1 150.9];
% F = @(c,x) c(1)*exp(-x/c(2))+c(3);     
% B = fminsearch(@(c) norm(y - F(c,x)), [-10; 1; 1]);                  % Estimate Parameters
% 
% x = x(1:100);
% y = y(1:100);
% 
% x = 1:1:size(y,2);
% 
% figure
%             plot(x, y, 'pg')
%             hold on
%             plot(x, F(B,x), '-r')
%             hold off
%             grid
%             xlabel('x')
%             ylabel('f(x)')
%             text(27, 105, sprintf('f(x) = %.1f\\cdote^{%.3f\\cdotx}%+.1f', B))
%             
%             figure
%             plot(F(B,x))
%             
%             A = B(1);
%             tau = B(2);
%             C = B(3);
%             
%             for i = 1:size(x,2)            
%                 exp_f(i) = A*exp(-i/tau)+C;
%             end
%             
%             figure
%             plot(exp_f)
%             
% %%
% x = [1, 2, 3, 4, 5];
% y = [10, 30, 80, 200, 500];
% 
% x = 1:1:size(y,2);
% x = x';
% y = y';
% 
% 
% initialParams = [-1, 1, 1];
% g = fittype('a*exp(-x/b)+c');
%             f = fit(x,y,g,'StartPoint', initialParams,'MaxIter', 6000,'MaxFunEvals',24000,'Display','final');
%             
%             figure
%             plot(f,x,y,'g')
% %%
% 

