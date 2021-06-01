clc
clear
close all
set(0,'DefaultAxesFontSize',15)

%%
addpath(genpath('rt_bench_may_2021'))

%% USE h5disp to view structure of h5 file 
h5disp('gen_decode_n2_len30.h5')

%{
'001 HmmGenerator  DECODING DISABLED 12  Spike'
'002 HmmGenerator  DECODING DISABLED 12  TrueState'
'003 HmmDecoder BENCHMARK VERSION 11  spikes inout'
'004 HmmDecoder BENCHMARK VERSION 11  state out' 
'005 Performance Measurement 16  Comp Time ns' 
'006 Performance Measurement 16  Realtime Period ns'
%}


channelKey = {'spikes',...
    'trueState',...
    'spike_thru',...
    'decodedState',...
    'comp',...
    'rt per'};

getIdx = @(keyName) find(contains(channelKey,keyName));
spikesIdx = getIdx('spikes')
trueStateIdx = getIdx('trueState')
spikesThruIdx = getIdx('spike_thru')
decodedStateIdx = getIdx('decodedState')

compIdx = getIdx('comp')
rtPeriodIdx = getIdx('rt per')


%%

%h5read('gen_decode_n2_len30.h5', "/Trial1/Synchronous Data/Channel Data")

nstates = [2,3,4,5,6,7];
bufflens = [30,100,300,1000];

genFilename = @(src,num1,num2) sprintf('%s_decode_n%i_len%i.h5',src, num1,num2);
readFun = @(src, num1, num2) h5read(genFilename(src,num1,num2), "/Trial1/Synchronous Data/Channel Data")

% readFun('gen',2,30)


%%  Read in all the data. note, depending on the size of the sweep, this may have to be done in slices

NSweep_gen = {};
NSweep_tdt = {};
NSweep_tdt_keys = {};
NSweep_tdt_ct = [];
NSweep_gen_ct = [];





for n = nstates
    NSweep_gen{end+1} = readFun('gen',n,300);
    NSweep_tdt{end+1} = readFun('tdt',n,300);
    
%     tdtInfo_ = h5info(genFilename('tdt',n,300));
%     channelNames = {tdtInfo_.Groups.Groups(3).Datasets.Name};
%     NSweep_tdt_keys{end+1} = channelNames;
%     NSweep_tdt_ct(end+1) = find(contains(channelNames,'Comp Time'))

    NSweep_gen_ct(end+1) = getIndexFromParams('gen',n,300,'Comp Time')
    NSweep_tdt_ct(end+1) = getIndexFromParams('tdt',n,300,'Comp Time')
    
end

LenSweep_gen = {};
LenSweep_tdt = {};
LenSweep_tdt_keys = {};

LenSweep_gen_ct = [];
LenSweep_tdt_ct = [];

for b = bufflens
    LenSweep_gen{end+1} = readFun('gen',2,b);
    LenSweep_tdt{end+1} = readFun('tdt',2,b);
    
    LenSweep_gen_ct(end+1) = getIndexFromParams('gen',2,b,'Comp Time')
    LenSweep_tdt_ct(end+1) =  getIndexFromParams('tdt',2,b,'Comp Time')
    
    
end


% return
%%

base_time = 1e-9;
%others: us = 1e-6; ms = 1e-3;
time_unit = 'ms';
time_conv = base_time/1e-3;

%%
gray = [1,1,1]*.7;
lightBlue = gray;
lightBlue(3)=1;
plotSubset = [1:1e3] + 1e3;

YL = [0,1];


posLeft = [50, 50, 400, 255];
posRight = [440, 50, 650, 255];

bigMS = 9;
smallMS = 4;

figure(1)
clf
hold on

i=1;
genData = NSweep_gen{i};
tdtData = NSweep_tdt{i};
gy =  genData(NSweep_gen_ct(i),:)*time_conv;
%     ctIdx = NSweep_tdt_ct(i);
ty = tdtData(NSweep_tdt_ct(i),:)*time_conv;

%     gQ = quantile(gy,[0.0, 0.1, 0.5, 0.9, 1.0]);
gQ = quantile(gy,0);
tQ = quantile(ty,0);
xVal = nstates(i);
xVal2 = xVal+.15;
plot(xVal,mean(gy),'k+','MarkerSize',bigMS,'LineWidth',2)
plot(xVal, gQ,'k+','MarkerSize',smallMS)

plot(xVal2, mean(ty), 'bo','MarkerSize',bigMS,'LineWidth',2)
plot(xVal2, tQ, 'bx','MarkerSize',smallMS)
    
for i = 1:length(nstates)
    
    genData = NSweep_gen{i};
    tdtData = NSweep_tdt{i};
    gy =  genData(NSweep_gen_ct(i),:)*time_conv;
%     ctIdx = NSweep_tdt_ct(i);
    ty = tdtData(NSweep_tdt_ct(i),:)*time_conv;
    
%     gQ = quantile(gy,[0.0, 0.1, 0.5, 0.9, 1.0]);
    gQ = quantile(gy,[0, 1.0]);
    tQ = quantile(ty,[0, 1.0]);


    
    xVal = nstates(i);
    xVal2 = xVal+.15;
    
    plot(xVal, gy(plotSubset),'.','Color',gray)
    plot(xVal, gQ,'k+','MarkerSize',smallMS)
    plot(xVal,mean(gy),'k+','MarkerSize',bigMS,'LineWidth',2)
    
    
    plot(xVal2, ty(plotSubset),'.','Color',lightBlue)
     plot(xVal2, mean(ty), 'bo','MarkerSize',bigMS,'LineWidth',2)
     plot(xVal2, tQ, 'bx','MarkerSize',smallMS)

     
end

legend('mean, from gen.','min / max, from gen.','mean, from TDT','min / max, from tdt.','location','Northwest')


% ylim(YL)
xlim([1,max(nstates)+1])
xlabel({'# HMM states',''})
ylabel('compute time [ms]')
set(gcf,'Position',posLeft+[0,300,0,0])


figure(2)
clf
hold on
plot(gy(plotSubset),'k','LineWidth',2)
plot(ty(plotSubset),'b','LineWidth',2)
ylim(YL)
ylabel('compute time [ms]')
xlabel({'time','[samples]'})
set(gcf,'Position',posRight+[0,300,0,0])
legend('from internally generated spikes', 'from TDT spikes','location','southwest')


%%


figure(3)
clf
hold on
for i = 1:length(bufflens)
    genData = LenSweep_gen{i};
    tdtData = LenSweep_tdt{i};
    
%     ctIdx = LenSweep_tdt_ct(i);

    
    
    gy =  genData(LenSweep_gen_ct(i),:)*time_conv;
    ty = tdtData(LenSweep_tdt_ct(i),:)*time_conv;
    
%       ty = tdtData(ctIdx,:)*time_conv;
    
%     gQ = quantile(gy,[0.0, 0.1, 0.5, 0.9, 1.0]);
    gQ = quantile(gy,[0, 1.0]);
    tQ = quantile(ty,[0, 1.0]);

    xVal = bufflens(i);
    xVal2 = xVal*1.15;
    
    plot(xVal, gy(plotSubset),'.','Color',gray)
    plot(xVal, gQ,'k+','MarkerSize',smallMS)
    plot(xVal,mean(gy),'k+','MarkerSize',bigMS,'LineWidth',2)
    
    
    plot(xVal2, ty(plotSubset),'.','Color',lightBlue)
     plot(xVal2, mean(ty), 'bo','MarkerSize',bigMS,'LineWidth',2)
     plot(xVal2, tQ, 'bx','MarkerSize',smallMS)
    
end
xlim([bufflens(1)/3, bufflens(end)*3])
ylim(YL)
ylabel('compute time [ms]')
set(gca,'XScale','log')
xlabel({'buffer length','[# samples]'})
set(gcf,'Position',posLeft)

figure(4)
clf
hold on
plot(gy(plotSubset),'k','LineWidth',2)
plot(ty(plotSubset),'b','LineWidth',2)
ylim(YL)
ylabel('compute time [ms]')
xlabel({'time','[samples]'})
set(gcf,'Position',posRight)


function [foundIndex] = getIndexFromParams(src,num1,num2, channelName)
    genFilename = @(src,num1,num2) sprintf('%s_decode_n%i_len%i.h5',src, num1,num2);

    hi = h5info(genFilename(src,num1,num2));
    channelNames = {hi.Groups.Groups(3).Datasets.Name};
% getChannelNamesFromParams = @(src,num1, num2, channelName) 
    foundIndex = find(contains(channelNames, channelName));
    
end







