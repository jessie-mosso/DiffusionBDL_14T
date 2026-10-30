clear;clc;close all; 

%% add FID-A master to Matlab path
LBall=12;
doprocessing=true;
dosavesum=true;
dosavefid=false;
docopyrefscan=true;

lcmodelfac=1;

%% path to raw data
rawpath=[pwd '/raw']; 
processedpath=[pwd '/processed']; 

filelist = dir([rawpath '/*ser.mat']); 

%% step 1 - do preprocessing with FIDA
if doprocessing
for file=1:length(filelist)
    disp(filelist(file).name)
    f1=figure;
    
    %% 1-Convert Bruker study to FID A structure
    out0=convertNicoStudyFidAformat(rawpath, filelist(file).name);
    %apply LB
    dw=out0.dwelltime;
    tt=[0:dw:dw*(out0.n-1)];
    out0lb=out0;
    fids0lb=out0lb.fids.*repmat(exp(-tt*pi*LBall).',1,out0.averages);
    out0lb.fids=fids0lb;
    out0lb.specs=fftshift(fft(out0lb.fids.',[],2),2).';
    figure(f1);
    subplot(2,2,1)
    for k=1:out0.averages
        plot(real(out0.specs(:,k)))
        hold on 
        xlim([2900,3050])
    end 
    title(['0-ini'])
    figure(f1);
    subplot(2,2,2)
    for k=1:out0lb.averages
        plot(real(out0lb.specs(:,k)))
        hold on 
        xlim([2900,3050])
    end 
    title(['1-ini - LB' num2str(LBall)])

    %% 2-align av.
%     [out1lb,fs,phs]=op_alignAverages_jm(out0lb,0.2,'y'); %requires prior LB
    [out1lb,fs,phs]=op_alignAverages_fd_jm(out0lb,7,8,0.5,'y');%5.2,9,0.5,'y');
    
    figure(f1);
    subplot(2,2,3)
    for k=1:out1lb.averages
        plot(real(out1lb.specs(:,k)))
        hold on 
        xlim([2900,3050])
    end
    title(['2-spec. reg. JNear - LB' num2str(LBall)])
    
    %remove LB
    out1=out1lb;
    fids1=out1lb.fids.*repmat(exp(tt*pi*LBall).',1,out1lb.averages);
    out1.fids=fids1;
    out1.specs=fftshift(fft(out1.fids.',[],2),2).';
       

    %% 3-outlier removal 
    [out2,metric,badAverages]=op_rmbadaverages_jm(out1,1.5,'f'); %performs 10Hz LB inside
%     out2=out1;
    
    %apply LB
    out2lb=out2;
    fids2lb=out2.fids.*repmat(exp(-tt*pi*LBall).',1,out2.averages);
    out2lb.fids=fids2lb;
    out2lb.specs=fftshift(fft(out2lb.fids.',[],2),2).';
    
    figure(f1);
    subplot(2,2,4)
    for k=1:out2lb.averages
        plot(real(out2lb.specs(:,k)))
        hold on 
        xlim([2900,3050])
    end
    title(['3-rmv motion corrupted av. - LB' num2str(LBall)])
    suptitle(filelist(file).name)

    %% 4-add all the info to the Matlab study structure
    load([rawpath '/' filelist(file).name]);
    study.fidaprocess.phs=phs;
    study.fidaprocess.fs=fs;
    study.fidaprocess.metric=metric;
    study.fidaprocess.badAverages=badAverages; 
    study.params.nt=size(out2.fids,2); 
    study.multiplicity=size(out2.fids,2);
    study.process.apodparam1=zeros(1,size(out2.fids,2));
    study.process.apodparam2=zeros(1,size(out2.fids,2));
    study.process.phasecorr0=zeros(1,size(out2.fids,2));
    study.process.phasecorr1=zeros(1,size(out2.fids,2));
    study.data.real=zeros(size(out2.fids,2),1,size(out2.fids,1));
    study.data.real(:,1,:)=real(out2.fids.');
    study.data.imag=zeros(size(out2.fids,2),1,size(out2.fids,1));
    study.data.imag(:,1,:)=imag(out2.fids.');
    if ~exist([processedpath '/'], 'dir')
       mkdir([processedpath '/'])
    end
    save([processedpath '/' filelist(file).name(1:end-4) '_processed.mat'],'study')   
end 
end 
    
%% STEP 2 - APPLY PREPROCESSING, SUM AND PHASE THE SUM 

if dosavesum
for file=1:length(filelist)
    load([processedpath '/' filelist(file).name(1:end-4) '_processed.mat'],'study')   
    
    fidmocor=squeeze(study.data.real)+1i*squeeze(study.data.imag);
    sumfid=sum(fidmocor);
    sumfid=sumfid.*lcmodelfac;

    ftcorr=fftshift(fft(sumfid,[],2),2); 
    figure; 
    plot(real(ftcorr))
    hold on 
    load([rawpath '/' filelist(file).name(1:end-4) '.mat'],'study')   
    fidini=squeeze(study.data.real)+1i*squeeze(study.data.imag);
    sumfidini=sum(fidini)./size(fidini,1); 
    sumfidini=[sumfidini(77:end),zeros(1,76)]; 
    ftini=fftshift(fft(sumfidini,[],2),2); 
    plot(real(ftini))
    legend('corr','ini')
    
    
    %% save
    study.data.real=zeros(1,1,study.np/2);
    study.data.imag=zeros(1,1,study.np/2);

    study.data.real(1,1,:)=real(sumfid); 
    study.data.imag(1,1,:)=imag(sumfid);

    study.multiplicity=1;
    study.process.lsfid=0;
    study.process.apodparam1=0;
    study.process.apodparam2=0;
    study.process.phasecorr0=0;
    study.process.phasecorr1=0;
    study.process.B0=zeros(1,study.np/2);
        
        
        
    filename=filelist(file).name(1:end-4);
    if ~exist([processedpath '/sum/'], 'dir')
        mkdir([processedpath '/sum/'])
    end
    save([processedpath '/sum/SUM_' filename '_processed.mat'],'study');
end
end 


%% save also fid files already summed on the scanner
if dosavefid
filelist = dir([rawpath '/*fid.mat']);
for file=1:length(filelist)
    load([rawpath '/' filelist(file).name(1:end-4) '.mat'],'study')   
    
    fidmocor=squeeze(study.data.real).'+1i*squeeze(study.data.imag).';
    grpdly=round(study.params.grpdly)+1;
    fidmocor=[fidmocor(grpdly:end),zeros(1,grpdly-1)];

    sumfid=fidmocor;%./study.params.nt;

    ftcorr=fftshift(fft(sumfid,[],2),2); 
    figure; 
    plot(real(ftcorr))
          
    
    %% save
    study.data.real=zeros(1,1,study.np/2);
    study.data.imag=zeros(1,1,study.np/2);

    study.data.real(1,1,:)=real(sumfid); 
    study.data.imag(1,1,:)=imag(sumfid);

    study.multiplicity=1;
    study.process.lsfid=0;
    study.process.apodparam1=0;
    study.process.apodparam2=0;
    study.process.phasecorr0=0;
    study.process.phasecorr1=0;
    study.process.B0=zeros(1,study.np/2);
        
        
        
    filename=filelist(file).name(1:end-4);
    if ~exist([processedpath '/sum/'], 'dir')
        mkdir([processedpath '/sum/'])
    end
    save([processedpath '/sum/SUM_' filename '.mat'],'study');
end
end 
%%
if docopyrefscan
    filelist = dir([rawpath '/*fidrefscan.mat']);
    for file=1:length(filelist)
    copyfile([rawpath '/' filelist(file).name(1:end-4) '.mat'],[processedpath '/sum/']) 
    end 
end 
    
    
   
