clear;clc;close all; 

goodindshamw4=[1:3];
goodindbdlw4=[4:13];


goodindshamw6=[1:3,6,9:11];
goodindbdlw6=[4:5,7:8,12:13];


week=4;

load(['DWSpupstdd60_week' num2str(week) '_qc.mat'])


bval=[50,100,500,1000,3000,5000,10000,15000];
 
bval2=bval(2:end)./1000; 
listtoplot=DWSpupstdd60.listttoplot;


x0_sticks=0.6; %inital guess for diff

%% Fit Dintra in Indiv Decays
Dsticks_week4=zeros(length(listtoplot),size(DWSpupstdd60.conctotmetrat,2),2);

for k=1:length(DWSpupstdd60.listttoplot)
%     figure; 
    for s=1:length(goodindshamw4)
        concsham0=squeeze(DWSpupstdd60.conctotmetrat(k,goodindshamw4(s),:));
        crlbsham0=squeeze(DWSpupstdd60.conctotmetrat(k,goodindshamw4(s),:).*DWSpupstdd60.crbtotmetrat(k,goodindshamw4(s),:));
        concsham=concsham0./concsham0(1);
        aa=concsham0;
        bb=concsham0(1); 
        saa=crlbsham0;
        sbb=crlbsham0(1); 
        crlbsham=aa./bb.*sqrt((saa./aa).^2+(sbb/bb).^2);
%         crlbsham=concsham.*sqrt((crlbsham0./concsham0).^2+(crlbsham0(1)/concsham0(1)).^2);

        kkk = ~isnan(concsham);
        weighted_s_sticks = @(diff,b) 1./(crlbsham(kkk).^2).*sqrt(pi./(4.*b.*diff)).*erf(sqrt(b.*diff));
        [pfit, resnorm, residual, exitflag, output, lambda, J] =lsqcurvefit(@(diff,b) weighted_s_sticks(diff,bval2(kkk).'),x0_sticks,bval2(kkk),1./(crlbsham(kkk).^2).*concsham(kkk),0,1);

        Dsticks_week4(k,goodindshamw4(s),:)=pfit;

        J = full(J);
        n = numel(concsham(kkk));        % number of data points
        p = numel(pfit);         % number of parameters
        dof = n - p;             % degrees of freedom
        
        resid    = concsham(kkk).'-sqrt(pi./(4.*bval2(kkk).*pfit)).*erf(sqrt(bval2(kkk).*pfit));    % residuals
        sigma2   = crlbsham(kkk).^2;                       % measurement variance per b
        nu       = numel(bval2(kkk)) - 1;                % degrees of freedom
        
        chi2     = sum( resid.'.^2 ./ sigma2 );        % weighted chi-square
        chi2_red = chi2 / nu                        % reduced chi-square

        % Estimate of error variance (mean squared residual)
        sigma2 = resnorm / dof;  % resnorm = sum(residual.^2)
        
        % Parameter covariance matrix
        Cov = sigma2 * inv(J' * J);
        Dsticks_week4_cov(k,goodindshamw4(s),:)=Cov(1,1);

        
    end 
    for b=1:length(goodindbdlw4)
        concbdl0=squeeze(DWSpupstdd60.conctotmetrat(k,goodindbdlw4(b),:));
        crlbbdl0=squeeze(DWSpupstdd60.conctotmetrat(k,goodindbdlw4(b),:).*DWSpupstdd60.crbtotmetrat(k,goodindbdlw4(b),:));
        concbdl=concbdl0./concbdl0(1);
        aa=concbdl0;
        bb=concbdl0(1); 
        saa=crlbbdl0;
        sbb=crlbbdl0(1); 
        crlbbdl=aa./bb.*sqrt((saa./aa).^2+(sbb/bb).^2);

        kkk = ~isnan(concbdl);
        
        weighted_s_sticks = @(diff,b) 1./(crlbbdl(kkk).^2).*sqrt(pi./(4.*b.*diff)).*erf(sqrt(b.*diff));
        [pfit, resnorm, residual, exitflag, output, lambda, J] =lsqcurvefit(@(diff,b) weighted_s_sticks(diff,bval2(kkk).'),x0_sticks,bval2(kkk),1./(crlbbdl(kkk).^2).*concbdl(kkk),0,1);

        Dsticks_week4(k,goodindbdlw4(b),:)=pfit;
        
        J = full(J);
        n = numel(concbdl(kkk));        % number of data points
        p = numel(pfit);         % number of parameters
        dof = n - p;             % degrees of freedom
        
        % Estimate of error variance (mean squared residual)
        sigma2 = resnorm / dof;  % resnorm = sum(residual.^2)
        
        % Parameter covariance matrix
        Cov = sigma2 * inv(J' * J);
        Dsticks_week4_cov(k,goodindbdlw4(b),:,:)=Cov(1,1);

    end 
end 
  
week=6;

load(['DWSpupstdd60_week' num2str(week) '_qc.mat'])


bval=[50,100,500,1000,3000,5000,10000,15000];
 
bval2=bval(2:end)./1000; 
listtoplot=DWSpupstdd60.listttoplot;



%% Fit Dintra in Indiv Decays
Dsticks_week6=zeros(length(listtoplot),size(DWSpupstdd60.conctotmetrat,2),2);

for k=1:length(DWSpupstdd60.listttoplot)
%     figure; 
    for s=1:length(goodindshamw6)
        concsham0=squeeze(DWSpupstdd60.conctotmetrat(k,goodindshamw6(s),:));
        crlbsham0=squeeze(DWSpupstdd60.conctotmetrat(k,goodindshamw6(s),:).*DWSpupstdd60.crbtotmetrat(k,goodindshamw6(s),:));
        concsham=concsham0./concsham0(1);
        aa=concsham0;
        bb=concsham0(1); 
        saa=crlbsham0;
        sbb=crlbsham0(1); 
        crlbsham=aa./bb.*sqrt((saa./aa).^2+(sbb/bb).^2);
        concsham=concsham./concsham(1);
        kkk = ~isnan(concsham);

        weighted_s_sticks = @(diff,b) 1./(crlbsham(kkk).^2).*sqrt(pi./(4.*b.*diff)).*erf(sqrt(b.*diff));
        [pfit, resnorm, residual, exitflag, output, lambda, J] =lsqcurvefit(@(diff,b) weighted_s_sticks(diff,bval2(kkk).'),x0_sticks,bval2(kkk),1./(crlbsham(kkk).^2).*concsham(kkk),0,1);

        Dsticks_week6(k,goodindshamw6(s),:)=pfit;

        J = full(J);
        n = numel(concsham(kkk));        % number of data points
        p = numel(pfit);         % number of parameters
        dof = n - p;             % degrees of freedom
        
        % Estimate of error variance (mean squared residual)
        sigma2 = resnorm / dof;  % resnorm = sum(residual.^2)
        
        % Parameter covariance matrix
        Cov = sigma2 * inv(J' * J);
        Dsticks_week6_cov(k,goodindshamw6(s),:)=Cov(1,1);

    end 
    for b=1:length(goodindbdlw6)
        concbdl0=squeeze(DWSpupstdd60.conctotmetrat(k,goodindbdlw6(b),:));
        crlbbdl0=squeeze(DWSpupstdd60.conctotmetrat(k,goodindbdlw6(b),:).*DWSpupstdd60.crbtotmetrat(k,goodindbdlw6(b),:));
        concbdl=concbdl0./concbdl0(1);
        aa=concbdl0;
        bb=concbdl0(1); 
        saa=crlbbdl0;
        sbb=crlbbdl0(1); 
        crlbbdl=aa./bb.*sqrt((saa./aa).^2+(sbb/bb).^2);

        kkk = ~isnan(concbdl);
        
        weighted_s_sticks = @(diff,b) 1./(crlbbdl(kkk).^2).*sqrt(pi./(4.*b.*diff)).*erf(sqrt(b.*diff));
        [pfit, resnorm, residual, exitflag, output, lambda, J] =lsqcurvefit(@(diff,b) weighted_s_sticks(diff,bval2(kkk).'),x0_sticks,bval2(kkk),1./(crlbbdl(kkk).^2).*concbdl(kkk),0,1);
        Dsticks_week6(k,goodindbdlw6(b),:)=pfit;

        J = full(J);
        n = numel(concsham(kkk));        % number of data points
        p = numel(pfit);         % number of parameters
        dof = n - p;             % degrees of freedom
        
        % Estimate of error variance (mean squared residual)
        sigma2 = resnorm / dof;  % resnorm = sum(residual.^2)
        
        % Parameter covariance matrix
        Cov = sigma2 * inv(J' * J);
        Dsticks_week6_cov(k,goodindbdlw6(b),:)=Cov(1,1);

    end 
end 

% 
prismindiv=NaN(10,24);

%Dapp
pos=1; 
for indmet=[1,2,4,5,6,7] 
    dsticks(1:3,(pos-1)*3+1)=Dsticks_week4(indmet,goodindshamw4,1);
    dsticks(4:10,(pos-1)*3+1)=Dsticks_week6(indmet,goodindshamw6,1);
    dsticks(1:10,(pos-1)*3+2)=Dsticks_week4(indmet,goodindbdlw4,1);
    dsticks(1:6,(pos-1)*3+3)=Dsticks_week6(indmet,goodindbdlw6,1);

    dsticks_var(1:3,(pos-1)*3+1)=Dsticks_week4_cov(indmet,goodindshamw4,1);
    dsticks_var(4:10,(pos-1)*3+1)=Dsticks_week6_cov(indmet,goodindshamw6,1);
    dsticks_var(1:10,(pos-1)*3+2)=Dsticks_week4_cov(indmet,goodindbdlw4,1);
    dsticks_var(1:6,(pos-1)*3+3)=Dsticks_week6_cov(indmet,goodindbdlw6,1);

    pos=pos+1;
end 

save([pwd '/FLAMEO_input_indiv_CRLBweightedfit/dsticks.mat'],'dsticks')
save([pwd '/FLAMEO_input_indiv_CRLBweightedfit/dsticks_var.mat'],'dsticks_var')
