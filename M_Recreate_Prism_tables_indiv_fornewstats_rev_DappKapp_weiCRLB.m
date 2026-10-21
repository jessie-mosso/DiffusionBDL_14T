clear;clc;close all; 

goodindshamw4=[1:3];
goodindbdlw4=[4:13];


goodindshamw6=[1:3,6,9:11];
goodindbdlw6=[4:5,7:8,12:13];


week=4;

load(['DWSpupstdd60_week' num2str(week) '_qc.mat'])


bval=[50,100,500,1000,3000,5000,10000,15000];
bmax=15; 

if bmax==10
    bval2=bval(2:end-1)./1000; 
else
    bval2=bval(2:end)./1000; 
end 
listtoplot=DWSpupstdd60.listttoplot;

ll=[0.05;0.5];
ul=[0.2;3];
x0_dk=[0.1,1]; %inital guess for diff

%% Fit Dintra in Indiv Decays
DappKapp_week4=zeros(length(listtoplot),size(DWSpupstdd60.conctotmetrat,2),2);
figure; 
for k=1:length(DWSpupstdd60.listttoplot)
%     figure; 
    for s=1:length(goodindshamw4)
        
        concsham0=squeeze(DWSpupstdd60.conctotmetrat(k,goodindshamw4(s),:));
        crlbsham0=squeeze(DWSpupstdd60.conctotmetrat(k,goodindshamw4(s),:).*DWSpupstdd60.crbtotmetrat(k,goodindshamw4(s),:));
        if bmax==10
            concsham0=concsham0(1:end-1);
            crlbsham0=crlbsham0(1:end-1);
        end

        concsham=concsham0./concsham0(1);
        aa=concsham0;
        bb=concsham0(1); 
        saa=crlbsham0;
        sbb=crlbsham0(1); 
        crlbsham=aa./bb.*sqrt((saa./aa).^2+(sbb/bb).^2);
        
        
        kkk = ~isnan(concsham);
        weighted_s_dappkapp = @(diff,b) 1./(crlbsham(kkk).^2).*exp(-b.*diff(1)+b.^2.*diff(1).^2.*diff(2)/6); %
        [pfit, resnorm, residual, exitflag, output, lambda, J] = lsqcurvefit(@(diff,b) weighted_s_dappkapp(diff,bval2(kkk).'),x0_dk,bval2(kkk),1./(crlbsham(kkk).^2).*concsham(kkk),ll,ul);

        DappKapp_week4(k,goodindshamw4(s),:)=pfit;


        errorbar(bval2(kkk),concsham(kkk),crlbsham(kkk))
        hold on 
        plot(bval2(kkk),exp(-bval2(kkk).*pfit(1)+bval2(kkk).^2.*pfit(1).^2.*pfit(2)/6),'LineWidth',2,'color','black')
        hold on 

        J = full(J);
        n = numel(concsham(kkk));        % number of data points
        p = numel(pfit);         % number of parameters
        dof = n - p;             % degrees of freedom
        
        % Estimate of error variance (mean squared residual)
        sigma2 = resnorm / dof;  % resnorm = sum(residual.^2)
        
        % Parameter covariance matrix
        Cov = sigma2 * inv(J' * J);
        DappKapp_week4_cov(k,goodindshamw4(s),:)=[Cov(1,1),Cov(2,2)];

        %Weighted linear LS
        % Design matrix: columns are [b, b.^2], no intercept
        X = [bval2(kkk).',bval2(kkk).'.^2];
        logy = log(concsham(kkk));
        
        % Linear least squares
%         c = X \ logy;     % c(1) = coeff of b, c(2) = coeff of b^2
        
        %weighted LLS 
        w = concsham(kkk).^2;              % approximate weights
        W = diag(w);
        c = (X' * W * X) \ (X' * W * logy);


        c1 = c(1);        % = -d1
        c2 = c(2);        % = d1^2 * d2 / 6

        d1 = -c1;
        d2 = 6 * c2 / d1^2;
        DappKapp_week4lls(k,goodindshamw4(s),:)=[d1,d2];
%         plot(bval2,concsham,'-o')
%         hold on 
%         plot(bval2,s_dappkapp(x_sol_dappkapp,bval2),'--');
%         hold on 
    end 
    for b=1:length(goodindbdlw4)
        concbdl0=squeeze(DWSpupstdd60.conctotmetrat(k,goodindbdlw4(b),:));
        crlbbdl0=squeeze(DWSpupstdd60.conctotmetrat(k,goodindbdlw4(b),:).*DWSpupstdd60.crbtotmetrat(k,goodindbdlw4(b),:));
        if bmax==10
            concbdl0=concbdl0(1:end-1);
            crlbbdl0=crlbbdl0(1:end-1);
        end

        concbdl=concbdl0./concbdl0(1);
        aa=concbdl0;
        bb=concbdl0(1); 
        saa=crlbbdl0;
        sbb=crlbbdl0(1); 
        crlbbdl=aa./bb.*sqrt((saa./aa).^2+(sbb/bb).^2);

        kkk = ~isnan(concbdl);
        weighted_s_dappkapp = @(diff,b) 1./(crlbbdl(kkk).^2).*exp(-b.*diff(1)+b.^2.*diff(1).^2.*diff(2)/6); %
        [pfit, resnorm, residual, exitflag, output, lambda, J] = lsqcurvefit(@(diff,b) weighted_s_dappkapp(diff,bval2(kkk).'),x0_dk,bval2(kkk),1./(crlbbdl(kkk).^2).*concbdl(kkk),ll,ul);


        DappKapp_week4(k,goodindbdlw4(b),:)=pfit;

        errorbar(bval2(kkk),concbdl(kkk),crlbbdl(kkk))
        hold on 
        plot(bval2(kkk),exp(-bval2(kkk).*pfit(1)+bval2(kkk).^2.*pfit(1).^2.*pfit(2)/6),'LineWidth',2,'color','black')
        hold on 
        J = full(J);
        n = numel(concbdl(kkk));        % number of data points
        p = numel(pfit);         % number of parameters
        dof = n - p;             % degrees of freedom
        
        % Estimate of error variance (mean squared residual)
        sigma2 = resnorm / dof;  % resnorm = sum(residual.^2)
        
        % Parameter covariance matrix
        Cov = sigma2 * inv(J' * J);
        DappKapp_week4_cov(k,goodindbdlw4(b),:,:)=[Cov(1,1),Cov(2,2)];

         %Weighted linear LS
        % Design matrix: columns are [b, b.^2], no intercept
        X = [bval2(kkk).',bval2(kkk).'.^2];
        logy = log(concbdl(kkk));
        
        % Linear least squares
%         c = X \ logy;     % c(1) = coeff of b, c(2) = coeff of b^2
        
        %weighted LLS 
        w = concbdl(kkk).^2;              % approximate weights
        W = diag(w);
        c = (X' * W * X) \ (X' * W * logy);


        c1 = c(1);        % = -d1
        c2 = c(2);        % = d1^2 * d2 / 6

        d1 = -c1;
        d2 = 6 * c2 / d1^2;
        DappKapp_week4lls(k,goodindbdlw4(b),:)=[d1,d2];
%         plot(bval2,concbdl,'-*')
%         hold on 
%         plot(bval2,s_dappkapp(x_sol_dappkapp,bval2),'--');
%         hold on 
    end 
end 
  

figure; 
for kk=1:2
subplot(1,2,kk)
for k=[1,2,4:7]
    h=plot(DappKapp_week4(k,:,kk),'-')
    cc=get(h,'color'); 
    hold on 
%     plot(DappKapp_week4lls(k,:,kk),'--','color',cc)
%     hold on 
end 
end 


week=6;

load(['DWSpupstdd60_week' num2str(week) '_qc.mat'])


bval=[50,100,500,1000,3000,5000,10000,15000];
 
bval2=bval(2:end)./1000; 
listtoplot=DWSpupstdd60.listttoplot;


x0_dk=[0.1,1]; %inital guess for diff

figure; 

%% Fit Dintra in Indiv Decays
DappKapp_week6=zeros(length(listtoplot),size(DWSpupstdd60.conctotmetrat,2),2);

for k=1:length(DWSpupstdd60.listttoplot)
%     figure; 
    for s=1:length(goodindshamw6)

        concsham0=squeeze(DWSpupstdd60.conctotmetrat(k,goodindshamw6(s),:));
        crlbsham0=squeeze(DWSpupstdd60.conctotmetrat(k,goodindshamw6(s),:).*DWSpupstdd60.crbtotmetrat(k,goodindshamw6(s),:));
        if bmax==10
            concsham0=concsham0(1:end-1);
            crlbsham0=crlbsham0(1:end-1);
        end

        concsham=concsham0./concsham0(1);
        aa=concsham0;
        bb=concsham0(1); 
        saa=crlbsham0;
        sbb=crlbsham0(1); 
        crlbsham=aa./bb.*sqrt((saa./aa).^2+(sbb/bb).^2);
        
        kkk = ~isnan(concsham);
        weighted_s_dappkapp = @(diff,b) 1./(crlbsham(kkk).^2).*exp(-b.*diff(1)+b.^2.*diff(1).^2.*diff(2)/6); %
        [pfit, resnorm, residual, exitflag, output, lambda, J] = lsqcurvefit(@(diff,b) weighted_s_dappkapp(diff,bval2(kkk).'),x0_dk,bval2(kkk),1./(crlbsham(kkk).^2).*concsham(kkk),ll,ul);

        DappKapp_week6(k,goodindshamw6(s),:)=pfit;

        errorbar(bval2(kkk),concsham(kkk),crlbsham(kkk))
        hold on 
        plot(bval2(kkk),exp(-bval2(kkk).*pfit(1)+bval2(kkk).^2.*pfit(1).^2.*pfit(2)/6),'LineWidth',2,'color','black')
        hold on 

        J = full(J);
        n = numel(concsham(kkk));        % number of data points
        p = numel(pfit);         % number of parameters
        dof = n - p;             % degrees of freedom
        
        % Estimate of error variance (mean squared residual)
        sigma2 = resnorm / dof;  % resnorm = sum(residual.^2)
        
        % Parameter covariance matrix
        Cov = sigma2 * inv(J' * J);
        DappKapp_week6_cov(k,goodindshamw6(s),:)=[Cov(1,1),Cov(2,2)];

    end 
    for b=1:length(goodindbdlw6)
        concbdl0=squeeze(DWSpupstdd60.conctotmetrat(k,goodindbdlw6(b),:));
        crlbbdl0=squeeze(DWSpupstdd60.conctotmetrat(k,goodindbdlw6(b),:).*DWSpupstdd60.crbtotmetrat(k,goodindbdlw6(b),:));
        if bmax==10
            concbdl0=concbdl0(1:end-1);
            crlbbdl0=crlbbdl0(1:end-1);
        end

        concbdl=concbdl0./concbdl0(1);
        aa=concbdl0;
        bb=concbdl0(1); 
        saa=crlbbdl0;
        sbb=crlbbdl0(1); 
        crlbbdl=aa./bb.*sqrt((saa./aa).^2+(sbb/bb).^2);

        kkk = ~isnan(concbdl);
        weighted_s_dappkapp = @(diff,b) 1./(crlbbdl(kkk).^2).*exp(-b.*diff(1)+b.^2.*diff(1).^2.*diff(2)/6); %
        [pfit, resnorm, residual, exitflag, output, lambda, J] = lsqcurvefit(@(diff,b) weighted_s_dappkapp(diff,bval2(kkk).'),x0_dk,bval2(kkk),1./(crlbbdl(kkk).^2).*concbdl(kkk),ll,ul);


        DappKapp_week6(k,goodindbdlw6(b),:)=pfit;
        
        errorbar(bval2(kkk),concbdl(kkk),crlbbdl(kkk))
        hold on 
        plot(bval2(kkk),exp(-bval2(kkk).*pfit(1)+bval2(kkk).^2.*pfit(1).^2.*pfit(2)/6),'LineWidth',2,'color','black')
        hold on 
      

        J = full(J);
        n = numel(concsham(kkk));        % number of data points
        p = numel(pfit);         % number of parameters
        dof = n - p;             % degrees of freedom
        
        % Estimate of error variance (mean squared residual)
        sigma2 = resnorm / dof;  % resnorm = sum(residual.^2)
        
        % Parameter covariance matrix
        Cov = sigma2 * inv(J' * J);
        DappKapp_week6_cov(k,goodindbdlw6(b),:)=[Cov(1,1),Cov(2,2)];

    end 
end 


prismindiv=NaN(10,24);

%Dapp
pos=1; 
for indmet=[1,2,4,5,6,7] 
    dapp(1:3,(pos-1)*3+1)=DappKapp_week4(indmet,goodindshamw4,1);
    dapp(4:10,(pos-1)*3+1)=DappKapp_week6(indmet,goodindshamw6,1);
    dapp(1:10,(pos-1)*3+2)=DappKapp_week4(indmet,goodindbdlw4,1);
    dapp(1:6,(pos-1)*3+3)=DappKapp_week6(indmet,goodindbdlw6,1);

    dapp_var(1:3,(pos-1)*3+1)=DappKapp_week4_cov(indmet,goodindshamw4,1);
    dapp_var(4:10,(pos-1)*3+1)=DappKapp_week6_cov(indmet,goodindshamw6,1);
    dapp_var(1:10,(pos-1)*3+2)=DappKapp_week4_cov(indmet,goodindbdlw4,1);
    dapp_var(1:6,(pos-1)*3+3)=DappKapp_week6_cov(indmet,goodindbdlw6,1);

    kapp(1:3,(pos-1)*3+1)=DappKapp_week4(indmet,goodindshamw4,2);
    kapp(4:10,(pos-1)*3+1)=DappKapp_week6(indmet,goodindshamw6,2);
    kapp(1:10,(pos-1)*3+2)=DappKapp_week4(indmet,goodindbdlw4,2);
    kapp(1:6,(pos-1)*3+3)=DappKapp_week6(indmet,goodindbdlw6,2);

    kapp_var(1:3,(pos-1)*3+1)=DappKapp_week4_cov(indmet,goodindshamw4,2);
    kapp_var(4:10,(pos-1)*3+1)=DappKapp_week6_cov(indmet,goodindshamw6,2);
    kapp_var(1:10,(pos-1)*3+2)=DappKapp_week4_cov(indmet,goodindbdlw4,2);
    kapp_var(1:6,(pos-1)*3+3)=DappKapp_week6_cov(indmet,goodindbdlw6,2);
    pos=pos+1;
end 

save([pwd '/FLAMEO_input_indiv_CRLBweightedfit/dapp.mat'],'dapp')
save([pwd '/FLAMEO_input_indiv_CRLBweightedfit/dapp_var.mat'],'dapp_var')

save([pwd '/FLAMEO_input_indiv_CRLBweightedfit/kapp.mat'],'kapp')
save([pwd '/FLAMEO_input_indiv_CRLBweightedfit/kapp_var.mat'],'kapp_var')