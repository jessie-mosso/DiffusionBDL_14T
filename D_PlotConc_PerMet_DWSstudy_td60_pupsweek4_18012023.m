%Read batch excel file, extract concentrations and plot

clear; clc;close all;

%% to adapt
indratqc=[1:15];%[1:9,12:15];
% indratqc=[1:15];
dosaveall=true;
bmin=2;

%%
datesw4={'2021-03-26','2021-03-22','2021-03-31','2021-03-29','2021-03-29','2021-03-31',...
    '2021-09-30','2021-09-30','2021-10-01','2022-08-17','2022-08-15','2022-08-15','2022-08-16',...
    '2022-08-18','2022-08-19'};
timeweek4={'','','am','am','pm','pm',...
    'am','pm','',...
    '','am','pm','','',''};

ratnames={['SHAM576'],['SHAM577'],['SHAM578'],...
    ['BDL579'],['BDL580'],['BDL581'],['BDL610'],['BDL611'],['BDL612'],...
    ['SHAM659'],['BDL667'],['BDL668'],['BDL669'],['BDL670'],['BDL671']};

week=4;


iii=1;
for rat=1:length(ratnames)
    ratlist(rat)=convertCharsToStrings([ratnames{rat} '_' datesw4{rat}]);
end 

bval=[50,100,500,1000,3000,5000,10000,15000];
 
% bval2=bval./1000; %remove b500
bval2=bval(bmin:end)./1000; %remove b500


listtoplot=["Gln","Glu","Glu+Gln","Ins","Tau","NAA+NAAG","Cr+PCr","GPC+PCho","NAA","NAAG","Cr","PCr"];
monotonicitytable=zeros(length(ratnames),length(listtoplot));

uplim=1.4;
downlim=0;
iii=1;
figure;
for met=1:length(listtoplot)%[4]%[1,2,4,6]%
    disp(listtoplot(met))
    
    
    kk=1;
    ll=1;
    
    subplot(2,6,met)
    % subplot(1,2,iii)
    iii=iii+1;
    for rat=indratqc%1:length(ratlist)
        ratname=char(ratlist(rat));

        excelsheet=readtable(['DWSpups_stelaser_cereb_week' num2str(week) '_allgroups_01032023.xlsx'],'Sheet',ratlist(rat));
 
        metlist=excelsheet.Metabolite;

        mettofind=strcmp(metlist,char(listtoplot(met)));
        indtot=1:1:40;
        indmet=indtot(mettofind==1);

        clear conctot1met;
        clear crbtot1met;
        
        %get conc for 1 met, 1 rat
        ii=1;
        for b=bmin:length(bval)%[1:length(bval2)]%remove b500
            conctot1met(ii)=excelsheet{indmet,1+(b-1)*5};
            crbtot1met(ii)=excelsheet{indmet,2+(b-1)*5};
            ii=ii+1;
        end
        
        conctotmetrat(met,rat,:)=conctot1met;
        crbtotmetrat(met,rat,:)=crbtot1met; 
        
        %abscrb
        abscrbtot1met=crbtot1met.*conctot1met;
        
        %normalized conc and crb         
        normconc=conctot1met./conctot1met(1); 
        normcrb=normconc.*sqrt((abscrbtot1met./conctot1met).^2+(abscrbtot1met(1)./conctot1met(1)).^2);
        
        % if rat==11
        %     h1=errorbar(bval2,normconc-0.2*rat,normcrb,'-','color','black');%,'DisplayName',ratname);
        % elseif rat==12
        %     h1=errorbar(bval2,normconc-0.2*rat,normcrb,'-','color','black');%,'DisplayName',ratname);
        if ratname(1)=='B'
            
            h1=errorbar(bval2,normconc,normcrb,'-','DisplayName',ratname);
        else
            % h1=errorbar(bval2,normconc,normcrb,'--','DisplayName',ratname);%-0.2*rat
        end 
        monotonicitytable(rat,met)=issorted(normconc(1:end),'descend');
        xlim([-1,bval2(end)+2])
        % ylim([downlim,uplim])
        hold on
        set(gca, 'box','off')
        title([char(listtoplot(met))])
        set(gca,'FontSize',20)    
        xlabel('bval')

    end 
end
%%
if dosaveall
    DWSpupstdd60.listttoplot=listtoplot;
    DWSpupstdd60.ratnames=ratnames(indratqc);
    DWSpupstdd60.bval2=bval2;
%     DWSadultstd60.indexesuse=[3:7,11,12];
    DWSpupstdd60.conctotmetrat=conctotmetrat(:,indratqc,:);
    DWSpupstdd60.crbtotmetrat=crbtotmetrat(:,indratqc,:);
    save('DWSpupstdd60_week4_qc_moremet.mat','DWSpupstdd60')
end 
