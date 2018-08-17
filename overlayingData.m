%% overlaying data
% f1=xlsread('20180701_cmax_halflife_LC50_forMATLAB.xlsx',1);
%[num, words, raw]=xlsread('20180701_cmax_halflife_LC50_forMATLAB.xlsx',1);
% f1=readtable('20180701_cmax_halflife_LC50_forMATLAB.xlsx');

[num, words, raw]=xlsread('20180701_cmax_halflife_LC50_forMATLAB.xlsx',1);

groupID = num(:,1);
cmax = num(:,2);
halfLife = num(:,3);
LC50_Ag = num(:,4);
LC50_Aa = num(:,5);
host = num(:,6);
drug = num(:,7);
route = num(:,8);
labelStr = words(:,9);

for ii=1:length(groupID)
%     figure(1); hold on
%     subplot(1,3,1)
%     if LC50_Aa(ii) >= 1000 && host(ii) == 8 && drug(ii) == 8 && route(ii) == 4
     if 1000 < LC50_Aa(ii) && LC50_Aa(ii) <= 1500  && drug(ii) == 8 && host(ii) == 8 && route(ii) == 5
        plot(halfLife(ii),cmax(ii),'k sq')
%         text(halfLife(ii),cmax(ii),num2str(ii))
%         str = labelStr{ii+1,1};
%         text(halfLife(ii),cmax(ii),str,'FontSize',11)
        hold on
    end        
end

for ii=1:length(groupID)
    figure(1); hold on
%     subplot(1,3,2)
    if LC50(ii)> 25/1000 && LC50(ii)<= 500/1000
        plot(halfLife(ii),cmax(ii),'k o');%text(halfLife(ii),cmax(ii),num2str(ii),'Color','r')
        hold on
    end        
end

for ii=1:length(groupID)
    figure(2); hold on
%     subplot(1,3,3)
    if LC50_Aa(ii) >= 500
        plot(halfLife(ii),cmax(ii),'k o');%text(halfLife(ii),cmax(ii),num2str(ii))
        hold on
    end        
end

for ii=1:length(groupID)
    figure(1); hold on
    if halfLife(ii) >= 10 && halfLife(ii) < 15
        text(halfLife(ii),cmax(ii),num2str(ii))
        hold on
    end        
end

for ii=1:length(groupID)
    figure(1); hold on
    if halfLife(ii) >= 15 && halfLife(ii) <30
        text(LC50(ii),cmax(ii),num2str(ii))
        hold on
    end        
end