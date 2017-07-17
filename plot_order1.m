function plot_order1( lambdas, Order1Sol, coupling, gainloss, phi, nonlinearity, guesses, V0, V1, c0, c1, cond, plotStyle,legendInfo,visibleFlag  )
 

% set Phase tag 
if phi == 0
    phi_tag = '\phi = 0';
else
    phi_tag = '\phi = \pi/6';
end

% set gainloss tag 
gamma = gainloss;
strgl = strcat('-gainloss-',num2str(gamma(2)));
gamma_tag = strcat('\gamma = ',num2str(gamma(2)));


% set Phase tag 
if phi == 0
    %phi_tag = '\phi = 0';
    phi_filename = '-untwisted';
    fig_titlestr = 'Untwisted';
else
    %phi_tag = strcat('\phi = \pi/',num2str(phi_denom));
    phi_filename = strcat('-twisted-pis');%,num2str(pi_denom));
    fig_titlestr = strcat('Twisted \phi = ',num2str(phi));
end



%h1 = figure(5);  hAx1 = axes;  hold on

figure, 
set(gcf,'Visible', visibleFlag); 

h1=subplot(2,2,1); hold on
for iguess=1:size(guesses,1)
    plot(h1,lambdas,Order1Sol(:,1,iguess),'.','Color', plotStyle)
    % set(h1, 'Tag','ReA');
end
ylabel('ReA^{(1)}'); xlabel('\lambda^{(0)}'); %grid on
%legend(legendInfo,'Location','SouthWest'); legend boxoff
xlim([lambdas(1) lambdas(length(lambdas))])
text(0.01,0.90,'a)','Units','normalized'); 
set(gca,'FontSize',14) % Axis thicks fontsize
hold off
%plot_appearance(h1,strcat('stationary-sol-ReA1','-gamma1',phi_filename))

%h2 = figure(6);  hAx2 = axes;  hold on   

h2=subplot(2,2,2); hold on
for iguess=1:size(guesses,1)
    plot(h2,lambdas,Order1Sol(:,3,iguess),'.','Color', plotStyle)
end
ylabel('ImA^{(1)}'); xlabel('\lambda^{(0)}'); %grid on
%legend(legendInfo,'Location','SouthWest'); legend boxoff
xlim([lambdas(1) lambdas(length(lambdas))])
text(0.01,0.90,'b)','Units','normalized'); 
set(gca,'FontSize',14) % Axis thicks fontsize
hold off
%plot_appearance(h2,strcat('stationary-sol-ImA1','-gamma1',phi_filename))


%h3 = figure(7);  hAx3 = axes;  hold on  

h3=subplot(2,2,3); hold on
for iguess=1:size(guesses,1)
    plot(h3,lambdas,Order1Sol(:,2,iguess),'.','Color', plotStyle)
end
ylabel('ReB^{(1)}'); xlabel('\lambda^{(0)}'); %grid on 
%legend(legendInfo,'Location','SouthWest');legend boxoff
xlim([lambdas(1) lambdas(length(lambdas))])
set(gca,'FontSize',14) % Axis thicks fontsize
text(0.01,0.90,'c)','Units','normalized'); 
hold off
%plot_appearance(h3,strcat('stationary-sol-ReB1','-gamma1',phi_filename))

%h4 = figure(8);  hAx4 = axes;  hold on 
h4=subplot(2,2,4); hold on
for iguess=1:size(guesses,1)
    plot(h4,lambdas,Order1Sol(:,4,iguess),'.','Color', plotStyle) 
end
ylabel('ImB^{(1)}'); xlabel('\lambda^{(0)}'); %grid on
%legend(legendInfo,'Location','SouthWest');legend boxoff
xlim([lambdas(1) lambdas(length(lambdas))])
set(gca,'FontSize',14) % Axis thicks fontsize
text(0.01,0.90,'d)','Units','normalized'); 
hold off

hold off

plot_appearance(gcf,strcat('stationary-sols-order1','-gamma1',phi_filename))


% %% Graph condition to have all pivots
% cosphi2 = 4*coupling*coupling*(cos(phi))^2*ones(1,length(lambdas));
% 
% figure(3); 
% hold on
% plot(lambdas,cond,'Color', plotStyle)
% plot(lambdas,cosphi2,'.','Color', plotStyle)
% %ylabel('');
% xlabel('\lambda^{(0)}'); 
% set(gca,'FontSize',14) % Axis thicks fontsize
% %text(0.01,0.90,'a)','Units','normalized'); 
% legend('4\sigma^2x_1^{(0)2}y_1^{(0)2}','4k^2cos^2(\phi)','Location','NorthEast')
% xlim([lambdas(1)-.2 lambdas(length(lambdas))])
% legend boxoff
% grid on
% box off
% hold off
% plot_appearance(gcf,'FullRankCondition')






end