function [F] = nonlinear_regime( lambdas, coupling, gamma, phi, nonlinearity, x0 , plotStyle, legendInfo,visibleFlag )
% usage [F] =  nonlinear_regime( lambdas, coupling, gamma, phi, nonlinearity, x0 , plotStyle, legendInfo,visibleFlag )

% solutions of order 1 in the nonlinear regime

% set Phase tag 
if phi == 0
    phi_tag = '\phi = 0';
else
    phi_tag = strcat('\phi = ',num2str(phi));
end

% set gainloss tag 
strgl = strcat('-gainloss-',num2str(gamma));
gamma_tag = strcat('\gamma = ',num2str(gamma));


% set Phase tag 
if phi == 0
    %phi_tag = '\phi = 0';
    phi_filename = '-untwisted';
    fig_titlestr = 'Untwisted';
else
    %phi_tag = strcat('\phi = \pi/',num2str(phi_denom));
    phi_filename = strcat('-twisted-pi-',num2str(phi));
    fig_titlestr = strcat('Twisted \phi = ', num2str(phi));
end

options = optimset('Display','off');

figure, hold on
set(gcf,'Visible', visibleFlag); 

ReA_order0 = sqrt(lambdas - sqrt(3))./(2*nonlinearity); 

F = zeros(length(lambdas), 4 , size(x0,1));

% loop over however many initial guesses  
for iguess = 1:size(x0,1)
    
    for i = 1:length(lambdas)

        [F(i,:,iguess),fval] = fsolve(@(x) root4d(x,lambdas(i), coupling, gamma, phi, nonlinearity), x0(iguess,:), options);

    end

    % lambda zero 
    [F0,fval] = fsolve(@(x) root4d(x,0, coupling,gamma,phi, nonlinearity), x0(iguess,:), options);

    set(0,'defaultlinelinewidth',2) % line thickness everywhere
    
    h1=subplot(2,2,1); 
    
    if iguess == 1
        % h1=figure(1);
        %hAx1 = axes;
        %hold on 
        %plot(hAx1,lambdas,F(:,1,iguess),'-.','Color', plotStyle(iguess,:))
        plot(lambdas,F(:,1,iguess),'.','Color', plotStyle)
        set(h1,'Tag','ReA');
        ylabel('ReA^{(0)}'); xlabel('\lambda^{(0)}'); %grid on
        text(.7,.6,{gamma_tag},'Units','normalized')
        %legend(legendInfo,'Location','SouthWest')
        %legend boxoff
        %set(gca,'XTickLabel',lambdas(1):1:lambdas(length(lambdas))) 
        xlim([lambdas(1)-.2 lambdas(length(lambdas))])
        set(gca,'FontSize',14) % Axis thicks fontsize
        text(0.01,0.90,'a)', 'Units','normalized'); 
    end 
    
    if iguess > 1
       h1 = findobj('Tag','ReA'); % get handle to object tagged as 'left'
       set(h1,'Nextplot','add')
       plot(h1,lambdas,F(:,1,iguess),'.','Color', plotStyle)
       %set(gca,'XTickLabel',lambdas(1):2:lambdas(length(lambdas)))  
    end   
    
    h2=subplot(2,2,2); 
    
    if iguess == 1
        %h2=figure(2);
        %hAx2 = axes;
        %hold on 
        plot(lambdas,F(:,3,iguess),'.','Color', plotStyle)
        set(h2,'Tag','ImA');
        ylabel('ImA^{(0)}'); xlabel('\lambda^{(0)}'); %grid on
        text(.7,.6,{gamma_tag},'Units','normalized')
        %legend(legendInfo,'Location','SouthWest')
        %legend boxoff
        %set(gca,'XTickLabel',lambdas(1):1:lambdas(length(lambdas)))
        xlim([lambdas(1)-.2 lambdas(length(lambdas))])
        set(gca,'FontSize',14) % Axis thicks fontsize
        text(0.01,0.90,'b)','Units','normalized'); 
    end 
    
    if iguess > 1
        h2 = findobj('Tag','ImA'); % get handle to object tagged as 'left'
        set(h2,'Nextplot','add')
        plot(h2,lambdas,F(:,3,iguess),'.','Color', plotStyle)
        %set(gca,'XTickLabel',lambdas(1):1:lambdas(length(lambdas)))  
    end
    
    h3=subplot(2,2,3);
    
    if iguess==1 
        %h3 = figure(3);
        %hAx3 =axes;
        %hold on
        plot(lambdas,F(:,2,iguess),'.','Color', plotStyle)
        set(h3,'Tag','ReB');
        ylabel('ReB^{(0)}'); xlabel('\lambda^{(0)}'); %grid on 
        text(.7,.6,{gamma_tag},'Units','normalized')
        %legend(legendInfo,'Location','SouthWest')
        %legend boxoff
        %set(gca,'XTickLabel',lambdas(1):1:lambdas(length(lambdas)))
        xlim([lambdas(1)-.2 lambdas(length(lambdas))])
        set(gca,'FontSize',14) % Axis thicks fontsize
        text(0.01,0.90,'c)','Units','normalized'); 
        
    end
    
    if iguess > 1
        h3 = findobj('Tag','ReB'); % get handle to object tagged as 'left'
        set(h3,'Nextplot','add')
        plot(h3,lambdas,F(:,2,iguess),'.','Color', plotStyle)
        %set(gca,'XTickLabel',lambdas(1):1:lambdas(length(lambdas)))  
    end

    h4=subplot(2,2,4);
    %h4 = figure(4), hold on
    if iguess==1
        plot(lambdas,F(:,4,iguess),'.','Color', plotStyle)
        set(h4,'Tag','ImB');
        ylabel('ImB^{(0)}'); xlabel('\lambda^{(0)}'); %grid on
        text(.7,.6,{gamma_tag},'Units','normalized')
        %legend(legendInfo,'Location','SouthWest')
        %legend boxoff
        %set(gca,'XTickLabel',lambdas(1):1:lambdas(length(lambdas)))  
        xlim([lambdas(1)-.2 lambdas(length(lambdas))])
        set(gca,'FontSize',14) % Axis thicks fontsize
        text(0.01,0.90,'d)','Units','normalized'); 
    end
    if iguess > 1
        h4 = findobj('Tag','ImB'); % get handle to object tagged as 'left'
        set(h4,'Nextplot','add')
        plot(h4,lambdas,F(:,4,iguess),'.','Color', plotStyle)
        %set(gca,'XTickLabel',lambdas(1):1:lambdas(length(lambdas)))  
    end
    
end
    
    %suptitle(fig_titlestr)
    hold off
    plot_appearance(gcf,strcat('stationary-sols-order0',strgl, phi_filename))
end