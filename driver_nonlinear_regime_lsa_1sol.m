%% Ring fiber with 6 cores un/twisted
% Claudia Castro-Castro
% Spring 2017
% Problem  
%
% Builds numerical solution of amplitudes for twisted nonlinear mulcores
% fiber and spectral plane and stability analysis
%
% Generate .avi file instability-propagation depending on \lambda and
% nonliarity, and varying the angle \phi
%
%
% Reference : Longhi, Stefano. "PT phase control in circular multi-core fibers." 
%             Optics Letters 41, no. 9 (2016): 1897-1900.
%
% $c(1) = Re(c_1),
%  c(2) = Re(c_2),
%  c(3) = Re(c_3),
%  c(4) = Re(c_4),
%  c(5) = Re(c_5),
%  c(6) = Re(c_6),
%  c(7) = Im(c_1),
%  c(8) = Im(c_2),
%  c(9) = Im(c_3),
%  c(10) = Im(c_4),
%  c(11) = Im(c_5),
%  c(12) = Im(c_6)$
%%
close all; clear all; clc
tic

% Pre-define color of curves in plots
% https://www.mathworks.com/help/matlab/graphics_transition/why-are-plot-lines-different-colors.html
plotStyle = [      0    0.4470    0.7410; % blue
              0.8500    0.3250    0.0980; % orange 
              0.9290    0.6940    0.1250; % yellow
              0.4940    0.1840    0.5560; % purple
              0.4660    0.6740    0.1880; % green 
              0.3010    0.7450    0.9330; % sky blue
              0.6350    0.0780    0.1840  % red
              1 0 1;
              0 1 1;
              1 0 0;
              0 1 0;
              0 0 1;
              0 0 0;
              1 1 0;
              0.4 0.6 0.7;
              0.2 0.8 0.3];
          

disp('-------------------------------------------------------------------')
disp(' Numerical solution of the Amplitudes of Linear field w/ rk4')
disp('-------------------------------------------------------------------')
disp(date);

% parameters

% Number of fibers
N = 6;

% coupling strength
coupling = 1;

% nonlinearity strength 
nonlinearities = -1:0.5:5;

% propagation 
lambdas = 0:0.1:5;  % define span of lambdas in [sqrt(3), 3.1]

fun = @root4d;

% gain/loss parameter
gainloss = [ 0 , 0.1 ];

%small parameter
epsilon = 0.05;

%idenom = 1;

pi_denom = [6];

% phis = pi/6;

% phis = [ pi/9:  pi/81 :pi/3 ];

phis = [0:pi/84:11*pi/24];

for idenom = 1:length(pi_denom)
    legendInfo{idenom} = ['\phi = \pi/',num2str(pi_denom(idenom))];
end

% guess to try out
guesses = [ 3, 3, -3, -3 ];      % + sqrt



% Entries of matrix M
V0 = zeros(1,length(lambdas));
V1 = zeros(1,length(lambdas));
c0 = zeros(1,length(lambdas));
c1 = zeros(1,length(lambdas)); 

lambda1vector = 0:0.01:4;
gamma1vector = zeros(length(lambdas),length(lambda1vector)); 

cond =  zeros(1,length(lambdas));

LAMBDA = zeros(1,length(lambdas));
g = zeros(length(nonlinearities),length(lambdas));


powerA = zeros(1,length(lambdas));
powerB = zeros(1,length(lambdas));

realpart_eigvals  = zeros(4,length(nonlinearities));
imagpart_eigvals  = zeros(4,length(nonlinearities));
positive_realpart_eigvals = zeros(4,1);


Movie_data(length(phis)) = struct('cdata',[],'colormap',[]);
mov = avifile('instability_movie.avi'); 


%hfig4=figure(4); hold on 


%for idenom = 1:length(pi_denom)

visibleFlag = 'off';
set(gcf,'Visible', visibleFlag); 

for iphi = 1:length(phis)
  
    
    for inl = 1:length(nonlinearities)
    %close all
    % phi = pi/pi_denom(idenom);
    
    for iguess = 1:size(guesses,1)
    %% Order 1 solution - twisted 
    
    igl = 1 ; % gain/loss = 0
    
    phi = phis(iphi);
    
    nonlinearity = nonlinearities(inl);
    
    
    % call routine which solves nonlinear regime order 1 solution  
    [F0] = nonlinear_regime(lambdas, coupling, gainloss(igl), phi, ...
                            nonlinearity, guesses, plotStyle(mod(inl,16)+1,:),legendInfo,visibleFlag);

    %% Order \epsilon solution 
    F1 = zeros(length(lambdas), 4, size(guesses,1));
    alpha = 1.5;
    beta = -2; 

    for i=1:length(lambdas)
        
        lambda0sample=lambdas(i);

        x10 = F0(i,1,iguess); % real part of A0
        x20 = F0(i,2,iguess); % real part of B0
        y10 = F0(i,3,iguess); % imaginary part of A0
        y20 = F0(i,4,iguess); % imaginary part of B0

        V0(i)= -lambda0sample + nonlinearity*(3*(x10^2) + y10^2);
        
        V1(i)= -lambda0sample + nonlinearity*(x10^2 + 3*(y10^2));
        
        c0(i)= 2*nonlinearity*x10*y10;
        c1(i)= 2*nonlinearity*x20*y20;
        
        cond(i)= 4*nonlinearity*nonlinearity*x10*x10*y10*y10;
        
        U0 = [x10;x20;y10;y20];
        
        coef1 =(x10*(c0(i)/2)/(coupling*cos(phi)))*alpha;
        coef2 =(x20*(c1(i)/2)/(coupling*cos(phi)))*beta;
        coef3 = y10*alpha;
        coef4 = y20*beta;
        
        %fff1= @(x,y) coef1*(x-y)+coef2*(-x-y)+coef3*(x-y)+coef4*(-x+y);
        
        lambda1 = 1;
        gamma1 =  -lambda1*( (-coef1-coef2+coef3+coef4)/(coef1-coef2+coef3-coef4)  );
        
        % for last gamma0 plot gamma1 as funcion of lambda 1 for this given
        % angle
        % if i==length(lambdas)
        %    lambda1vector = 0:0.01:4;
        gamma1vector(i,:) = -lambda1vector.*( (-coef1-coef2+coef3+coef4)/(coef1-coef2+coef3-coef4)  );


        PSI = [lambda1 ,    0   ,   gamma1,     0;
                  0    , lambda1,      0   ,  - gamma1;         
              -gamma1  ,    0   ,   lambda1,     0; 
                  0    ,  gamma1,      0   ,  lambda1;];


        % Build matrix M
        MM = [-lambda0sample + nonlinearity*(3*x10^2 + y10^2), 2*coupling*cos(phi), 2*nonlinearity*x10*y10,  0  ;
               2*coupling*cos(phi), -lambda0sample + nonlinearity*(3*x20^2 + y20^2), 0 , 2*nonlinearity*x20*y20 ;
               2*nonlinearity*x10*y10, 0, -lambda0sample + nonlinearity*(x10^2 + 3*y10^2), 2*coupling*cos(phi)  ;
               0, 2*nonlinearity*x20*y20, 2*coupling*cos(phi), -lambda0sample + nonlinearity*(x20^2 + 3*y20^2)  ];

        % get pseudoinverse of M via SVD
        [U,S,V]=svd(MM);
        pS = zeros(4,4);

        % pseudo inverse of Sigma matrix
        for ii = 1:3
            pS(ii,ii) = 1/S(ii,ii);
        end

        % pseudoinverse of M
        pMM = V*(pS*U');

        % Show four properties of Pseudoinverse matrices
        % MM*(pMM*MM) == MM
        % pMM*(MM*pMM) == pMM
        % (MM*pMM)' == MM*pMM
        % (pMM*MM)' == pMM*MM

        % solution U1
        U1 = pMM*(PSI*U0);

        % store data into array to be plotted next 
        F1(i,:,iguess) = U1';
        
        x11 = F1(i,1,iguess); % real part of A1 
        x21 = F1(i,2,iguess); % real part of B1
        y11 = F1(i,3,iguess); % imaginary part of A1
        y21 = F1(i,4,iguess); % imaginary part of B1
        
        X1= x10 + epsilon*x11;
        X2= x20 + epsilon*x21;
        Y1= y10 + epsilon*y11;
        Y2= y20 + epsilon*y21;
        
        powerA(i) = X1^2 + Y1^2;
        powerB(i) = X2^2 + Y2^2;
        

        
        GAMMA =                     epsilon*gamma1;
        LAMBDA(i) = lambda0sample + epsilon*lambda1;
        
        Mlambda = [ GAMMA + 2*nonlinearity*X1*Y1, - LAMBDA(i) + nonlinearity*(X1^2 + 3*Y1^2),0, 2*coupling*cos(phi) ;
                    LAMBDA(i) - nonlinearity*(3*X1^2 + Y1^2), GAMMA - 2*nonlinearity*X1*Y1, -2*coupling*cos(phi), 0 ;
                    0, 2*coupling*cos(phi), -GAMMA + 2*nonlinearity*X2*Y2, - LAMBDA(i) + nonlinearity*(X2^2 + 3*Y2^2);
                   -2*coupling*cos(phi), 0 , LAMBDA(i) - nonlinearity*(3*X2^2 + Y2^2), -GAMMA - 2*nonlinearity*X2*Y2 ];
                
        % instability gain (real part of eigenvalue with largest positive real part
        realpart_eigvals(:,inl)  = real(eig(Mlambda));
        imagpart_eigvals(:,inl)  = imag(eig(Mlambda));
        
        positive_realpart_eigvals  = real(eig(Mlambda));
        
        for ii=1:4 
            if positive_realpart_eigvals(ii) < 0
              positive_realpart_eigvals(ii) = 0;
            end
        end
        
        g(inl,i) = max(positive_realpart_eigvals);
      

    end
    figure,
    set(gcf,'Visible', visibleFlag); 
    plot(lambdas, powerA, lambdas, powerB,'--')
    legend('|A|^2','|B|^2')
    xlabel('\lambda')



    end

    hold off

    % plot order 1 solutions dependinf on lambda^{(0)}
    plot_order1( lambdas, F1, coupling, gainloss, phi, nonlinearity, guesses, V0, V1, c0, c1, cond, plotStyle(mod(inl,16)+1,:), legendInfo, visibleFlag )

    %%



    figure, 
    set(gcf,'Visible', visibleFlag); 
    plot(realpart_eigvals(:,inl),imagpart_eigvals(:,inl),'o','Color', plotStyle(mod(inl,16)+1,:))
    ylabel('imag(\nu)'); xlabel('real(\nu)'); 
    title('Spectral plane')


 
    close all

    end % end loop over nonlinearities



hold off


%% plot instability propagation 
hformovie = figure; 
axis tight manual
waterfall(LAMBDA,nonlinearities,g)  %,'Color', plotStyle(iphi,:)
drawnow
zlabel('p'); xlabel('\lambda'); ylabel('\sigma') 
xlim([LAMBDA(1) LAMBDA(length(LAMBDA))])
axis tight manual
text(0.25,0.90,0.9, strcat('\phi = ', num2str(phis(iphi))),'Units','normalized'); 
%set(gca,'YTick',[0:pi/8:pi/2])
%set(gca,'YTickLabel', {'0', '\pi/8', '\pi/4', '3\pi/8'})
colormap gray
set(gcf,'color','white')
%view([-53 +26])

Movie_data(iphi) = getframe(hformovie);

frame = getframe(hformovie);

mov = addframe(mov,frame); 

plot_appearance(gcf,strcat('instability-propagation-phi',num2str(iphi)))



% spectral plane for last value of \lambda
figure, 
plot(realpart_eigvals,imagpart_eigvals,'o')
plot_appearance(gcf,strcat('spectral-plane-nl',num2str(nonlinearity)))

end % end loop over phis

mov = close(mov);
toc 

%% store instability data on disk 
save('MovieStruc.mat','Movie_data')


% Create AVI file on disk.
movie2avi(Movie_data, 'instability.avi', 'compression', 'None');



