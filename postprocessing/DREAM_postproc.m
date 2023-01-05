

% ------------------------------------------------------------------------------------------------------------
% ------------------------------------------------------------------------------------------------------------
% ------------------------------------------------ PRE-PROCESSING -------------------------------------------
% ------------------------------------------------------------------------------------------------------------
% ------------------------------------------------------------------------------------------------------------

% First assemble all chain in one matrix
ParSet = GenParSet(chain); DREAMPar.N = size(chain,3);

% Take the last 25% of the posterior samples -- assume that these samples
% are posterior samples (double check that R_stat < 1.2 for all parameters)
Pars = ParSet ( floor ( 0.75 * size(ParSet,1) ) : size(ParSet,1), 1 : DREAMPar.d );

% How many posterior parameter samples?
N_Pars = size(Pars,1);

% Define function handle
f_handle = eval(['@(x)',char(Func_name),'(x)']);

% If not now check whether model produces simulation or not?
sim_out = [];

% Now see whether we use real data or not
if exist('Meas_info')
    % Check if field Y exists
    if isfield(Meas_info,'Y'),
        % How many elements does Meas_info have?
        Meas_info.N = size(Meas_info.Y,1);
        % Did we store runs in fx or not?
        if isfield(DREAMPar,'modout'),
            % If yes, then take from fx
            sim_out = fx ( floor ( 0.75 * size(fx,1) ) : size(fx,1), 1 : Meas_info.N );
        else
            sim_out = NaN ( N_Pars , Meas_info.N );
            % Initialize waitbar
            h = waitbar(0,'Running posterior simulations - Please wait...');
            % Loop over each sample
            for qq = 1 : N_Pars,
                sim_out(qq,1:Meas_info.N) = f_handle(Pars(qq,1:DREAMPar.d));
                % Update the waitbar --> to see simulation progress on screen
                waitbar(qq/N_Pars,h);
            end;
            % Now close waitbar
            close(h);
        end
    end
else
    % No simulations and no summary metrics
    fx_post = []; FX_post = [];
end;

% Not ABC
if ~isfield(DREAMPar,'ABC');
    % And sim_out exists
    if ~isempty(sim_out),
        % must be posterior simulations
        fx_post = sim_out; FX_post = [];
    end
elseif isfield(DREAMPar,'ABC');
    % If field "S" of Meas_info exists --> summary metrics as prior
    if isfield(Meas_info,'S'),
        % sim_out are model simulations
        fx_post = sim_out;
        % Now compute summary statistics from fx (model simulations)
        h = waitbar(0,'Calculating posterior summary metrics - Please wait...');
        for qq = 1 : N_Pars,
            FX_post(qq,:) = DREAMPar.prior_handle(fx_post(qq,:));
            % Update the waitbar --> to see simulation progress on screen
            waitbar(qq/N_Pars,h);
        end;
        % Now close waitbar
        close(h);
        
    else
        % Field "S" of Meas_info does not exist --> summary metrics as likelihood
        fx_post = []; FX_post = sim_out;
    end;
end;

% Now determine the size of fx_post (columns is number observations)
Meas_info.N = size(fx_post,2);

% ------------------------------------------------------------------------------------------------------------
% ------------------------------------------------------------------------------------------------------------
% -------------------------------------------- END OF PRE-PROCESSING -----------------------------------------
% ------------------------------------------------------------------------------------------------------------
% ------------------------------------------------------------------------------------------------------------

% Find the maximum aposteriori parameter values (last column of ParSet are log-density values!)
[~,idx] = max(ParSet(:,end)); idx = idx(1);

% Print those to screen
MAP = ParSet(idx,1:DREAMPar.d)

% Calculate the mean posterior value of each parameter
MEAN = mean(Pars)

% Calculate the posterior standard deviation of the parameters
STD = std(Pars)

% Calculate the DREAMPar.d-dimensional parameter correlation matrix (R-values)
CORR = corrcoef(Pars)

% Set figure number
fig_number = 1;

% ------------------------------------------------------------------------------------------------------------
% -------------------------------- EVOLUTION OF R_STATISTIC OF GELMAN AND RUBIN ------------------------------
% ------------------------------------------------------------------------------------------------------------

% Now plot the R_statistic for each parameter
figure(fig_number),
% Update figure number
fig_number = fig_number + 1;
% First print the R-statistic of Gelman and Rubin (each parameter a different color)
plot(output.R_stat(:,1),output.R_stat(:,2:DREAMPar.d+1)); hold on;
% Add labels
xlabel('迭代的种群数目','fontsize',14,'fontname','Times');
ylabel('比例缩小因子 R_{stat}','fontsize',14,'fontname','Times');
% Add title
title('采样的链收敛情况','fontsize',14,'fontname','Times');
% Now add the theoretical convergence value of 1.2 as horizontal line
plot([0 output.R_stat(end,1)],[1.2 1.2],'k--','linewidth',2);
% Set the the axes
axis([0 output.R_stat(end,1) 0.8 5]);
% Add a legend
evalstr = strcat('legend(''参数.1''');
% Each parameter a different color
for j = 2:DREAMPar.d,
    % Add the other parameters
    evalstr = strcat(evalstr,',''参数. ',num2str(j),'''');
end;
% And now conclude with a closing bracket
evalstr = strcat(evalstr,');');
% Now evaluate the legend
eval(evalstr);

% ------------------------------------------------------------------------------------------------------------
% -------------------------------- HISTOGRAMS OF MARGINAL DENSITIES OF PARAMETERS ----------------------------
% ------------------------------------------------------------------------------------------------------------

% Plot the histograms (marginal density) of each parameter;
% What lay out of marginal distributions is desired subplot(r,t)
% r = 3; t = 2;
% % How many figures do we need to create with this layout?
% N_fig = ceil( DREAMPar.d / (r * t) ); counter = 1; j = 1;
% % Open new figure
% figure(fig_number);
% % Now plot each parameter
% while counter <= DREAMPar.d
%     % Check whether to open a new figure?
%     if j == (r * t) + 1,
%         % Update fig_number
%         fig_number = fig_number + 1;
%         % Open new figure
%         figure(fig_number);
%         % Reset j to 1
%         j = 1;
%     end;
%     % Now create histogram
%     [N,X] = hist(Pars(:,counter));
%     % And plot histogram in red
%     subplot(r,t,j),bar(X,N/sum(N),'r'); hold on; % --> can be scaled to 1 if using "trapz(X,N)" instead of "sum(N)"!
%     if j == 1,
%         % Add title
%         title('每个参数的边缘分布直方图','fontsize',14,'fontweight','bold','fontname','Times');
%     end;
%     % Add x-labels
%     evalstr = strcat('参数',{' '},num2str(counter)); xlabel(evalstr,'fontsize',14,'fontweight','bold','fontname','Times');
%     % Then add y-label (only if j == 1 or j = r;
%     if j == 1 | ( min(abs(j - ([1:r]*t+1))) == 0 ),
%         ylabel('边缘密度','fontsize',14,'fontweight','bold','fontname','Times');
%     end;
%     % Now determine the min and max X values of the plot
%     minX = min(X); maxX = max(X); minY = 0; maxY = max(N/sum(N));
%     % Now determine appropriate scales
%     deltaX = 0.1*(maxX - minX);
%     % Calculate x_min and x_max
%     x_min = minX - deltaX; x_max = maxX + deltaX;
%     % Now determine the min and max Y values of the plot
%     y_min = 0; y_max = 1.1*maxY;
%     % Lets add the MAP value
%     plot(MAP(counter),0.98*y_max,'bx','Markersize',15,'linewidth',3);
%     % Adjust the axis
%     axis([x_min x_max y_min y_max]);
%     % Check if counter = 1,
%     if counter == 1, % --> add a title for first figure
%         % Add title
%         title('每个参数的边缘分布直方图','fontsize',14,'fontweight','bold','fontname','Times');
%     end;
%     % Now update the counter
%     counter = counter + 1;
%     
%     % Update j
%     j = j + 1;
% end;
% 
% % Update fig_number
% fig_number = fig_number + 1;

% ------------------------------------------------------------------------------------------------------------
% -------------------------------------- MARGINAL DENSITIES OF PARAMETERS ------------------------------------
% ------------------------------------------------------------------------------------------------------------

% Plot the histograms (marginal density) of each parameter;
% What lay out of marginal distributions is desired subplot(r,t)
% r = 3; t = 2;
% % How many figures do we need to create with this layout?
% N_fig = ceil( DREAMPar.d / (r * t) ); counter = 1; j = 1;
% % Open new figure
% figure(fig_number);
% % Now plot each parameter
% while counter <= DREAMPar.d
%     % Check whether to open a new figure?
%     if j == (r * t) + 1,
%         % Update fig_number
%         fig_number = fig_number + 1;
%         % Open new figure
%         figure(fig_number);
%         % Reset j to 1
%         j = 1;
%     end;
%     % Now create density
%     [N,X]=density(Pars(:,counter),[]);
%     % And plot density in red
%     subplot(r,t,j),plot(X,N,'r-','linewidth',2); hold on;
%     if j == 1,
%         % Add title
%         title('每个参数的边缘后验密度图','fontsize',14,'fontweight','bold','fontname','Times');
%     end;
%     % Add x-labels
%     evalstr = strcat('参数',{' '},num2str(counter)); xlabel(evalstr,'fontsize',14,'fontweight','bold','fontname','Times');
%     % Then add y-label (only if j == 1 or j = r;
%     if j == 1 | ( min(abs(j - ([1:r]*t+1))) == 0 ),
%         ylabel('边缘密度','fontsize',14,'fontweight','bold','fontname','Times');
%     end;
%     % Now determine the min and max X values of the plot
%     minX = min(X); maxX = max(X); minY = 0; maxY = max(N);
%     % Now determine appropriate scales
%     deltaX = 0.1*(maxX - minX);
%     % Calculate x_min and x_max
%     x_min = minX - deltaX; x_max = maxX + deltaX;
%     % Now determine the min and max Y values of the plot
%     y_min = 0; y_max = 1.1*maxY;
%     % Lets add the MAP value
%     plot(MAP(counter),0.98*y_max,'bx','Markersize',15,'linewidth',3);
%     % Adjust the axis
%     axis([x_min x_max y_min y_max]);
%     % Check if counter = 1,
%     if counter == 1, % --> add a title for first figure
%         % Add title
%         title('每个参数的边缘后验密度图','fontsize',14,'fontweight','bold','fontname','Times');
%     end;
%     % Now update the counter
%     counter = counter + 1;
%     
%     % Update j
%     j = j + 1;
% end;
% 
% % Update fig_number
% fig_number = fig_number + 1;


% ------------------------------------------------------------------------------------------------------------
% -------------------------------- CONVERGENCE OF INDIVIDUAL CHAINS TO TARGET DISTRIUBUTION ------------------
% ------------------------------------------------------------------------------------------------------------

% Define colors for different chains
symbol = {'ys','rx','g+','ko','c<'};

% Now loop over each parameter
for j = 1:DREAMPar.d,
    % Open new figures
    figure(fig_number);
    % Update fig_number
    fig_number = fig_number + 1;
    % How many elements does the chain have
    Nseq = size(chain,1)-1;
    % Now plot a number of chains
    for i = 1:min(DREAMPar.N,5);
        plot([0:Nseq],chain(1:end,j,i),char(symbol(i)),'markersize',3,'linewidth',3); if i == 1; hold on; end;
    end
    % Add an axis
    if isfield(Par_info,'min'),
        % Use scaling with prior parameter ranges
        axis([0 Nseq Par_info.min(j) Par_info.max(j)]);
    else
        % Ranges have not been defined -- need to derive them from ParSet
        min_j = min(ParSet(:,j)); max_j = max(ParSet(:,j));
        % Now make the ranges a little wider
        if min_j < 0,
            min_j = 1.1*min_j;
        else
            min_j = 0.9*min_j;
        end;
        if max_j > 0,
            max_j = 1.1*max_j;
        else
            max_j = 0.9*max_j;
        end;
        % And scale the figure
        axis([0 Nseq min_j max_j]);
    end;
    % Lets add the MAP value
%     plot( Nseq , MAP(j),'bx','Markersize',15,'linewidth',3);
    % Add a legend
    evalstr = strcat('legend(''Markov Chain. 1''');
    % Each parameter a different color
    for jj = 2:min(DREAMPar.N,5),
        % Add the other parameters
        evalstr = strcat(evalstr,',''Markov Chain.',{' '},num2str(jj),'''');
    end;
    % And now conclude with a closing bracket
    evalstr = strcat(evalstr,');');
    % Now evaluate the legend
    eval(char(evalstr));
    % Add a title
    xlabel('No.iterations','fontsize',20,'fontname','Times');
    % Then add a y-label
    ylabel('Stiffness parameters','fontsize',20,'fontname','Times');
    set(gca,'fontsize',20);
    set(gca,'ylim',[0.5 1.5],'ytick',[0.5:0.2:1.5]);
    set(gca,'xlim',[0 15000],'xtick',[0:5000:15000]);
%     evalstr = strcat('Stiffness parameter ',{' '},num2str(j)); ylabel(evalstr,'fontsize',14,'fontweight','bold','fontname','bold','fontname','times');
    % Then add title
%     title('链的收敛图','fontsize',14,'fontweight','bold','fontname','bold','fontname','times');
end;
