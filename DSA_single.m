% Channel parameter estimation and time slot allocation scheme
% Based on 2012 paper by Q. Liang et al. and 2014 paper by Sharma and Sahoo

% Simulation variables
%--------------------------------------------------------------------------
Length = 1000000;                    % number of samples in each channel of spectrum occupancy data
channels = 1;
L0 = 14;                         % Vacancy event rate
L1 = 10;                         % Occupancy event rate (lambda)
threshold = 0.9;                  % interference threshold (probability of successful transmission)
theta = (-1)*log(threshold);
eta = 1 - threshold;
tau = 1;
Tcs = 10000;                        % Continuous sensing window size
Trs = 90000;                       % Random sensing window size
Tw = Tcs + Trs;                 % Total window size
nWin = floor(Length/Tw);        % number of sampling windows
m = Trs / 5;                         % starting number of samples per window

% Generate spectrum occupancy array
%--------------------------------------------------------------------------
A = spectrum_occ_exp(channels, Length, L0, L1);
% A = spectrum_occ_poiss(channels, Length, L1, L0);
%--------------------------------------------------------------------------

transTot = zeros(nWin, 1);
util = zeros(nWin, 1);
interfTot = zeros(nWin, 1);
interfRate = zeros(nWin, 1);
nSamples = zeros(nWin, 1);
efficiencyRS = zeros(nWin, 1);
efficiency = zeros(nWin, 1);

occupied = zeros(nWin, 1);       % number of occupied samples per continuous sensing window
vacant = zeros(nWin, 1);         % number of vacant samples per continuous sensing window
occupiedTot = zeros(nWin, 1);
vacantTot = zeros(nWin, 1);
u = zeros(nWin, 1);             % estimation of channel utilization
H = zeros(nWin, Tcs);           % cumulative hazard function


for k = 1:nWin
    % Continuous sensing for first window
    startCS = (k-1)*Tw + 1;
    stopCS = (k-1)*Tw + Tcs;
        
    transmit = zeros(1, Trs);         % segments where secondary user successfully transmits
    interfere = zeros(1, Trs);        % segments where secondary user collides with primary user
    schedule = zeros(1, Trs + 100);         % transmit grant schedule for secondary user

    countsIdle = zeros(1, Tcs);
    countsOcc = zeros(1, Tcs);
    t1 = 0;
    t2 = 0;
    z = zeros(1, Tcs);
    for i = 1:Tcs
        T = (k-1)*Tw + i;
        z(i) = A(T);
        if z(i) == 0
            t2 = 0;
            t1 = t1 + 1;
            vacant(k) = vacant(k) + 1;
            if ((i+1) > Tcs) || (A(T+1) == 1)
               countsIdle(t1) = countsIdle(t1) + 1; 
            end
        elseif z(i) == 1
            t1 = 0;
            t2 = t2 + 1;
            occupied(k) = occupied(k) + 1;
            if ((i+1) > Tcs) || (A(T+1) == 0)
               countsOcc(t2) = countsOcc(t2) + 1; 
            end
        end
    end
    
    occupied(k) = sum(z);
    vacant(k) = Tcs - occupied(k);
    u(k) = occupied(k)/Tcs;
    n = length(find(countsIdle));
    periodsIdle = sum(countsIdle);
    pdf = countsIdle./periodsIdle;
    cdf = cumsum(pdf);
    
    meanI = 0;
    for i = 1:Tcs
        meanI = meanI + (i * pdf(i));
    end
    fri = (1 - cdf)./meanI;
    FRI = cumsum(fri);

    
    for i = 1:Tcs
        if FRI(i) <= eta
            ymax = i; 
        else
            break
        end
    end

    %=============================================================================
    % Cumulative Hazard Function
    %=============================================================================
    Ti = [];
    for i = 1:Tcs
        Ti = [Ti, i*ones(1, countsIdle(i))];
    end

    h = zeros(1, Tcs);
    j = 1;
    for t = 1:n
        temp = 0;    
        while t >= Ti(j)
            temp = temp + 1/(periodsIdle - j + 1);
            j = j + 1;
            if j > periodsIdle
                j = periodsIdle;
                break
            end
        end
        h(t) = temp;
        if t == 1
            H(k, t) = temp;
        elseif t == n   
            H(k, t:Tcs) = H(k, t-1) + temp;
        else
            H(k, t) = H(k, t-1) + temp;
        end    
    end
    %-----------------------------------------------------------------------------

    % Generate bit mask for random sampling period
    sampler = [generate_intervals2(Trs, m), 0] ;
    
    % Apply random sensing to rest of window
    for i = 1:channels
        t = 0;
        for j = 1:Trs
            T = (k-1)*Tw + Tcs + j;
            if sampler(j) == 1
                current = A(T);
                if current == 1
                    sampler(j+1) = 1;
                elseif current == 0
                    if ((j<2)&&(z(Tcs)==1))
                        %--------------------------------------------------
                        % Hall Algorithm 1
                        %--------------------------------------------------
%                         T0 = 1 + tau;
%                         if T0 > Trs
%                             T0 = Trs; 
%                         end
%                         if H(k, T0) - H(k, 1) < theta
%                             schedule(i, (j + 1) : (j + tau)) = 1;
%                         end
                        %--------------------------------------------------
                        % Hall Algorithm 2
                        %--------------------------------------------------
                        tau = 1;
                        while (H(k, 1 + tau) - H(k, 1)) < theta
                            tau = tau + 1;
                            if (1 + tau) > Ti(periodsIdle)
                               break
                            end
                        end
                        tau = tau - 1;
                        schedule(i, (j + 1) : (j + 1 + tau)) = 1;
                        %--------------------------------------------------
                    elseif j > 1
                        if (sampler(j-1) == 1) && (A(T-1) == 1)
                            %----------------------------------------------
                            % Hall Algorithm 1
                            %----------------------------------------------
%                             T0 = 1 + tau;
%                             if T0 > Trs
%                                 T0 = Trs; 
%                             end
%                             if H(k, T0) - H(k, 1) < theta
%                                 schedule(i, (j + 1) : (j + tau)) = 1;
%                             end
                            %----------------------------------------------
                            % Hall Algorithm 2
                            %----------------------------------------------
                            tau = 1;
                            while (H(k, 1 + tau) - H(k, 1)) < theta
                                tau = tau + 1;
                                if (1 + tau) > Ti(periodsIdle)
                                   break
                                end
                            end
                            tau = tau - 1;
                            schedule(i, (j + 1) : (j + 1 + tau)) = 1;
                            %----------------------------------------------
                        end
                    else
                        %----------------------------------------------------------
                        % RIBS
                        %----------------------------------------------------------
                        schedule(i, j:j+ymax) = 1;   
                    end
                end               
            end
            if schedule(i, j) == 1
                if A(i, j) == 0
                    transmit(i, j) = 1; 
                elseif A(i, j) == 1
                    interfere(i, j) = 1;
                end
            end
        end
    end
    
    % Calculate metrics
    occupiedTot(k) = sum(A(:, ((k-1)*Tw + 1):(k*Tw)));
    vacantTot(k) = Tw - occupiedTot(k);
    
    temp1 = sum(transmit, 2);
    temp2 = temp1./vacantTot(k);
    transTot(k) = mean(temp1, 1);
    util(k) = 100 * mean(temp2, 1);
    
    temp3 = sum(interfere, 2);
    temp4 = temp3./Tw;
    interfTot(k) = mean(temp3, 1);
    interfRate(k) = 100 * mean(temp4, 1);
    
    nSamples(k) = sum(sampler);
    efficiencyRS(k) = 100 * (Trs-nSamples(k))./Trs;
    efficiency(k) = 100 * (Tw - (nSamples(k) + Tcs))./Tw;
    
    %==========================================================================
    % Maximum Likelihood Estimation
    %==========================================================================
    % L = 0;
    % S1 = u^(z(1))*(1 - u)^(1-z(1));
    % for i = 2:m2(1)
    %     S2 = u^z(i)*(1 - u)^(1-z(i));
    %     S3 = (-1)^(z(i)+z(i-1)).*u^(1-z(i-1)).*(1-u)^z(i-1);
    %     S4 = exp(-1*(thetVec.*delT(i))/u);
    %     L = L + log((S2 + S3.*S4));
    % end
    % lnL = L + log(S1);
    % [Lmax, indMax] = max(lnL);
    % theta0(1) = thetVec(indMax);
    % theta1(1) = (1-u)*theta0(1)/u;
    % 
    % pdf0 = exppdf(1:Tw, 1/theta0(1));
    % F0 = cumsum(pdf0);
    % pdf1 = exppdf(1:Tw, 1/theta1(1));
    % F1 = cumsum(pdf1);
    % fRI = (1 -F1)*theta0(1);
    % FRI = cumsum(fRI);
    %--------------------------------------------------------------------------

    %==========================================================================
    % RIBS stuff
    %==========================================================================
    % n0 = find(countsIdle, 1, 'last');
    % n1 = find(countsOcc, 1, 'last');
    % periodsIdle = sum(countsIdle);
    % periodsOcc = sum(countsOcc);
    % pdf0 = countsIdle./periodsIdle;
    % pdf1 = countsOcc./periodsOcc;
    % FI = cumsum(pdf0);     % Cumulative Distribution function for channel idle time
    % FO = cumsum(pdf1);
    % 
    % meanI = 0;
    % for i = 1:n0
    %     meanI = meanI + (i * pdf0(i));
    % end
    % 
    % fRI = (1 -FI)/meanI;
    % FRI = cumsum(fRI);
    % for i = 1:n0
    %    if FRI(i) <= eta
    %       ymax = i; 
    %    else
    %        break
    %    end
    % end
    %--------------------------------------------------------------------------
end


