% Survival analysis-based dynamic spectrum access algorithm
% Modified spectrum segment allocation method
% Channel mean idle time varies over length of sampling period
%
% Based on 2017 journal and conference paper by T.A. Hall et al.
%--------------------------------------------------------------------------

% Training algorithm with spectrum occupancy data representative of channel
% characteristics

% Simulation variables
Length = 1000000;                    % number of samples in each channel of spectrum occupancy data
channels = 1;                   % number of channels in test occupancy matrix
Tw = 5000;                      % window length
nWin = floor(Length/Tw);            % number of windows in sample array
nTrain = 10;                        % number of windows before retraining
t = 0;                            % time marker
tau = 5;                          % transmit duration requested
threshold = 0.95;                  % interference threshold (probability of successful transmission)
theta = (-1)*log(threshold);

transTot = zeros(nWin, 1);
util = zeros(nWin, 1);
interfTot = zeros(nWin, 1);
interfRate = zeros(nWin, 1);

% Occupancy data
%--------------------------------------------------------------------------
m = 1:nWin;
% P1 = 50;                                % Occupancy event rate (lambda)

%--------------------------------------------------------------------------
% Sinusoidal channel parameter variations
P1 = 20 - 15*sin(2*pi*m/100);             % Occupancy event rate
P2 = 20 + 15*sin(2*pi*m/100);            % Vacancy event rate
%--------------------------------------------------------------------------
% % Linear channel parameter variations
% P2 = 2 + abs(100-m);
% P1(1:100) = 2 + m(1:100);
% P1(101:200) = 102 - abs(100 - m(101:200));
%=============================================================================
% Variant 1: Time-varying mean vacancy, exponential
%=============================================================================
% M = [];
% for i = 1:nWin
%    M =  [M, spectrum_occ_exp(channels, Tw, P2(i), P1)];
% end
%=============================================================================
% Variant 2: Time-varying mean vacancy, dual poisson processes
%=============================================================================
M = [];
for i = 1:nWin
   M =  [M, spectrum_occ_poiss(channels, Tw, P1(i), P2(i))];
end
%=============================================================================

for m = 1:nWin
    start = (m-1)*Tw + 1;
    stop = m*Tw;
    window = M( : , start:stop);
    transmit = zeros(channels, Tw);         % segments where secondary user successfully transmits
    interfere = zeros(channels, Tw);        % segments where secondary user collides with primary user
    schedule = zeros(channels, Tw + 100);         % transmit grant schedule for secondary user
    
    occupied = sum(window);
    vacant = Tw - occupied;
    
    % Initial training on first window
    if m == 1
        %=====================================================================================================
        % Train DSA algorithm
        %=====================================================================================================
        % Generate array with number of instances of each length of idle period in
        % training vector
        counts = occupancy(window);
        n = length(find(counts));
        H = zeros(1, Tw);

        % Calculate survival/hazard function
        periodsIdle = sum(counts);
        pdf = counts./periodsIdle;
        cdf = cumsum(pdf);
        ccdf = cumsum(pdf, 'reverse');
        %=============================================================================
        % Initial Cumulative Hazard Function
        %=============================================================================
        Ti = [];
        for i = 1:Tw
            Ti = [Ti, i*ones(1, counts(i))];
        end

        h = zeros(1, Tw);
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
                H(t) = temp;
            elseif t == n   
                H(t:Tw) = H(t-1) + temp;
            else
                H(t) = H(t-1) + temp;
            end    
        end
        %-----------------------------------------------------------------------------
    else                
        %==================================================================
        % Retrain every nTrain windows
        %==================================================================
%         if rem(m, nTrain) == 1
%             % Calculate survival/hazard function
%             periodsIdle = sum(counts);
%             n = length(find(counts));
%               
%             % Periodic Rebuild of Cumulative Hazard Function
%             %--------------------------------------------------------------
%             Ti = [];
%             for i = 1:Tw
%                 Ti = [Ti, i*ones(1, counts(i))];
%             end
% 
%             h = zeros(1, Tw);
%             H = zeros(1, Tw);
%             j = 1;
%             for t = 1:n
%                 temp = 0;    
%                 while t >= Ti(j)
%                     temp = temp + 1/(periodsIdle - j + 1);
%                     j = j + 1;
%                     if j > periodsIdle
%                         j = periodsIdle;
%                         break
%                     end
%                 end
%                 h(t) = temp;
%                 if t == 1
%                     H(t) = temp;
%                 elseif t == n   
%                     H(t:Tw) = H(t-1) + temp;
%                 else
%                     H(t) = H(t-1) + temp;
%                 end    
%             end
%         end
        %==================================================================
        % Retrain every window
        %==================================================================
        % Calculate survival/hazard function
        periodsIdle = sum(counts);
        n = length(find(counts));
        H = zeros(1, Tw);

        % Rebuild Cumulative Hazard Function
        %------------------------------------------------------------------
        Ti = [];
        for i = 1:Tw
            Ti = [Ti, i*ones(1, counts(i))];
        end
        h = zeros(1, Tw);
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
                H(t) = temp;
            elseif t == n   
                H(t:Tw) = H(t-1) + temp;
            else
                H(t) = H(t-1) + temp;
            end    
        end
        %-----------------------------------------------------------------------------
              
        % Scan test matrix of occupancy data and grant or deny transmission
        % requests
        counts = zeros(channels, Tw);
        for i = 1:channels
            t = 0;
            for j = 1:Tw
                sample = window(i, j);
                if sample == 0
                    t = t + 1;
                    if ((j+1) > Tw) || (window(i, (j+1)) == 1)
                        counts(i, t) = counts(i, t) + 1;
                    end
                    if schedule(i, j) == 1
                        transmit(i, j) = 1;
                    end
                    %=============================================================
                    % Algorithm 1
                    %=============================================================
                    T = t + tau;
                    if T > Length
                        T = Length; 
                    end
                    if H(T) - H(t) < theta
                        schedule(i, (j + 1) : (j + tau)) = 1;
                    end                    
                    %=============================================================
                    % Algorithm 2
                    %=============================================================
%                     if ((j > 1) && (M(i, (j-1)) == 1)) 
%                         tau = 1;
%                         while (H(t + tau) - H(t)) < theta
%                             tau = tau + 1;
%                             if (t + tau) > Ti(periodsIdle)
%                                break
%                             end
%                         end
%                         tau = tau - 1;
%                         schedule(i, (j + 1) : (j + 1 + tau)) = 1;
%                     end
                    %------------------------------------------------------------- 
                elseif sample == 1
                    t = 0;
                    if schedule(i, j) == 1
                        interfere(i, j) = 1;
                    end
                end
            end
        end
    end

    % Calculate metrics
    transTot(m) = sum(transmit, 2);
    util(m) = 100*transTot(m)./vacant;
    interfTot(m) = sum(interfere, 2);
    interfRate(m) = 100*interfTot(m)./(Tw);

end

UTIL = mean(util);
INTERF = mean(interfRate);
