% Survival analysis-based dynamic spectrum access algorithm
% Modified spectrum segment allocation method
% Sweeping length of training sequence
%
% Based on 2017 journal and conference paper by T.A. Hall et al.
%--------------------------------------------------------------------------

% Training algorithm with spectrum occupancy data representative of channel
% characteristics

% Simulation variables
Length = 1000000;                    % number of samples in each channel of spectrum occupancy data
channels = 1;                   % number of channels in test occupancy matrix
t = 0;                            % time marker
tau = 1;                          % transmit duration requested
threshold = 0.9;                  % interference threshold (probability of successful transmission)
theta = (-1)*log(threshold);

startM = 1000;
stopM = 20000;
numSweeps = 20;
transTot = zeros(numSweeps, 1);
util = zeros(numSweeps, 1);
interfTot = zeros(numSweeps, 1);
interfRate = zeros(numSweeps, 1);

% Occupancy data
P1 = 10;             % Occupancy event rate (lambda)
P2 = 10;             % Vacancy event rate
%=============================================================================
% Variant 1: Randomly generated occupancy, exponential
%=============================================================================
% trainer = spectrum_occ_exp(1, Length, P2, P1);          % training array for DSA algorithm
% M = spectrum_occ_exp(channels, Length, P2, P1);         % test matrix of occupancy data
%=============================================================================
% Variant 2: Randomly generated occupancy, dual poisson processes
%=============================================================================
M = spectrum_occ_poiss(channels, Length, P1, P2);
trainer = spectrum_occ_poiss(1, Length, P1, P2);
%=============================================================================
% Variant 3: Periodic spectrum occupancy
%=============================================================================
% duty1st = 0.3;
% period1st = 10;
% trainer = [ones(1, period1st * duty1st), zeros(1, period1st - (period1st * duty1st))];
% trainer = repmat(trainer, 1, Length/period1st);
% M = [ones(channels, period1st * duty1st), zeros(channels, period1st - (period1st * duty1st))];
% M = repmat(M, 1, Length/period1st);
%----------------------------------------------------------------------------
occupied = sum(M);
vacant = Length - occupied;
H = zeros(numSweeps, Length);

for m = startM:1000:stopM
    x = m/1000;
    %=====================================================================================================
    % Train DSA algorithm
    %=====================================================================================================
    % Generate array with number of instances of each length of idle period in
    % training vector
    
    counts = occupancy(trainer(1:m));
    transmit = zeros(channels, Length);         % segments where secondary user successfully transmits
    interfere = zeros(channels, Length);        % segments where secondary user collides with primary user
    schedule = zeros(channels, Length + 100);         % transmit grant schedule for secondary user

    n = length(find(counts));

    % Calculate survival/hazard function
    periodsIdle = sum(counts);
    pdf = counts./periodsIdle;
    cdf = cumsum(pdf);
    ccdf = cumsum(pdf, 'reverse');
    %=============================================================================
    % Cumulative Hazard Function
    %=============================================================================
    Ti = [];
    for i = 1:m
        Ti = [Ti, i*ones(1, counts(i))];
    end

    h = zeros(1, m);
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
            H(x, t) = temp;
        elseif t == n   
            H(x, t:Length) = H(x, t-1) + temp;
        else
            H(x, t) = H(x, t-1) + temp;
        end    
    end
    %-----------------------------------------------------------------------------
    
    % Scan test matrix of occupancy data and grant or deny transmission
    % requests
    for i = 1:channels
        t = 0;
        for j = 1:Length
            sample = M(i, j);
            if sample == 0
                t = t + 1;
                if schedule(i, j) == 1
                    transmit(i, j) = 1;
                end
                %=============================================================
                % Algorithm 1
                %=============================================================
                tau = 1;
                T = t + tau;
                if T > Length
                    T = Length; 
                end
                if H(x, T) - H(x, t) < theta
                    schedule(i, (j + 1) : (j + tau)) = 1;
                end                    
                %=============================================================
                % Algorithm 2
                %=============================================================
%                 if ((j > 1) && (M(i, (j-1)) == 1)) 
%                     tau = 1;
%                     while (H(x, t + tau) - H(x, t)) < theta
%                         tau = tau + 1;
%                         if (t + tau) > Ti(periodsIdle)
%                            break
%                         end
%                     end
%                     tau = tau - 1;
%                     schedule(i, (j + 1) : (j + 1 + tau)) = 1;
%                 end
                %------------------------------------------------------------- 
            elseif sample == 1
                t = 0;
                if schedule(i, j) == 1
                    interfere(i, j) = 1;
                end
            end
        end
    end

    % Calculate metrics
    transTot(x) = sum(transmit, 2);
    util(x) = 100*transTot(x)./vacant;
    interfTot(x) = sum(interfere, 2);
    interfRate(x) = 100*interfTot(x)./(Length);

end
