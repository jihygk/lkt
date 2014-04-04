%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Leybourne, Kim and Taylor (2007) full sample
% 'Detecting Multiple Changes in Persistence'
% Studies in Nonlinear Dynamics & Econometrics, 2007, 11(3)
%
% jhk 2013-05-13
% Iterates through LKT-M procedure and finds all I(0) subsamples.
% Requires matlabpool()
%
% USAGE
% -------------
% function [lkt_full_table] = ...
%               lkt_full(y,dates,sig,para) 
%
% INPUTS
% --------
% 
% y:    tx1 vector time-series data
% dates: 3x1 array specifying data frequency/sample period
%       [datetype: 1 = annual, 12 = monthly etc;
%       year of first obs.;
%       month of first obs.]
% sig: {1,2,3} to continue testing at {10,5,1}% significance only
% para: 3x1 array of parameters
%        [T (tau: recommended 0.2);
%        k_max (max number of ADF lags: recommended 4);
%        c (local to unity parameter: recommended -10)]
%
% OUTPUTS
% --------
% lkt_full_table: table listing sequential M-tests conducted until the
%                  end of the sample, finding all significant I(0) periods.
% [ lkt_out lkt_start lkt_end
%                               ...] (see lkt_M.m)
%
% CALLS
% -----
% - lkt_full
%   - lkt_M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Core procedure
% Loop through all I(1) periods

function [full_table,submean_out] = lkt_full(y,dates,sig,para)



% arg checks

if nargin < 4;
    para = [0.2 4 -10];
    if nargin < 3;
        sig = 1;
        if nargin < 2;
           dates = [1 1 1];
        end
    end
else
%     para = varargin{3};
%     sub = varargin{2};
%     dates = varargin{1};
end

full_table = [];
timeperiods = 1:length(y);
times_all = timeperiods;
i0 = []; i1 = [];
datetype = dates(1);
date_year = dates(2);
date_month = dates(3);

k_max = para(2); % k_max

%% Main loop

while numel(timeperiods) > 0;
        % first run
    if length(timeperiods) == length(y);
        t_st = 1;
        t_end = length(y);
        sub = [t_st;t_end]; dates=[datetype;date_year;date_month];

        [lkt_out,lkt_start,lkt_end,~] = lkt_M(y,para,sub,dates);
    else
        %non-first run
        t_st = timeperiods(1);
        for j = 1:length(timeperiods)-1;
            if j == length(timeperiods)-1;
                t_end = timeperiods(end);
                break
            end
            if timeperiods(j+1) ~= timeperiods(j)+1;
                t_end = timeperiods(j);
                break
            end
        end

        sub = [t_st;t_end]; dates=[datetype;date_year;date_month];



        if t_end - t_st + 1 >= 30;

        [lkt_out,lkt_start,lkt_end,~] = lkt_M(y,para,sub,dates);
        else
            lkt_out = [-1 -1 -1];
            lkt_start = [-1 -1 -1];
            lkt_end = [-1 -1 -1];
        end
    end

            if lkt_out(3) >= sig;
                %significant case
                %remove sig part
                i0 = [i0;[lkt_start(1):lkt_end(1)]'];
                 timeperiods = timeperiods(~ismember(timeperiods, ...
                    [lkt_start(1):lkt_end(1)]' ));
                %add results 
                results_temp = [sub(1) sub(2) lkt_out lkt_start lkt_end];
                full_table = [full_table;results_temp];
                
            else
                %remove insig part
                i1 = [i1;[sub(1):sub(2)]'];
                timeperiods = timeperiods(~ismember(timeperiods,[sub(1):sub(2)]'));
                results_temp = [sub(1) sub(2) lkt_out lkt_start lkt_end];
                full_table = [full_table;results_temp];
            end

end


% subsample means

    temp = zeros(length(y),1);
        for n = 1:size(full_table,1)
          if full_table(n,5) >= 1
         submean = mean(y(full_table(n,6):full_table(n,9)));
         temp(full_table(n,6):full_table(n,9)) = submean;
          end
         for o = 1:length(temp)
             if temp(o) == 0
                temp(o) = NaN;
            end
        end
        submean_out = temp;
        end
        
        
end



