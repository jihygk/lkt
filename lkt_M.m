%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Leybourne, Kim and Taylor (2007) M-test
% 'Detecting Multiple Changes in Persistence'
% Studies in Nonlinear Dynamics & Econometrics, 2007, 11(3)
%
% jhk 2013-05-13
% Runs one iteration of the LKT-M procedure. Finds one I(0) subsample.
% Requires matlabpool()
%
% USAGE
% -------------
% 1) lkt_out = lkt_M(y)
% 2) [lkt_out,lkt_start,lkt_end,lkt_table] = lkt_M(y,[para],[sub],[date])
%
% INPUTS
% ------
% y: nx1 time-series to perform LKT-M test
%
% OPTIONAL INPUTS
% -----------------
% para: 3x1 array of parameters
%        [T (tau: recommended 0.2);
%        k_max (max number of ADF lags: recommended 4);
%        c (local to unity parameter: recommended -10)]
% sub: 2x1 array of subsample specification
%        [subsample start index; 
%        subsample end index]
% date: 3x1 array of date data
%        [date format (12 = monthly, 1 = annual, etc.);
%        date year (corresponding to y(1));
%        date month (corresponding to y(1))]
%
% OUTPUTS
% ------
% lkt_out = [M-value, k (ADF lags), significance*];
% lkt_start = [index, year, month] of beginning of detected I(0) period
% lkt_end = [index, year, month] of end of detected I(0) period
% * 1 = 10%, 2 = 5%, 3 = 1%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [lkt_out,varargout] = lkt_M(y,varargin)

% arg checks

if nargin < 4;
    date = [1 1 1];  
    if nargin < 3;
        sub = [1 length(y)];
        if nargin < 2;
        para = [0.2 4 -10];
        end
    end
else
    para = varargin{1};
    sub = varargin{2};
    date = varargin{3};
end



% Hardcoded critical values at 10/5/1% (LKT, 2007)
% for tau = 0.2, c = -10
cri = [20 -4.736 -5.369 -7.53;
30 -4.391 -5.024 -6.618;
40 -4.082 -4.515 -5.662;
50 -3.954 -4.351 -5.292;
60 -3.883 -4.24 -5.133;
70 -3.803 -4.143 -4.974;
80 -3.744 -4.083 -4.811;
90 -3.735 -4.088 -4.781;
100 -3.718 -4.049 -4.734;
120 -3.699 -4.026 -4.669;
140 -3.669 -3.992 -4.614;
160 -3.668 -3.963 -4.558;
180 -3.657 -3.973 -4.512;
200 -3.662 -3.964 -4.536;
250 -3.645 -3.917 -4.513;
300 -3.646 -3.926 -4.466;
350 -3.642 -3.925 -4.463;
400 -3.627 -3.9 -4.438;
10^5 -3.616 -3.885 -4.421];


sub_s = sub(1);
sub_e = sub(2);
datetype = date(1);
date_y = date(2);
date_m = date(3);
tau = para(1);
k_max = para(2);
c = para(3);
y_orig = y;
y = y(sub_s:sub_e);
alpha = 1 + c/size(y,1);

%% Main looping for DF_G 

% default values
% lkt_M = 9999;
% %lkt_temp = 9999;
% lkt_end = -1;
% lkt_start = -1;
% lkt_k = -1;
% lkt_sig = -1;
% max_c = 0;


% generate index for parfor
indices = zeros(length(y)*length(y),2);
kk = 0;
for ii = 1:size(y,1);
    for jj = 1:size(y,1);
        kk = kk+1;
        indices(kk,1) = ii;
        indices(kk,2) = jj;
    end
end
ll = indices(:,1) + ceil(tau*size(y,1)) <= indices(:,2) ;
indices = indices(ll,:);
clear ii jj kk ll
parend = length(indices);
lkt_results = zeros(length(indices),4);

index_lambda = indices(:,1);
index_tau = indices(:,2);


% main parfor loop
parfor parind = 1:parend; 
    [lkt_temp, k_star] = lkt_dfg4(y,index_lambda(parind),index_tau(parind),alpha,k_max);
    lkt_results(parind,:) = [lkt_temp k_star index_lambda(parind) index_tau(parind)];
end



[~,lkt_inf] = sortrows(lkt_results,1);
[lkt_M] = lkt_results(lkt_inf(1),1);
[lkt_k] = lkt_results(lkt_inf(1),2);
[lkt_start] = lkt_results(lkt_inf(1),3);
[lkt_end] = lkt_results(lkt_inf(1),4);



%% Critical value
crit = interp1(cri(:,1),cri(:,2:4),length(y));
lkt_sig = sum(lkt_M < crit);


%% Dates and return values
%[lkt_sy,lkt_sm,lkt_ey,lkt_em] = lkt_date(datetype,date_y,date_m,lkt_start,lkt_end,sub_s);
%function [d1,d2]=lkt_date2(T,frac,start_y,start_m,index_s,index_e,sub_st)
[aa,bb] = lkt_date2(size(y_orig,1),datetype,date_y,date_m,lkt_start,lkt_end,sub_s);
lkt_sy = aa(1); lkt_sm = aa(2); lkt_ey = bb(1); lkt_em = bb(2); % dates
lkt_start = lkt_start - 1 + sub_s;
lkt_end = lkt_end - 1 + sub_s;

% function outputs
lkt_out = [lkt_M,lkt_k,lkt_sig];
varargout{1} = [lkt_start lkt_sy lkt_sm];
varargout{2} = [lkt_end lkt_ey lkt_em];
varargout{3} = lkt_results;

end

%% lkt_dfg

function [dfg,k_star] = lkt_dfg4(y,y_start,y_end,alpha,k_max)
% GLS fitted
y_dfg = y(y_start:y_end);

z_LT = [1;ones(y_end-y_start,1).*(1-alpha)];
y_LT = [y_dfg-[0;y_dfg(1:end-1)].*alpha];

% beta
b = ((z_LT'*z_LT)^-1)*z_LT'*y_LT;
% fitted
yd = y_dfg - ones(y_end-y_start+1,1)*b;
yd_1 = [NaN;yd(1:end-1)];

% ADF test

% LHS of (6.4)
ADF_lhs = yd - yd_1;

% RHS of (6.4): make Deltas of lag
ADF_rhs2 = zeros(size(y_dfg,1),k_max);
ADF_rhs2(:,1) = yd_1;
for c = 2:k_max
    try
        ADF_rhs2(:,c) = [ones(c-1,1).*NaN;y_dfg(1:end-c+1)] - [ones(c,1).*NaN;y_dfg(1:end-c)];
    catch dimension;
        ADF_rhs2(:,c) = NaN;        
    end
end
kindex = zeros(k_max,2);

% %lag selection: AICc
% for j = k_max:-1:1
%     [ADF_beta, ADF_se, ADF_r] = lkt_ols(ADF_lhs(j+1:end),ADF_rhs2(j+1:end,1:j));
%     
%     rr = ADF_r'*ADF_r;
%     yyr = size(ADF_lhs(j+1:end),1);
%     ll = -yyr/2*(log(rr/yyr));
%     AICc = (-2*ll/yyr)+(2*j/yyr)+(2*j*(j+1)/(yyr-j-1))/yyr;
%     kindex(j,1) = AICc;
%     kindex(j,2) = ADF_beta(1)/ADF_se(1);
%         
% end
% [~, k_star] = min(kindex(:,1));
% dfg = kindex(k_star,2);

% %lagselection: BIC
% for j = k_max:-1:1
%     [ADF_beta, ADF_se, ADF_r] = lkt_ols(ADF_lhs(j+1:end),ADF_rhs2(j+1:end,1:j));
%     
%     rr = ADF_r'*ADF_r;
%     yyr = size(ADF_lhs(j+1:end),1);
%     ll = -yyr/2*(log(rr/yyr));
%     BIC = -2.*ll + j*log(yyr);
%     kindex(j,1) = BIC;
%     kindex(j,2) = ADF_beta(1)/ADF_se(1);
% end
% [~, k_star] = min(kindex(:,1));
% dfg = kindex(k_star,2);

%lag selection:  Ng-Perron (1995)
for j = k_max:-1:1
    
    [ADF_beta,ADF_se,~] = lkt_ols(ADF_lhs(j+1:end),ADF_rhs2(j+1:end,1:j));
    
    if j == 1;
        k_star = 1;
        dfg = ADF_beta./ADF_se;
        break
    end
    
    if ADF_beta(j)/ADF_se(j) > 1.65;
        k_star = j;
        dfg = ADF_beta(1)/ADF_se(1);
        break
    end
    
end


end

%% lkt_ols
function [betas, stderrors, resids, yhat, varcov] = lkt_ols(y,x)

betas = (x'*x) \ (x' * y);
yhat = x*betas;     
resids = y - yhat;    
[T, k] = size(y);
residvar = resids'*resids/(T-k);
varcov = residvar*((x'*x)^-1);
stderrors = sqrt(diag(varcov));

end

%% lkt_date

function [d1,d2]=lkt_date2(T,frac,start_y,start_m,index_s,index_e,sub_st)
if start_m > frac; error('Invalid starting month: %d, with fraction: %d.',start_m,frac); end
y = floor(start_y+(start_m-1)/frac+(0:frac*T-1)'./frac);
f = mod((0:frac*T-1)'+start_m,frac); f2 = f == 0; f(f2) = frac;
date = [y f]; date = date(sub_st:end,:);
d1 = date(index_s,:); d2 = date(index_e,:);
end


