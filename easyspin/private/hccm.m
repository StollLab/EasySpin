function C = hccm(J,res,mode)
%
% HCCM  Heteroscedasticity Consistent Covariance Matrix (HCCM)
%
%   C = HCCM(J,V)
%   C = HCCM(J,res,mode)
%
%   Computes the heteroscedasticity consistent covariance matrix (HCCM) of
%   a given LSQ problem given by the Jacobian matrix (J) and the covariance
%   matrix of the data (V). If the residual (res) is specified, the
%   covariance matrix is estimated using some of the methods specified in
%   (mode). The HCCM matrices are valid for both heteroscedastic and
%   homoscedastic residual vectors. 
%
%  Input:
%    J         Jacobian (NxM matrix), with
%                  N = number of data points / variables
%                  M = number of parameters
%    res       vector of residuals (Nx1 array)
%    mode      HCCM estimator (string)
%                 'HC0' - White, H. (1980)
%                 'HC1' - MacKinnon and White, (1985)
%                 'HC2' - MacKinnon and White, (1985)
%                 'HC3' - Davidson and MacKinnon, (1993)
%                 'HC4' - Cribari-Neto, (2004)
%                 'HC5' - Cribari-Neto, (2007)
%
%  Output:
%    C     heteroscedasticity consistent covariance matrix (MxM matrix)
%
% REFERENCES: 
% [1] White, H. (1980). A heteroskedasticity-consistent covariance matrix
%     estimator and a direct test for heteroskedasticity. Econometrica, 48(4), 817-838
%     https://doi.org/10.2307/1912934
%
% [2] MacKinnon and White, (1985). Some heteroskedasticity-consistent covariance
%     matrix estimators with improved finite sample properties. Journal of Econometrics, 
%     29 (1985), pp. 305-325.
%     https://doi.org/10.1016/0304-4076(85)90158-7
%
% [3] Davidson and MacKinnon, (1993). Estimation and Inference in Econometrics
%     Oxford University Press, New York. 
%
% [4] Cribari-Neto, F. (2004). Asymptotic inference under heteroskedasticity of
%     unknown form. Computational Statistics & Data Analysis, 45(1), 215-233
%     https://doi.org/10.1016/s0167-9473(02)00366-3
%
% [5] Cribari-Neto, F., Souza, T. C., & Vasconcellos, K. L. P. (2007). Inference
%     under heteroskedasticity and leveraged data. Communications in Statistics –
%     Theory and Methods, 36(10), 1877-1888.
%     https://doi.org/10.1080/03610920601126589

% Number of variables
N = size(J,1);
% Number of parameters
M = size(J,2);

if nargin<3 || isempty(mode) || all(size(res)>1)
    diagV = res;
else
    diagV = [];
end

pinvJTJ = pinv(J.'*J);
% Get leverage (diagonal of hat matrix)
%H = J*pinvJTJ*J.'; h = diag(H);  % consumes too much memory if 
h = sum((J*pinvJTJ).*J,2);

if isempty(diagV)
    % Select estimation method using established nomenclature
    % to estimate the diagonal of the data covariance matrix
    switch upper(mode)
        
        case 'HC0' % White,(1980),[1]
            diagV = res.^2;
            
        case 'HC1' % MacKinnon and White,(1985),[2]
            diagV = N/(N-M)*res.^2;
            
        case 'HC2' % MacKinnon and White,(1985),[2]
            % Estimate the diagonal of the data covariance matrix
            diagV = res.^2./(1-h);
            
        case 'HC3' % Davidson and MacKinnon,(1993),[3]
            diagV = (res./(1-h)).^2;
            
        case 'HC4' % Cribari-Neto,(2004),[4]
            delta = min(4,N*h/M);  % discount factor
            diagV = res.^2./((1 - h).^delta);
            
        case 'HC5' % Cribari-Neto,(2007),[5]
            M = 0.7;
            alpha = min(max(4,M*max(h)/mean(h)),h/mean(h));  % inflation factor
            diagV = res.^2./(sqrt((1 - h).^alpha));
            
        otherwise
            error('HCCM estimation mode not found.')
    end
end

% Heteroscedasticity Consistent Covariance Matrix (HCCM) estimator
C = pinvJTJ*J.'.*diagV(:).'*J*pinvJTJ;

end
