function r = alpha_noise(alpha,beta,gamma,delta,varargin)
%输入形式：r=alpha_noise(alpha,beta,gamma,delta,M,N,......)
%输出为 M*N*......维的r噪声
%例如：r=alpha_noise(alpha,beta,gamma,delta,1,1000）
%输出为1*1000长度的噪声
%
% ALPHA,BETA,GAMMA and DELTA must be scalars which fall in the following 
% ranges :
%    0 < ALPHA <= 2 特征指数，用来表示概率密度函数拖尾的厚度，当ALPHA逼近0，代表脉冲尖峰出现的越频繁；当ALPHA=2，退化成Gaussian分布
%    -1 <= BETA <=1 对称系数，代表关于DELTA的倾斜程度，当BETA=0，对称，反之不对称。
%    0 < GAMMA < inf 分散系数，越大代表脉冲尖峰的强度越大
%    -inf < DELTA < inf  位置参数
%
% References:
% [1] J.M. Chambers, C.L. Mallows and B.W. Stuck (1976) 
%     "A Method for Simulating Stable Random Variables"  
%     JASA, Vol. 71, No. 354. pages 340-344  
%
% [2] Aleksander Weron and Rafal Weron (1995)
%     "Computer Simulation of Levy alpha-Stable Variables and Processes" 
%     Lec. Notes in Physics, 457, pages 379-392
%

if nargin < 4
    error('stats:stblrnd:TooFewInputs','Requires at least four input arguments.'); 
end

%varargin=[1,length(x)];

% Check parameters
if alpha <= 0 || alpha > 2 || ~isscalar(alpha)
    error('stats:stblrnd:BadInputs',' "alpha" must be a scalar which lies in the interval (0,2]');
end
if abs(beta) > 1 || ~isscalar(beta)
    error('stats:stblrnd:BadInputs',' "beta" must be a scalar which lies in the interval [-1,1]');
end
if gamma < 0 || ~isscalar(gamma)
    error('stats:stblrnd:BadInputs',' "gamma" must be a non-negative scalar');
end
if ~isscalar(delta)
    error('stats:stblrnd:BadInputs',' "delta" must be a scalar');
end


% Get output size
[err, sizeOut] = genOutsize(4,alpha,beta,gamma,delta,varargin{:});
if err > 0
    error('stats:stblrnd:InputSizeMismatch','Size information is inconsistent.');
end


%---Generate sample----

% See if parameters reduce to a special case, if so be quick, if not 
% perform general algorithm

if alpha == 2                  % Gaussian distribution 
    r = sqrt(2) * randn(sizeOut);

elseif alpha==1 && beta == 0   % Cauchy distribution
    r = tan( pi/2 * (2*rand(sizeOut) - 1) ); 

elseif alpha == .5 && abs(beta) == 1 % Levy distribution (a.k.a. Pearson V)
    r = beta ./ randn(sizeOut).^2;

elseif beta == 0                % Symmetric alpha-stable
    V = pi/2 * (2*rand(sizeOut) - 1); 
    W = -log(rand(sizeOut));          
    r = sin(alpha * V) ./ ( cos(V).^(1/alpha) ) .* ...
        ( cos( V.*(1-alpha) ) ./ W ).^( (1-alpha)/alpha ); 

elseif alpha ~= 1                % General case, alpha not 1
    V = pi/2 * (2*rand(sizeOut) - 1); 
    W = - log( rand(sizeOut) );       
    const = beta * tan(pi*alpha/2);
    B = atan( const );
    S = (1 + const * const).^(1/(2*alpha));
    r = S * sin( alpha*V + B ) ./ ( cos(V) ).^(1/alpha) .* ...
       ( cos( (1-alpha) * V - B ) ./ W ).^((1-alpha)/alpha);

else                             % General case, alpha = 1
    V = pi/2 * (2*rand(sizeOut) - 1); 
    W = - log( rand(sizeOut) );          
    piover2 = pi/2;
    sclshftV =  piover2 + beta * V ; 
    r = 1/piover2 * ( sclshftV .* tan(V) - ...
        beta * log( (piover2 * W .* cos(V) ) ./ sclshftV ) );      
          
end
    
% Scale and shift
if alpha ~= 1
   r = gamma * r + delta;
else
   r = gamma * r + (2/pi) * beta * gamma * log(gamma) + delta;  
end

end


%====  function to find output size ======%
function [err, commonSize, numElements] = genOutsize(nparams,varargin)
%varargin=[1,length(x)];
    try
        tmp = 0;
        for argnum = 1:nparams
            tmp = tmp + varargin{argnum};
        end
        if nargin > nparams+1
            tmp = tmp + zeros(varargin{nparams+1:end});
        end
        err = 0;
        commonSize = size(tmp);
        numElements = numel(tmp);

    catch
        err = 1;
        commonSize = [];
        numElements = 0;
    end
end