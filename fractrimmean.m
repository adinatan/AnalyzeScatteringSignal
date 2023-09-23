
function [m, s] = fractrimmean(x,percent,dim)
%FRACTRIMMEAN Fractional trimmed mean.
%   Calculates the M fractional trimmed mean (truncated mean) of the values 
%   in X. Excluding the highest and lowest PERCENT/2 (%) of the data. It is
%   an advanced trimming method. Given a sample of n ordered observations 
%   x_(1) <= x_(2) <= ... <= x_(n), the trimmed mean is a specific example 
%   of weighted mean or L-estimator of location. It will apply zero weights
%   in the extremes of the sample and positive weights in the center. This 
%   weighting scheme is reasonable when the underlying distribution is 
%   symmetric but the sample from it may be contaminated by outliers who 
%   can bias the averaging. So, it is robust to the effects of outliers and
%   with small variance.
%
%   Not to be wrong. It is necessary to be careful in distinguishing 
%   between the PERCENT of both sides of the sample distribution and the 
%   TRIMMED PERCENT, which is PERCENT/2.   
%
%   Thus, a trimmed mean is calculated by discarding a certain percentage
%   of the lowest and the highest scores and then computing the mean of the
%   remaining scores. The median is the mean trimmed 100% and the
%   arithmetic mean is the mean trimmed 0%. Obviously, it is less 
%   susceptible to the effects of extreme scores than is the arithmetic 
%   mean.
%  
%   As Horn and Kafadar (2002) stated. A trimmed mean requires a 
%   prespecified proportion of deleted observations, say p, often stated as
%   a percentage. The p-trimmed mean (trimmed percent) is defined as,
%
%     T(p) = ((1 - r)(x_(g+1) + x_(n-g)) + i=g+2_S_n-g-1 x_(i))/n(1 - 2p)
%     
%   where g=[np] and r=np - g, the integer and fractional remainder parts
%   of the product of n and p, respectively. S is the sum operator. If p*n 
%   is an integer, then r=0 and there are no partially weighted extreme 
%   values. The denominator part n(1 - 2p) is the retained number of 
%   observations after trimming p(%) on each side. Then,
%
%                  T(p) = (i=g+1_S_n-g x_(i))/(n(1 - 2p))
%
%   Thus, as in Lim and Loh (1996, p.292), the 20% fractional trimmed mean 
%   (40%) at n=4 is, with g=[4*0.2]=[0.8]=0.0 and r=0.8-0.0=0.8:
%
%   T(0.2) = ((1 - 0.8)(x_(0.0+1) + x_(4-0.0)) + i=0.0+2_S_4-0.0-1 x_(i))/4(1-2*0.2) 
%          = ((0.2)(x_(1) + x_(4)) + i=2_S_3 x_(i))/4(0.6)
%          = ((0.2)(x_(1) + x_(4)) + (x_(2) + x_(3)))/2.4
%          = (0.2*x_(1) + x_(2) + x_(3) + 0.2*x_(4))/2.4
%
%   at n=5 is, with g=[5*0.2]=[1.0]=1.0 and r=1.0-1.0=0.0:
%
%               T(0.2) = (i=1+1_S_5-1 x_(i))/(5(1 - 2*0.2))
%                      = (i=2_S_4 x_(i))/(5(1-0.4))
%                      = (i=2_S_4 x_(i))/5(0.6)
%                      = (i=2_S_4 x_(i))/3
%                      = (x_(2) + x_(3) + x_(4))/3
%
%   For a matrix input, if required M can be a row vector containing the 
%   trimmed mean of each column of X (DIM = 1,default) or a column vector of 
%   the trimmed mean of each row of X (DIM = 2). PERCENT is a scalar between 
%   0 and 100.
%
%   FRACTRIMMEAN treats NaNs values as missing values, and removes them.
%
%   As a note. Some computer packages' trimmed mean function does not do 
%   advanced trimming but only simple trimming. This is, discarding a
%   certain bottom and top end order statistics considered as 'percentage' 
%   and then computing the mean of the remaining central ones. As you can
%   see in the follow example. If you have the x_(i) ordered statistics,
%
%                       1,2,2,4,5,5,6,6,7,7,8,10,13
%
%   If you call the corresponding trimmed mean M for a given percentage P:
%
%     -----------------------------------------------------
%      P           M                   mean of
%     -----------------------------------------------------
%     0-7        5.8462       [1,2,2,4,5,5,6,6,7,7,8,10,13]
%     8-23       5.6364       [2,2,4,5,5,6,6,7,7,8,10]
%     24-38      5.5556       [2,4,5,5,6,6,7,7,8]
%     39-53      5.7143       [4,5,5,6,6,7,7]
%     54-69      5.8000       [5,5,6,6,7]
%     70-84      5.6667       [5,6,6]
%     85-100     6.0000       [6]
%     -----------------------------------------------------
%
%    It can see that such a simple trimmed mean, involves to remove 
%    according to the calculus:
%
%                v = (ceil(([LP:1:HP]/100)*n))/2,
%
%    values on both sample ends: a series of two unique values.
%    
%    From the order statistics x_(i),
%
%                m = mean(x_(v+1):x_(n-v)).
%
%   Where:LP, lower percent;HP, higher percent;n, sample size.
%
%   Syntax: function fractrimmean(x,percent,dim)
%
%     Inputs:
%           x � data vector or matrix
%     percent - scalar between 0 and 100
%         dim - works along the dimension of X [columns=1(default);rows=2]
%
%     Output:
%           m - trimmed mean
%
%   Example 1. Suppose we have the sample: 50,98,82,23,46,40,63,52,92,54;
%      and we want to compute the 13% trimmed mean (26%). Step-by-step it 
%      is,
%   
%   x_(i)=23,40,46,50,52,54,63,82,92,98; n=10; g=[10*0.13]=[1.3]=1.0, and
%   r=1.3-1.0=0.3
%
%   T(0.13) = ((1 - 0.3)(x_(1+1) + x_(10-1) + i=1+2_S_10-1-1 x_(i))/10(1-2(0.13))
%           = ((0.7)(x_(2) + x_(9) + i=3_S_8 x_(i))/10(1-0.26)
%           = ((0.7)(x_(2) + x_(9) + i=3_S_8 x_(i)/7.4
%           = ((0.7)(40 + 92) + (46+50+52+54+63+82))/7.4
%           = ((0.7)(132) + 347)/7.4
%           = (92.4 + 347)/7.4
%           = 439.4/7.4
%           = 59.3784
%
%   By using this Matlab function,
%
%   Input data:
%   x=[50,98,82,23,46,40,63,52,92,54];
%
%   Calling on Matlab the function: 
%                fractrimmean(x,26)
%
%   Answer is:
%
%             59.3784
%
%   Example 2. Given the marix,
%
%           1     4     2     6     8
%           3     4     9     0    10
%           8   NaN     2     4   NaN
%           1     3   NaN     7     3
%
%   we want to compute the 15% trimmed row means (30%).
%
%   Input data:
%   x=[1 4 2 6 8;3 4 9 0 10;8 NaN 2 4 NaN;1 3 NaN 7 3];
%
%   Calling on Matlab the function: 
%                fractrimmean(x,30,2)
%
%   Answer is:
%
%             4.0714
%             5.2857
%             4.5238
%             3.2857
%
%   Example 3. From the example given by Horn and Kafadar (2002). The
%      20%-trimmed mean (40%) is,
%
%   Input data:
%   x=[46.7 72.9 79.6 83.6 80.7 60.3 79.0 64.8 49.6 54.7 71.8 49.1 103.9
%   51.6 81.6];
%
%   Calling on Matlab the function: 
%                fractrimmean(x,40)
%
%   Answer is:
%
%            68.3778
%
%   Created by A. Trujillo-Ortiz and R. Hernandez-Walls
%             Facultad de Ciencias Marinas
%             Universidad Autonoma de Baja California
%             Apdo. Postal 453
%             Ensenada, Baja California
%             Mexico.
%             atrujo@uabc.edu.mx
%
%   Copyright (C)  Sept. 14, 2011. 
%   
%   To cite this file, this would be an appropriate format:
%   Trujillo-Ortiz, A. and R. Hernandez-Walls. (2011). fractrimmean:
%      Fractional trimmed mean. A MATLAB file. [WWW document]. URL
%      http://www.mathworks.com/matlabcentral/fileexchange/32953-fractrimmean
% 
%  ---We thank Drs. Paul S. Horn (Department of Mathematical Sciences,
%  University of Cincinnati) and Karen Kafadar (Department of Statistics,
%  Indiana University Bloomington) for provided us a hard copy of their 
%  paper to be possible this work.---
%
%  References:
%  Horn, P. S. and Kafadar, K. (2002), Trimming and Winsorization.
%             In: Encyclopedia of Environmetrics, Vol. 4 (Abdel H. 
%             El-Shaarawi and Walter W. Piegorsch, editors), Chichester:
%             John Wiley and Sons., pp. 2264-2267.
%  Lim, T. S. and Loh, W. Y. (1996), A comparison of tests of equality of
%             variances. Computational Statistics and Data Analysis, 
%             22:287-301.
%
if nargin < 2,
    error('stats:fractrimmean:TooFewInputs', 'FRACTRIMMEAN requires two input arguments.');
elseif percent > 100 || percent < 0,
    error('stats:fractrimmean:InvalidPercent', 'PERCENT must be between 0 and 100.');
end
if nargin == 2, 
    dim = 1; %default 
end
p = percent/200;
[sr sc] = size(x);
if sr == 1 || sc == 1,
    x = sort(x(~isnan(x)));
    if p == 0,
        m = mean(x);
        s = std(x);
    elseif p == 0.5;
        m = median(x);
         s = std(x);
    else
        n = length(x);
        k = n*p;
        r = k-floor(k); %fractional remainder part of k
        g = k-r; %integer part of k
        C = n*abs(1-2*p); %retained number of observations after trimming p on each side
        if C >= 1,
            A = (1-r)*(x(g+1)+x(n-g)); %sum of the bottom and top end trimming
            B = sum(x(g+2:n-g-1)); %sum of the values between the bottom and top end trimming
            m = (A + B)/C;
            s = std(x(g+2:n-g-1));
        elseif C < 1,
            if x(g+1) == x(n-g);
                m = median(x);
                 s = std(x);
            end
        end
    end
else
    if dim == 1; %column mode
        if p == 0,
            m = nanmean(x,dim);
            s = std(x,dim,'omitnan');
            
        elseif p == 0.5;
            m = nanmedian(x,dim);
            s = std(x,dim,'omitnan');
        else
            for j = 1:sc,
                d(j).N = (~isnan(x(:,j)));
                d(j).o = x(:,j);
                d(j).x = sort(d(j).o(d(j).N));
                d(j).n = length(d(j).x);
                d(j).k = d(j).n * p;
                d(j).r = d(j).k - floor(d(j).k); %fractional remainder part of k
                d(j).g = d(j).k - d(j).r; %integer part of k
                d(j).C = d(j).n *(1-2*p); %retained number of observations after trimming p on each side
                if d(j).C >= 1,
                    d(j).A = (1 - d(j).r)*(d(j).x(d(j).g + 1) + d(j).x(d(j).n - d(j).g)); %sum of the bottom and top end trimming
                    d(j).B = sum(d(j).x(d(j).g + 2:d(j).n - d(j).g - 1)); %sum of the values between the bottom and top end trimming
                    dA=cat(1,d.A);dB=cat(1,d.B);dC=cat(1,d.C);
                    m = [(dA + dB)./dC]';
                    s =   std(d(j).x(d(j).g + 2:d(j).n - d(j).g - 1));
                elseif (C < 1),
                    if d(j).x(d(j).g + 1)==d(j).x(d(j).n - d(j).g);
                        m = nanmedian(x,dim);
                         s = std(x,dim,'omitnan');
                    end
                end
            end
        end
    elseif dim == 2; %row mode
        if p == 0,
            m = nanmean(x,dim);
            s = std(x,dim,'omitnan');
        elseif p == 0.5;
            m = nanmedian(x,dim);
            s = std(x,dim,'omitnan');
        else
            for i = 1:sr,
                d(i).N = (~isnan(x(i,:)));
                d(i).o = x(i,:);
                d(i).x = sort(d(i).o(d(i).N));
                d(i).n = length(d(i).x);
                d(i).k = d(i).n * p;
                d(i).r = d(i).k - floor(d(i).k); %fractional remainder part of k
                d(i).g = d(i).k - d(i).r; %integer part of k
                d(i).C = d(i).n *(1-2*p); %retained number of observations after trimming p on each side
                if d(i).C >= 1,
                    d(i).A = (1 - d(i).r)*(d(i).x(d(i).g + 1) + d(i).x(d(i).n - d(i).g)); %sum of the bottom and top end trimming
                    d(i).B = sum(d(i).x(d(i).g + 2:d(i).n - d(i).g - 1)); %sum of the values between the bottom and top end trimming
                    dA=cat(1,d.A);dB=cat(1,d.B);dC=cat(1,d.C);
                    m = [(dA + dB)./dC];
                    s = std(d(i).x(d(i).g + 2:d(i).n - d(i).g - 1));
                elseif (C < 1),
                    if d(i).x(d(i).g + 1)==d(i).x(d(i).n - d(i).g);
                        m = nanmedian(x,dim);
                        s = std(x,dim,'omitnan');
                    end
                end
            end
        end
    else dim > 2; %no mode
        x
    end
end
return,
