function [ X_sub, Y_sub ] = subsample( X_in, Y_in, nSubSamples, type )
% subsample a curve to get equal point spacing
% [ X_sub, Y_sub ] = subsample( X_in, Y_in, nSubSamples, type )


% check input data
if(~isequal(size(X_in),size(Y_in)))
    error('input vectors must be of same length');
end
if(size(X_in,2)~=1)
    error('input vectors must be row vectors');
end

% define linear or logarithmic metric
switch type
    case 'linear'
        x_metric = @(x) x;
        y_metric = @(x) x;
    case 'semilogy'
        x_metric = @(x) x;
        y_metric = @(x) log10(x);
    case 'semilogx'
        x_metric = @(x) log10(x);
        y_metric = @(x) x;
    case 'loglog'
        x_metric = @(x) log10(x);
        y_metric = @(x) log10(x);
    otherwise
        error('unknown subsampling type');
end

% initialization
X_sub = NaN(nSubSamples,1);
Y_sub = NaN(nSubSamples,1);
nSamples = length(X_in);

% calculate total length of curve
curve_length = 0;
for ii=1:nSamples-1
    curve_length = curve_length + sqrt( x_metric(abs(X_in(ii+1)-X_in(ii))).^2 + y_metric(abs(Y_in(ii+1)-Y_in(ii))).^2 );
end
average_dist = curve_length/(nSubSamples-1);

% subsampling
X_sub(1,1) = X_in(1,1);     % take first point as start
Y_sub(1,1) = Y_in(1,1);
sub_index = 2;
temp_dist = 0;
for ii=1:nSamples-1
    temp_dist = temp_dist + sqrt( x_metric(abs(X_in(ii+1)-X_in(ii))).^2 + y_metric(abs(Y_in(ii+1)-Y_in(ii))).^2 );
    if(temp_dist >= average_dist)
        temp_dist = 0;                      % reset temp distance
        X_sub(sub_index,1) = X_in(ii,1);    % take point
        Y_sub(sub_index,1) = Y_in(ii,1);
        sub_index = sub_index + 1;          % increase subsample index
        if(sub_index > nSubSamples)
            break
        end
    end
end

% always use last point
X_sub(nSubSamples,1) = X_in(nSamples,1);
Y_sub(nSubSamples,1) = Y_in(nSamples,1);


end

