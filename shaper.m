function y = shaper(f,v)

% y=shaper(f,v) receives a vector v which contains a reshaped version of the
% matrix that contains the unit vectors along which each one of the
% muscles exerts its force. Also it receives a vector f, which contains the
% magnitude of the contraction forces exerted by each of the 10 muscles.
% The output y consists of the forces to be applied in each of the 11 arm
% insertion points.

% This shaper is to be used with versions 10a and 10b.

v = reshape(v,10,3);    % each row corresponds to a muscle

forces = bsxfun(@times,f,v)';    
                   % we use a single dimension to address this vector, so
                   % we want each vector to be in a column

% y(1:6) = forces(1:6);   % shoulder to arm muscles
% y(5) = y(5) + forces(7); % that damned shared insertion point
% y(7) = -forces(9);        % brachialis muscle on arm
% y(8) = -forces(10);      % triceps on arm
% y(9) = forces(8);       % biceps on forearm
% y(10) = forces(9);        % brachialis on forearm
% y(11) = forces(10);     % triceps on forearm

y(1:18) = forces(1:18);     % shoulder to arm muscles
y(13:15) = y(13:15) + forces(19:21);  % the shared insertion point
y(19:21) = -forces(25:27);   % brachialis muscle on arm
y(22:24) = -forces(28:30);      % triceps on arm
y(25:27) = forces(22:24);       % biceps on forearm
y(28:30) = forces(25:27);        % brachialis on forearm
y(31:33) = forces(28:30);     % triceps on forearm