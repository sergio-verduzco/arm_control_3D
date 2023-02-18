function y = shaper2(f,v)

% y=shaper(f,v) receives a vector v which contains a reshaped version of the
% matrix that contains the unit vectors along which each one of the
% muscles exerts its force. Also it receives a vector f, which contains the
% magnitude of the contraction forces exerted by each of the 11 muscles.
% The output y consists of the forces to be applied in each of the 13 arm
% insertion points.

% This shaper is to be used with version 10c.

v = reshape(v,11,3);    % each row corresponds to a muscle

forces = bsxfun(@times,f,v)';    
                   % we use a single dimension to address this vector, so
                   % we want each vector to be in a column

y(1:24) = forces(1:24); % shoulder to arm muscles (muscles 1-8 of (46))
y(25:27) = -forces(28:30);   % triceps on arm
y(28:30) = -forces(31:33);   % brachialis on arm
y(31:33) = forces(25:27);   % biceps on forearm
y(34:39) = forces(28:33);   % triceps and brachialis on forearm

