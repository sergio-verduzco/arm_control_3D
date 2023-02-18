function y = piece5(a,c,p1,p2)

% y = piece5(a,c,p1,p2) is a version of piecewise5 where the required data
% is passed as arguments, and four values are returned. The first value
% indicates whether the muscle should be wrapped on the line between p1 and
% p2. The rest of the values are the point where the muscle intersects the
% restricted bending line.
% This program is called from geometry4 and params10c.

% all the required data
% a = [-10 -5 -12];
% c = [-1 -1 -5];
% p1 = [-3 -4 -7];
% p2 = [-2 -4 2];


% preliminary calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = p2 - p1;    % b points from p1 to p2 
b = b/norm(b); % b should be a unit vector;
r = p1 - sum(p1.*b)*b;    % shortest line from origin to bending line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the calculations about the side +++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
v = p1 - a - sum((p1-a).*b)*b; % shortest vector from a to bending line
s1 = b(2)*v(3) - v(2)*b(3); % s is the r vector from the back of (49)
s2 = v(1)*b(3) - b(1)*v(3);
s3 = b(1)*v(2) - v(1)*b(2);
gam = sum((c-r).*[s1 s2 s3]); % gamma value from back of (49)

if gam > 0 % if muscle should wrap around bending line
    
    % calculating k and p3 ==============================================
    %====================================================================
    ar = a-r;   % I shift so b can be considered to intersect the origin
    cr = c-r;
    f1 = sum(ar.*ar);  % these are the calculations from (50)
    f2 = sum(cr.*cr);
    f3 = sum(ar.*b);
    f4 = sum(cr.*b);
    ap = f2 + f3^2 - f1 - f4^2;
    bp = 2*(f1*f4 - f2*f3 + f3*f4*(f4-f3));
    cp = f2*f3^2 - f1*f4^2;
    
    % taking the root with smallest sum of distances
    % (wish I had figured a better way to decide this)
    k1 = (-bp + sqrt(bp^2 - 4*ap*cp))/(2*ap);
    k2 = (-bp - sqrt(bp^2 - 4*ap*cp))/(2*ap);
    n1 = norm(ar-k1*b)+norm(cr-k1*b);
    n2 = norm(ar-k2*b)+norm(cr-k2*b);
    if n1 < n2
        k = k1;
    else
        k = k2;
    end
    
    % limiting p3 between p1 and p2
    [biggs ind] = max(abs(b)); % this is because b may have zero entries
    k1 = (p1(ind) - r(ind))/b(ind);
    k2 = (p2(ind) - r(ind))/b(ind);
    k = min(k,max(k1,k2));
    k = max(k,min(k1,k2)); % k is now restricted to the segment p1 to p2
    p3 = r + k*b;   % the point of intersection between muscle and bending line
    y = [1 p3(1) p3(2) p3(3)];
    %====================================================================
else  % muscle doesn't wrap
    y = [0 0 0 0];
end
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
