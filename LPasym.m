function y = LPasym(t,e)

% y = LPasym(t,e) receives a time t, and a vector e. It uses two global
% variables in order to obtain a discrete approximation to a low-pass
% filter (like a simple RC filter: T v' = v_in - v).
% The two global variables are the filtered value of the input signal
% (called e_low), and the time when e_low was last updated.
% This version implements an "asymmetric filter", meaning that when e is
% larger than e_low, low-pass filtering will be in effect, but when e_low
% is larger than e, higher frequencies will be allowed.

global last_t e_low LPfactor

if t > last_t
    dif = t - last_t;
    last_t = t;
    
    e_low = e_low + dif*LPfactor*max(e - e_low,0) + ...
            10*dif*LPfactor*min(e - e_low,0);
end

y = e_low;