function y = lowpass(t,e)

% y = lowpass(t,e) receives a time t, and a vector e. It uses two global
% variables in order to obtain a discrete approximation to a low-pass
% filter (like a simple RC filter: T v' = v_in - v).
% The two global variables are the filtered value of the input signal
% (called e_low), and the time when e_low was last updated.

global last_t e_low LPfactor

if t > last_t
    dif = t - last_t;
    last_t = t;
    
    e_low = e_low + dif*LPfactor*(e - e_low);
end

y = e_low;