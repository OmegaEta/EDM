function isDetected = detected_merge(states,w)
%

w = exp(w);
idx = w > 0;
w = w(idx);
states = states(idx);
isDetected = [states.isDetected]*w;


end

