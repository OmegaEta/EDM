function [GGIW,lik]=misdetectionGGIW(GGIW)
% If `detectionBern` or `detectionPPP` did not generate any measurements, the shape components are not updated.
lik=(GGIW.b/(GGIW.b+1))^GGIW.a;
lik=log(lik);

GGIW.b=GGIW.b+1;
end