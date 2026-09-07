function [weighted,unweighted] = pg_phase_resultant(z,subset,scale)
%PG_PHASE_RESULTANT Instantaneous inter-sensor alignment, not temporal PLV.
z=z(:,subset,:); a=abs(z); u=z./max(a,eps);
r=a./max(scale(:,subset,:),eps); w=r.^2./(1+r.^2);
den=sum(w,2); c=abs(sum(w.*u,2))./max(den,eps);
valid=den>1e-12;
if any(valid(:)), weighted=mean(c(valid)); else, weighted=NaN; end
c0=abs(mean(u,2)); active=sum(a,2)>eps;
if any(active(:)), unweighted=mean(c0(active)); else, unweighted=NaN; end
end
