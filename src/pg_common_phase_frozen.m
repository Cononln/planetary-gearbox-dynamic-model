function out = pg_common_phase_frozen(z,fs,h,fmesh,zr,referenceScale)
%PG_COMMON_PHASE_FROZEN Common phase with no path fitting or encoder input.
assert(size(z,3)==numel(h) && size(z,1)>2);
assert(all(isfinite(z(:))) && fs>0 && fmesh>0 && zr>0);
t=(0:size(z,1)-1)'/fs;
nominal=2*pi*fmesh*t;
a=abs(z);
if nargin<6, referenceScale=median(a,1); end
r=a./max(referenceScale,eps);
w=r.^2./(1+r.^2);
num=zeros(size(z,1)-1,1); den=num;
for j=1:size(z,2)
    for k=1:numel(h)
        inc=angle(z(2:end,j,k).*conj(z(1:end-1,j,k))* ...
            exp(-1i*h(k)*2*pi*fmesh/fs))/h(k);
        wt=min(w(2:end,j,k),w(1:end-1,j,k));
        num=num+wt.*inc; den=den+wt;
    end
end
inc=smoothEdge(num./max(den,eps),round(.0015*fs));
delta=smoothEdge([0;cumsum(inc)],round(.012*fs));
out.meshPhase=nominal+delta;
out.carrierPhase=out.meshPhase/zr;
out.frequencyHz=smoothEdge(gradient(out.meshPhase)*fs/(2*pi),round(.1*fs));
out.monotonic=all(diff(out.carrierPhase)>0);
out.nominalFrequencyHz=fmesh;
end

function y=smoothEdge(x,n)
n=2*floor(n/2)+1; p=floor(n/2);
y=conv([repmat(x(1),p,1);x;repmat(x(end),p,1)],ones(n,1)/n,'valid');
end
