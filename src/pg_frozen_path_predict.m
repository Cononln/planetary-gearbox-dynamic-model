function [corrected,phasor] = pg_frozen_path_predict(model,z,theta,offset)
%PG_FROZEN_PATH_PREDICT Frozen coefficients and a single angular offset.
if nargin<4, offset=0; end
assert(isscalar(offset) && isfinite(offset));
assert(size(z,1)==numel(theta) && size(z,2)==model.nSensor && size(z,3)==model.nHarmonic);
assert(all(isfinite(z(:))) && all(isfinite(theta)));
phasor=complex(ones(size(z))); corrected=z;
if model.identity, return; end
B=exp(1i*(theta(:)+offset)*model.orders);
for k=1:size(z,3)
    for e=2:numel(model.subset)
        p=B*model.coefficient{e-1,k};
        ok=abs(p)>1e-10; q=complex(ones(size(p))); q(ok)=p(ok)./abs(p(ok));
        j=model.subset(e); phasor(:,j,k)=q;
        corrected(:,j,k)=z(:,j,k).*conj(q);
    end
end
end
