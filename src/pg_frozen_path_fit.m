function model = pg_frozen_path_fit(z,theta,subset,orders,penalty,ridge)
%PG_FROZEN_PATH_FIT Healthy-only periodic relative phasor calibration.
% Fixed star topology avoids per-record unwrap/gauge changes at prediction.
assert(size(z,1)==numel(theta) && all(isfinite(z(:))));
assert(all(isfinite(theta)) && numel(orders)==numel(penalty));
assert(numel(unique(subset))==numel(subset) && all(subset>=1 & subset<=size(z,2)));
model.subset=subset; model.orders=orders(:).';
model.penalty=penalty(:).'; model.ridge=ridge;
model.nSensor=size(z,2); model.nHarmonic=size(z,3);
model.scale=max(median(abs(z),1),eps);
model.coefficient=cell(numel(subset)-1,size(z,3));
model.identity=isempty(orders) || numel(subset)<2;
model.trainingSamples=size(z,1);
if model.identity, return; end
% Subsample only for fitting; no learned parameters depend on test data.
ix=1:8:size(z,1); B=exp(1i*theta(ix)*model.orders);
for k=1:size(z,3)
    for e=2:numel(subset)
        i=subset(1); j=subset(e);
        q=z(ix,j,k).*conj(z(ix,i,k)); q=q./max(abs(q),eps);
        ai=abs(z(ix,i,k))/model.scale(1,i,k);
        aj=abs(z(ix,j,k))/model.scale(1,j,k);
        w=sqrt(ai.^2./(1+ai.^2).*aj.^2./(1+aj.^2));
        if sum(w)<=eps
            coef=complex(zeros(numel(orders),1));
            zero=find(orders==0,1); assert(~isempty(zero)); coef(zero)=1;
        else
            A=B.*sqrt(w); y=q.*sqrt(w); N=A'*A;
            coef=(N+ridge*trace(N)/size(N,1)*diag(penalty))\(A'*y);
        end
        model.coefficient{e-1,k}=coef;
    end
end
end
