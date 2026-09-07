function [offset,gridScore] = pg_frozen_path_offset(model,z,theta)
%PG_FROZEN_PATH_OFFSET One record-wide nuisance offset from calibration only.
if model.identity || all(model.orders==0)
    offset=0; gridScore=NaN; return
end
ix=1:40:size(z,1); zz=z(ix,:,:); tt=theta(ix);
grid=(0:3:357)'*pi/180; gridScore=zeros(size(grid));
for i=1:numel(grid), gridScore(i)=objective(grid(i)); end
[~,best]=min(gridScore); step=3*pi/180;
offset=fminbnd(@objective,grid(best)-step,grid(best)+step,optimset('TolX',1e-5));
offset=mod(offset,2*pi);
    function loss=objective(a)
        y=pg_frozen_path_predict(model,zz,tt,a);
        loss=-pg_phase_resultant(y,model.subset,model.scale);
    end
end
