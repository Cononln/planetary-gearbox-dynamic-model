function tests=test_frozen_phase_reference
tests=functiontests(localfunctions);
end

function testCommonKnownSpeed(testCase)
fs=1000; t=(0:9999)'/fs; phi=2*pi*60*t+.3*sin(2*pi*.4*t);
z=complex(zeros(numel(t),3,4));
for h=1:4, for j=1:3, z(:,j,h)=exp(1i*(h*phi+.4*j)); end, end
out=pg_common_phase_frozen(z,fs,1:4,60,84);
ix=201:numel(t)-200; err=out.meshPhase(ix)-phi(ix); err=err-mean(err);
verifyLessThan(testCase,rms(err),.002);
verifyTrue(testCase,out.monotonic);
end

function testFrozenUnknownOffset(testCase)
[z,theta]=fixture(0); model=pg_frozen_path_fit(z,theta,1:3,-12:12,ones(1,25),1e-5);
[zt,tt]=fixture(.73); original=model;
[offset,~]=pg_frozen_path_offset(model,zt,tt);
[y,p]=pg_frozen_path_predict(model,zt,tt,offset);
verifyGreaterThan(testCase,pg_phase_resultant(y,1:3,model.scale),.999);
verifyLessThan(testCase,max(abs(abs(y(:))-abs(zt(:)))),1e-12);
verifyLessThan(testCase,max(abs(abs(p(:))-1)),1e-12);
verifyEqual(testCase,model,original);
end

function testIdentityAndMissingSignal(testCase)
[z,theta]=fixture(0); model=pg_frozen_path_fit(z,theta,1:3,[],[],.002);
y=pg_frozen_path_predict(model,z,theta,1);
verifyEqual(testCase,y,z);
zz=zeros(size(z)); m=pg_frozen_path_fit(zz,theta,1:3,-3:3,ones(1,7),.002);
[yy,p]=pg_frozen_path_predict(m,zz,theta,0);
verifyTrue(testCase,all(isfinite(p(:)))); verifyEqual(testCase,yy,zz);
verifyTrue(testCase,isnan(pg_phase_resultant(zz,1:3,m.scale)));
end

function testRetainedDifferentialAndEnergy(testCase)
[z,theta]=fixture(0); raw=z-mean(z,2);
model=pg_frozen_path_fit(z,theta,[1,2],-6:6,ones(1,13),.002);
y=pg_frozen_path_predict(model,z,theta,.2);
verifyEqual(testCase,y(:,3,:),z(:,3,:));
verifyEqual(testCase,z-mean(z,2),raw);
e=sqrt(abs(mean(y,2)).^2+mean(abs(y-mean(y,2)).^2,2));
verifyLessThan(testCase,max(abs(e(:)-reshape(sqrt(mean(abs(z).^2,2)),[],1))),1e-12);
end

function testStrictThreefoldAmbiguity(testCase)
[z,theta]=fixture(0); m=pg_frozen_path_fit(z,theta,1:3,-12:3:12,ones(1,9),.002);
[~,a]=pg_frozen_path_predict(m,z,theta,.1);
[~,b]=pg_frozen_path_predict(m,z,theta,.1+2*pi/3);
verifyLessThan(testCase,max(abs(a(:)-b(:))),1e-11);
end

function [z,theta]=fixture(offset)
t=(0:11999)'/1000; theta=2*pi*2*t;
z=complex(zeros(numel(t),3,4));
for h=1:4
    for j=1:3
        path=.4*(j-1)+.3*(j-1)*sin(theta+offset+.3*h)+ ...
            .12*(j-1)*cos(3*(theta+offset));
        z(:,j,h)=(1+.12*cos(theta+.2*j)).*exp(1i*(h*84*theta+path));
    end
end
end
