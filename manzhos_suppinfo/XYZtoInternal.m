function value = XYZtoInternal(x) % takes and retuns matrices Npts x Ncoords
C = x(:,1:3);
O = x(:,4:6);
H1 = x(:,7:9);
H2 = x(:,10:12);

r1vec = H1-C;
r1 = sqrt(sum(r1vec.^2,2));
r2vec = H2-C;
r2 = sqrt(sum(r2vec.^2,2));
r3vec = O-C;
r3 = sqrt(sum(r3vec.^2,2));

costheta1 =sum(r1vec.*r3vec,2)./(r1.*r3);
costheta2 =sum(r2vec.*r3vec,2)./(r2.*r3);
r1r3normvec = cross(r1vec,r3vec);
r1r3norm = sqrt(sum(r1r3normvec.^2,2));
r2r3normvec = cross(r2vec,r3vec);
r2r3norm = sqrt(sum(r2r3normvec.^2,2));
cosphi = sum(r1r3normvec.*r2r3normvec,2)./(r1r3norm.*r2r3norm);

value = [r3 r1 r2 acos([costheta1 costheta2]) acos(cosphi)];
ix=sign(sum(r1r3normvec.*r2vec,2));
value(ix>0,6)=2*pi-value(ix>0,6);
end