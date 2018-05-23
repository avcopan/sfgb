function value = XYZtoRij(x) % takes and retuns matrices Npts x Ncoords
C = x(:,1:3);
O = x(:,4:6);
H1 = x(:,7:9);
H2 = x(:,10:12);
CO = sqrt(sum((C-O).^2,2));  % norm(C-O);
CH1 = sqrt(sum((C-H1).^2,2));
CH2 =sqrt(sum((C-H2).^2,2));
OH1 = sqrt(sum((H1-O).^2,2));
OH2 = sqrt(sum((H2-O).^2,2));
H1H2 = sqrt(sum((H1-H2).^2,2));
value = [CO CH1 CH2 OH1 OH2 H1H2];
end
