function value = internalToXYZ(internal)
[Nsample ndim]=size(internal);
ndimxyz=12;
x=zeros(Nsample,ndimxyz);
x(:,5)=internal(:,1);                      % y(O)

x(:,7)=-sin(internal(:,4)).*internal(:,2); % x(H1)
x(:,8)=cos(internal(:,4)).*internal(:,2);  % y(H1)

x(:,11)=cos(internal(:,5)).*internal(:,3); % y(H2)

ksi=sin(internal(:,5)).*internal(:,3);     % distance from H2 to the y axis
x(:,10)=ksi.*cos(pi-internal(:,6));        % x(H2)
x(:,12)=ksi.*sin(pi-internal(:,6));        % z(H2)
value=x;
end
