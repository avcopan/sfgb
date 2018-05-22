clear all;
format shortG 
ndimxyz=12;
prec = 'double'

Vmax=15000               % cm-1
M = 15000;    Nsets=5;   % number of pts and number of sets to add into the collocation eq.
Ng = 15000;              % no.  of Gaussian basis functions
LHS = zeros(Ng,Ng,prec);
RHS = zeros(Ng,Ng,prec);
Nsample=2000000;         % number of Sobol points to use (should be large enough to select M*Ng collocation points after filter by PES values 

% reference zpe and frequencies
reference1997 = [5775.286288 1166.759293 1250.571992 1500.078837 1746.647599 2327.284157 2421.75511 2497.994011 2666.349468 2719.701761 2781.320106 2842.646376 2905.866486 2999.074342 3001.286221 3238.793038 3471.79816 ...
3484.413661 3586.533395 3674.870945 3742.096933 3825.05382 3887.222494 3937.668388 3940.375294 3995.732256 4033.154131 4058.896989 4085.553871 4163.882298 4164.451521 4194.304944 4250.487881 4253.081961 4336.660056 ...
4397.447747 4467.406823 4496.09804 4528.159465 4572.701793 4624.010791 4647.059591 4729.516266 4734.426384 4749.211937 4843.570876 4926.211024 4956.069902 4980.29346 4982.879596 5044.46328 5092.744198 5109.200325 5141.992836 ...
5154.220956 5176.075426 5194.832083 5210.57628 5244.906361 5275.119492 5317.714432 5319.735697 5327.846279 5358.648304 5389.556595 5411.848972 5418.792628 5433.632609 5461.101406 5489.869272 5494.224067 5531.813471 5543.668628 ...
5553.112447 5626.783873 5658.515495 5669.972685 5675.33035 5717.663885 5775.291654 5825.944142 5834.566287 5886.520125 5889.062861 5919.497273 5937.789733 6012.100863 6055.216617 6098.395261 6107.534823 6176.111069 6196.637335 ...
6197.871114 6220.671105 6243.272509 6277.561719 6287.482644 6321.818082 6332.397857 6441.198937]'; 
      
% equilibrium geometry
r1eq = 1.10064; % Ang
r2eq = 1.10064;
r3eq = 1.20296;
theta1eq = 121.65*pi/180;
theta2eq = 121.65*pi/180;
phieq = pi;

% Cartesian equilibrium coordinates in a.u.
Ceq = [0.0  0.0  0.0]*1.88973; 			
Oeq = [0.0  r3eq  0.0]*1.88973;
H1eq = [-sin(theta1eq)*r1eq  cos(theta1eq)*r1eq  0.0]*1.88973;
H2eq = [ sin(theta2eq)*r2eq  cos(theta2eq)*r2eq  0.0]*1.88973;
xeq=[Ceq Oeq H1eq H2eq];
m=[12.0*ones(1,3) 15.99491*ones(1,3)  1.00794*ones(1,3) 1.00794*ones(1,3)];

% Rij equilibrium coordinates in a.u.
COeq = norm(Ceq-Oeq);
CH1eq = norm(Ceq-H1eq);
CH2eq = norm(Ceq-H2eq);
OH1eq = norm(Oeq-H1eq);
OH2eq = norm(Oeq-H2eq);
H1H2eq = norm(H1eq-H2eq);
Req = [COeq CH1eq CH2eq OH1eq OH2eq H1H2eq];

% sample in internal coordinates, order: CO, CH1, CH2, H1CO H2CO, dihedral
try
    internal=dlmread('Sobol.dat');
    internal=internal(1:Nsample,:);
    message = 'existing Sobol.dat can be used for this Nsample'
catch
    message = 'existing Sobol.dat too short, recomputing Sobol sequence'
    internal=Sobol(6,Nsample);
    dlmwrite('Sobol.dat',internal);
end;
minmaxinternal= [1.03 1.50; 0.84 1.69; 0.84 1.69; 83 162; 83 162; 105 255];       % in Ang and degrees, for Vmax=15000
minmaxinternal(1:3,:)=minmaxinternal(1:3,:)*1.88973;                              % in Bohr
minmaxinternal(4:6,:)=minmaxinternal(4:6,:)*pi/180;                               % in radian
for i=1:6,
    internal(:,i)= ones(Nsample,1)*minmaxinternal(i,1)+ones(Nsample,1)*(minmaxinternal(i,2)-minmaxinternal(i,1)).*internal(:,i);
end;
% convert internals to Cartesians
internal_eq = [[r3eq r1eq r2eq]*1.88973 theta1eq theta2eq phieq];
internal = [internal_eq; internal];                                               % add min
x = internalToXYZ(internal);
Rij = XYZtoRij(x); 

Vall = 0.5*(PESH2CO(Rij/1.88973)+PESH2CO(Rij(:,[1 3 2 5 4 6])/1.88973));          % compute PES pts taking symmetr into account

MeanAbsFreq = [];
MeanRes = [];
FreqAccuracy = [];
for accum=1:Nsets,  % accmulates square matrices for eigs to alleviate memory limitations due to a large no. of pts
   
xcoords=[];
Rijcoords=[];
bondcoords=[];
V=[];
accepted = 0;
i=0;
while (accepted<M), 
    i=i+1;
    Vi = Vall(i+(accum-1)*M);    
    if (Vi < Vmax)&&(Vi<Vmax*1.2*rand),
        xcoords=[xcoords; x(i+(accum-1)*M,:)];
        Rijcoords=[Rijcoords; Rij(i+(accum-1)*M,:)];
        bondcoords=[bondcoords; internal(i+(accum-1)*M,:)];
        V=[V; Vi];
        accepted = accepted+1;
    end;        
end;
MinMaxV=minmax(V')

% basis data
if accum==1,
    G = bondcoords(1:Ng,:);   % centers of the Gaussian basis functions
end;
F = zeros(M, Ng, prec); 
F=basisGaussianH2CO(G, bondcoords);
[Npts Nfns] = size(F)

% Compute KEO(basis) and matrix TFik=KEO(fk(Ri)) using finite difference 2nd derivatives in Cartesian coordinates
step=1e-5;
TF=zeros(Npts, Nfns, prec);
for i=1:ndimxyz,
    dx=zeros(1,ndimxyz);
    dx(i)=step;    
    dx=ones(Npts,1)*dx;
    TF=TF-0.5/m(i)*(-30*F)/(12*step*step);
    Rijstep=XYZtoInternal(xcoords+2*dx);
    TF=TF-0.5/m(i)*(-basisGaussianH2CO(G, Rijstep))/(12*step*step);
    Rijstep=XYZtoInternal(xcoords+dx);
    TF=TF-0.5/m(i)*(16*basisGaussianH2CO(G, Rijstep))/(12*step*step);
    Rijstep=XYZtoInternal(xcoords-dx);
    TF=TF-0.5/m(i)*(16*basisGaussianH2CO(G, Rijstep))/(12*step*step);
    Rijstep=XYZtoInternal(xcoords-2*dx);
    TF=TF-0.5/m(i)*(-basisGaussianH2CO(G, Rijstep))/(12*step*step);  
end;
clear FD Rijstep

TF=TF*219474.63/1822.8886154; % divide by amu in au to bring KEO to au then x219474 to bring energy to cm-1 scale

L = 50; % number of levels
% square the matrix equation: accumulate squared matrices to reduce memory cost from big point sets
LHS = LHS+F'*(TF+diag(V)*F);
RHS = RHS+F'*F;
[Vec E] = eigs(LHS, RHS, L, 'sm');  
[Esorted ix] = sort(real(diag(E)));
Vec = Vec(:,ix(1:L));
E = real(E(:,ix(1:L)));
format shortG
freqs=[Esorted(1); Esorted(2:L)-Esorted(1)];
FreqAccuracy = [FreqAccuracy freqs freqs-reference1997(1:L)]

% check SE residual
MeanRes = [MeanRes mean(sqrt(sum(( (TF+diag(V)*F)*Vec-F*(Vec*E) ).^2,1)))]
MeanAbsFreq = [MeanAbsFreq mean(abs(freqs-reference1997(1:L)))]

end; % for accum











