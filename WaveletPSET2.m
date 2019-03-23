%% Michael Lendino Wavelets PSET 2
clc;
clear all;
%% Graphing unit cells for 2D and 3D lattices
%see function at bottom unitlattice()
%% Plotting points in lattices and cosets of lattices
%creating all lattices for this problem
V = [1 1; -1 1];
V1 = [1 -1/2; 0 sqrt(3)/2];
Vrecip = (V^-1)';
V1recip = (V1^-1)';
M = [2 0; 1 2];
U = V*M;
U1 = V1*M;
Urecip = (U^-1)';
U1recip = (U1^-1)';
%creating vector for lexicrographic indexing over the grid from (-3,3) in
%both x and y
n = -3:3;
[X,Y] = meshgrid(n);
%just need X and Y as a col vec
colX = X(:);
colY = Y(:);
vec2 = [transpose(colX); transpose(colY)];
%generating the coordinates for the lattice r=Vn
r = V*vec2;
rV = Vrecip*vec2;
rU = U*vec2;
rUrecip = Urecip*vec2;
r1 = V1*vec2;
rV1 = V1recip*vec2;
rU1 = U1*vec2;
rU1recip = U1recip*vec2;
%selecting the vector r in the lattice of V but not in the lattice of the
%unit cell (returns data in first component but not in second component)
rr = setdiff(rV.',r.','rows').';
rVU = setdiff(r.',rU.','rows').';
rVU1 = setdiff(r1.',rU1.','rows').';
rVUrecip = setdiff(rV.', rUrecip.', 'rows').';
rVU1recip = setdiff(rV1.', rU1recip.', 'rows').';

figure('Name','Quincunx Lattice','NumberTitle','off');
%plotting lattice
plot(r(1,:),r(2,:), 'x');
hold on;
%plotting U=VM
plot(rU(1,:),rU(2,:), 'o');
grid on;
title('Quincunx lattice and lattice corresponding to U=VM');
plot(rVU(1,:),rVU(1,:), '*');
%Above figures confirm that the lattice generated by U is a sublattice of
%the lattice generated by V (quincunx)

%now we wish to plot the lattices corresponding to the other V matrix, the
%other U matrix, and the points in the coset in the pset
figure('Name','Hex Lattice','NumberTitle','off');
%plotting lattice
plot(r1(1,:),r1(2,:), 'x');
hold on;
%plotting U=VM
plot(rU1(1,:),rU1(2,:), 'o');
grid on;
title('Hex lattice and lattice corresponding to U=VM');
plot(rVU1(1,:),rVU1(1,:), '*');
%Above figures confirm that the lattice generated by U is a sublattice of
%the lattice generated by V (hex)

%now we wish to plot the reciprocal lattice of the quincunx lattice and the reciprocal lattice
%of U, and the points in the coset in the pset
figure('Name','Reciprocal lattice for Quincunx','NumberTitle','off');
%plotting lattice
plot(rV(1,:),rV(2,:), 'x');
hold on;
%plotting U=VM
plot(rUrecip(1,:),rUrecip(2,:), 'o');
grid on;
title('Reciprocal lattices for Quincunx');
plot(rVUrecip(1,:),rVUrecip(1,:), '*');

%Above figures confirm that the reciprocal lattice generated by V is a sublattice of
%the reciprocal lattice generated by U (quincunx)

%now we wish to plot the reciprocal lattice of the hex lattice and the reciprocal lattice
%of U, and the points in the coset in the pset
figure('Name','Reciprocal lattice for Hex','NumberTitle','off');
%plotting lattice
plot(rV1(1,:),rV1(2,:), 'x');
hold on;
%plotting U=VM
plot(rU1recip(1,:),rU1recip(2,:), 'o');
grid on;
title('Reciprocal lattices for Hex');
plot(rVU1recip(1,:),rVU1recip(1,:), '*');

%Above figures confirm that the reciprocal lattice generated by V is a sublattice of
%the reciprocal lattice generated by U (hex)
%% Plotting unit cells
%we wish to plot the unit cells of U and V for the quincunx and the hex
%cases and their reciprocal cells

%plotting unit cell for V quincunx case
unitlattice(V);

%plotting unit cell for U quincunx case
unitlattice(U);

%plotting unit cell for reciprocal lattice for V quincunx case
unitlattice(Vrecip);

%plotting unit cell for reciprocal lattice for U quincunx case
unitlattice(Urecip);

%plotting unit cell for V hex case
unitlattice(V1);

%plotting unit cell for U hex case
unitlattice(U1);

%plotting unit cell for reciprocal lattice for V hex case
unitlattice(V1recip);

%plotting unit cell for reciprocal lattice for U hex case
unitlattice(U1recip);

%% Face centered orthorhombic (FCO) lattice
Vfco = [1 0 1; -1 -1 1; 0 -1 0];
Mfco = [2 0 0; 1 2 0; 0 0 1];
Ufco = Vfco*Mfco;
Vfcorecip = (Vfco^-1)';
Ufcorecip = (Ufco^-1)';
%repeat problem 2 (over the lexicographic range of -2 to 2 plotting V and U and confiriming the lattice generated
%by U is a sublattice of the lattice generated by V, and that the
%reciprocal lattice generated by V is a sublattice of the reciprocal
%lattice generated by U and by plotting the relevant cosets)
n1 = -2:2;
[X1,Y1,Z1] = meshgrid(n1);
%just need X, Y, and Z as a col vec for lexicographic indexing
colX1 = X1(:);
colY1 = Y1(:);
colZ1 = Z1(:);
vec21 = [transpose(colX1); transpose(colY1); transpose(colZ1)];

rVfco = Vfco*vec21;
rUfco = Ufco*vec21;
rVfcorecip = Vfcorecip*vec21;
rUfcorecip = Ufcorecip*vec21;

rVUfco = setdiff(rVfco.',rUfco.','rows').';
rVUfcorecip = setdiff(rVfcorecip.',rUfcorecip.','rows').';

figure('Name','Face Centered Orthorhombic Lattice','NumberTitle','off');
%plotting lattice
plot3(rVfco(1,:,:),rVfco(2,:,:),rVfco(3,:,:), 'x');
hold on;
%plotting U=VM
plot3(rUfco(1,:,:),rUfco(2,:,:), rUfco(3,:,:),'o');
grid on;
title('Face Centered Orthorhombic Lattice and lattice corresponding to U=VM');
plot3(rVUfco(1,:,:),rVUfco(2,:,:), rVUfco(3,:,:),'*');
%Above figure confirm that the lattice generated by U is a sublattice of
%the lattice generated by V (FCO)

figure('Name','Face Centered Orthorhombic Lattice Reciprocal Lattice','NumberTitle','off');
%plotting lattice
plot3(rVfcorecip(1,:,:),rVfcorecip(2,:,:),rVfcorecip(3,:,:), 'x');
hold on;
%plotting U=VM
plot3(rUfcorecip(1,:,:),rUfcorecip(2,:,:), rUfcorecip(3,:,:),'o');
grid on;
title('Face Centered Orthorhombic Lattice and lattice corresponding to U=VM');
plot3(rVUfcorecip(1,:,:),rVUfcorecip(2,:,:), rVUfcorecip(3,:,:),'*');
%Above figure confirm that the reciprocal lattice generated by V is a sublattice of
%the reciprocal lattice generated by U (FCO) 

%% Unit cell FCO
%Plot the unit cells and reciprocal unit cells for the face centered
%orthorhombic lattice 

%plotting unit cell for V FCO case
unitlattice(Vfco);

%plotting unit cell for U FCO case
unitlattice(Ufco);

%plotting unit cell for reciprocal lattice for V FCO case
unitlattice(Vfcorecip);

%plotting unit cell for reciprocal lattice for U FCO case
unitlattice(Ufcorecip);
%% Matrix relations and lattices
Uprime = [14 2; 2 2];
Vprime = [3 1; 0 1];
M6 = [4 0; 2 2];
%note that Uprime = Vprime*M6, also the smith form of M is M=E1muE2 where
E1 = [2 1; 1 0];
mu = diag([2 4]);
E2 = [1 1; 0 -1];
%Confirm these relationships including the fact that E1 and E2 are
%unimodular
uni1 = det(E1);
uni2 = det(E2);
M61 = E1*mu*E2;
%E1 and E2 are unimodular because they are both square integer matirces
%having determinants \pm 1, in this case -1. Moreover, we have confirmed
%that M=E1muE2

n = -5:5;
[X2,Y2] = meshgrid(n);
%just need X and Y as a col vec
colX2 = X2(:);
colY2 = Y2(:);
vec22 = [transpose(colX2); transpose(colY2)];

rUprime = Uprime*vec22;
rVprime = Vprime*vec22;
%create images of the lattices generate by Uprime and Vprime
figure('Name','Lattice generated by U` ','NumberTitle','off');
plot(rUprime(1,:),rUprime(2,:), 'x');
hold on;

figure('Name','Lattice generated by V` ','NumberTitle','off');
plot(rVprime(1,:),rVprime(2,:), 'x');
hold on;

%Renormalize U' and V' via U=U'E2^-1 and V=V'E1
U6 = Uprime*(E2^-1);
V6 = Vprime*E1;
U6recip = (U6^-1)';
V6recip = (V6^-1)';
%Check that U=Vmu, then create images of the lattice generated by U and the
%lattice generated by V
U6check = V6*mu;
%confirmed as we get U6

rU6 = U6*vec22;
rV6 = V6*vec22;
rU6recip = U6recip*vec22;
rV6recip = V6recip*vec22;
figure('Name','Lattice generated by renormalized U` and lattice generated by renormalized V` ','NumberTitle','off');
plot(rU6(1,:),rU6(2,:), 'o');
hold on;
plot(rV6(1,:),rV6(2,:), 'x');
hold on;
%we also wish to plot * at points in the coset of Lu corresponding to the
%lexicographic offset index [1 3]^T
x = [1;3];
Vx = V6*x;
Y = rU6 +Vx;
plot(Y(1,:),Y(2,:), '*');
title('Lattice generated by renormalized U` and lattice generated by renormalized V`')
%To find the index set we do a cartesian product of the two sets of points
%until the lattice repeats so looking at where the lattice of U overlaps
%with the lattice of V and starting at the bottom left corner of this
%overlap, calling it 0 and incrementing by one moving to the left and upwards, we see the index set will be {0,1}x{0,1,2,3} or 
%{(0,0),(0,1),(0,2),(0,3),(1,0),(1,1),(1,2),(1,3)}
%Now for U and V (i.e. what we just plotted) we wish to plot their unit
%cells, their reciprocal lattices, and the unit cells for their reciprocal
%lattices)

%plotting unit cell for U
unitlattice(U6);

%plotting unit cell for V
unitlattice(V6);

%plotting reciprocal unit cell for U
unitlattice(U6recip);

%plotting reciprocal unit cell for V
unitlattice(V6recip);

figure('Name','Reciprocal lattice generated by renormalized U` and reciprocal lattice generated by renormalized V` ','NumberTitle','off');
plot(rU6recip(1,:),rU6recip(2,:), 'o');
hold on;
plot(rV6recip(1,:),rV6recip(2,:), 'x');
hold on;
title('Reciprocal lattice generated by renormalized U` and reciprocal lattice generated by renormalized V`')

%Find an intermediate lattice Lw such that Lu is a sublattice of Lw and Lw
%is a sublattice of Lv, plot points in that lattice and its recirpocal
%lattice, check the tower of reciprocal lattices are related in the reverse
%sense, to find this intermediate lattice we do what we do in the next
%problem and pick some lattice in the decomposition so in this case our M
%matrix is mu so we can see det(mu) = 8 = 2x2x2 so we can write mu as a
%product of diagonal matrices each with determinant 2 so z= [1 0; 0 2] and a
%matrix from the generating sequence of matrices would be like U6*z
w = [1 0; 0 2];
ww = V6*w;
wrecip = (ww^-1)';
rw = ww*vec22;
rwrecip = wrecip*vec22;
%plot points in intermediate lattice
figure('Name','Intermediate lattice ','NumberTitle','off');
plot(rw(1,:),rw(2,:), 'x');
hold on;

%now we have to check that the tower of reciprocal lattices are related in
%the reverse sense, so plot points in reciprocal lattice of intermediate
%lattice, and reciprocal lattice of other two lattices
figure('Name','Reciprocal Lattice Tower ','NumberTitle','off');
plot(rwrecip(1,:),rwrecip(2,:), '*');
hold on;
plot(rU6recip(1,:),rU6recip(2,:), 'o');
hold on;
plot(rV6recip(1,:),rV6recip(2,:), 'x');
title('Reciprocal lattice tower')

%% Functions

function unitlattice(V)

dimV = size(V);

if dimV(1) == 2;
side1 = [-0.5 0.5];

x1 = -1/2;
x2 = 1/2;
y1 = -1/2;
y2 = 1/2;
z1 = -1/2;
z2 = 1/2;
x = [x1, x2, x2, x1, x1];
y = [y1, y1, y2, y2, y1];
a = [x;y];
unit = zeros(2,5);

 for i = 1:5
    unit(:,i) = V*a(:,i); 
 end

figure('Name','Unit Cell in 2D','NumberTitle','off');
title('Unit Cell in 2D');
plot(unit(1,:),unit(2,:));
grid on;
xlim([-2 2]);
ylim([-2 2]);

else 
x1 = -1/2;
x2 = 1/2;
y1 = -1/2;
y2 = 1/2;
z1 = -1/2;
z2 = 1/2;

x1 = [x1,x2,x2,x1,x1,x2,x2,x1,x1,x1,x1,x2,x2,x1,x1,x2,x2,x2];
y1 = [y1,y1,y2,y2,y2,y2,y1,y1,y2,y2,y1,y1,y1,y1,y1,y1,y2,y2];
z1 = [z1,z1,z1,z1,z2,z2,z2,z2,z2,z1,z1,z1,z2,z2,z1,z1,z1,z2];
b = [x1;y1;z1];
unit1 = zeros(3,18);

for i = 1:18
    unit1(:,i) = V*b(:,i);
end
    
figure('Name','Unit Cell in 3D','NumberTitle','off');
title('Unit Cell in 3D');
plot3(unit1(1,:), unit1(2,:), unit1(3,:));
grid on;
xlim([-2 2]);
ylim([-2 2]);
zlim([-2 2]);
end
end



