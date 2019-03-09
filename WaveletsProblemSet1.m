%% Michael Lendino Wavelets PSET 1 
clc;
clear all;
%% Computing Scaling and Wavelet Functions from Filter Coefficients
%Write code in MATLAB to set up the H matrix for arbitrary N
%we wish to generate the H matrix populated with appropriate
%indices n, this outputs the matrix for a given N, and indicates the
%indices that are within the support of Psi; see function at bottom of code
N = 3;
H = Hgen(N);

%Write code to compute phi(t), psi(t) for t dyadic integers of the form
%m/2^p for a prescribed p

%% Checking Paraunitary Property 
%There is a theoretical time domain relationship that connect h1[n], to
%f0[n] and f1[n] to h0[n] Check this.
[h,h1,f0,f1] = wfilters('db5');

a = fliplr(h);
test = a - f0;

b = fliplr(h1);
test1 = b - f1;

%Write a MATLAB function to compute A(z)B(z) using yor representation for
%such matrices (see function section)
%Form Hac and form the polyphase matrices E(z) and R(z)
e00 = h(1:2:end);
e01 = h(2:2:end);
e10 = h1(1:2:end);
e11 = h1(2:2:end);

r00 = f0(1:2:end);
r01 = f1(1:2:end);
r10 = f0(2:2:end);
r11 = f1(2:2:end);

E = {e00, e01; e10, e11};
R = {r00, r01; r10, r11};

P = matmulconv(R,E);
%computation of P(z) = I 
%computation of the max absolute error, thank you to karol for helping me
%figure out some dimension issues I was having in both cases
set = [zeros(1,4),1,zeros(1,4)];
I = {set, [zeros(1,9)]; [zeros(1,9)],set};
error = zeros(2);

for i = 1:2
    for j = 1:2
        errortemp = cell2mat(P(i,j)) - cell2mat(I(i,j));
        error(i,j) = sum(abs(errortemp));
    end
end

maxabserrorP = max(max(error));

Hac = {h,h1; f0, f1};
Hacpara = {fliplr(cell2mat(Hac(1,1))), fliplr(cell2mat(Hac(2,1))); fliplr(cell2mat(Hac(1,2))), fliplr(cell2mat((Hac(2,2))))};
%checking that Hac is paraunitary
hope = matmulconv(Hacpara,Hac);

set2 = [zeros(1,9),2,zeros(1,9)];
II = {set2, [zeros(1,19)]; [zeros(1,19)],set2};
errorH = zeros(2);

for i = 1:2
    for j = 1:2
        errortempH = cell2mat(hope(i,j)) - cell2mat(II(i,j));
        errorH(i,j) = sum(abs(errortempH));
    end
end

maxabserrorH = max(max(errorH));
%% Structurally Lossless Realization of 2x2 Lossless system

detE = conv(e00, e11) - conv(e01,e10); %we see the delay is four! with beta is 1

%Now we wish to implement the algorithm; after some algebra we find that
%ThetaN = arctan(-(e10 +e11)/(e00+e01)) where N=4 then we keep going for
%ThetaN-1...
%initialize arrays for the loop
%thetaN will be our vector where each theta sub i value will be stored in
%the algorithm, theta5 stored in theta 5
thetaN = zeros(4,1);

%Only four iterations because computing H0 is not done with the algorithm
for i=1:4
    
thetaN(i) = atan(-(E{2,1}(1) + E{2,2}(1))/(E{1,1}(1)+E{1,2}(1)));
% compute givens rotation
R = [cos(thetaN(i)), sin(thetaN(i)); -sin(thetaN(i)), cos(thetaN(i))];
% tranpose givens rotation
R_t = transpose(R);
% convert to 2x2 cell array
R_t = mat2cell(R_t,[1,1],[1,1]);
Gc = matmulconv(R_t,E);
%makes them into matrices to remove the last element
G11 = cell2mat(Gc(1,1));
G12 = cell2mat(Gc(1,2));
G21 = cell2mat(Gc(2,1));
G22 = cell2mat(Gc(2,2));

G21 = G21(2:end);
G22 = G22(2:end);
% recombine to form E matrix of reduced order
E = {G11, G12;G21,G22};
end
% to determine sign, we look at R 
if R(1,1) == R(2,2)
    Sign = 1;
else
    Sign = -1;
end
% to determine alpha, we look at the norm of the columns of R, since we
% have a product of unitary matrices, the resulting matrix will be unitary
% so its columns will be orthonormal, so we check the columns and expect
% them to have norm 1, if they have anything other than norm 1, then we
% know that will be alpha
norm_col1 = norm(R(:,1));
norm_col2 = norm(R(:,2));
% thus, alpha = 1
theta5 = asin(E{1,2}(1))

%% Functions 
%Function for 1(c) that generates the H matrix populated with appropriate
%indices n

function H = Hgen(N)
    n = 0:2*N-1;
    H = 2*n'-n;
    i = find(H>(2*N-1));
    j = find(H<0);
    H(i) = nan;
    H(j) = nan;
end

function X = matmulconv(A,B)
N = size(A);
[X{1:N,1:N}] = deal(zeros(1));
%i is the row, j is the column, k is the entry
    for i = 1:2
        for j = 1:2
            %create an empty vector to store the result, thisll help us
            %with the different sized matrices
            cij = [];
            for k = 1:2
                %we want to make these variables the ones we are convolving
                %in the moment
                temp1 =  A{i,k}
                temp2 =  B{k,j}
                Z = conv(temp1,temp2);
                
                c = size(cij,2);
                z = size(Z,2);
                %comparing sizes, if theyre the same, add Z to cij, else 0
                %pad appropriate one
                if c > z
                    %we are going to have to zero pad Z to make them the
                    %same size
                    cij = cij + [Z, zeros(1,c-z)];
                elseif c < z
                    cij = Z + [cij, zeros(1,z-c)];
                else
                    cij = cij + Z;
                end    
              X(i,j) = {cij}
            end
        end
    end

end


