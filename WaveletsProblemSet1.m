%% Michael Lendino Wavelets PSET 1 unfinished but its not due yet so its okay
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

%Compute the four filters for db5 and plot the results.
[H0,H1,F0,F1] = wfilters('db5'); 
subplot(2,2,1)
stem(H0)
title('H0 Decomposition Lowpass Filter')
subplot(2,2,2)
stem(H1)
title('H1 Decomposition Highpass Filter')
subplot(2,2,3)
stem(F0)
title('F0 Reconstruction Lowpass Filter')
subplot(2,2,4)
stem(F1)
title('F1 Reconstruction Highpass Filter')
xlabel(['The four filters for ','db5'])
%Confirming time domain relationship
a = fliplr(H0);
test = a - F0;

b = fliplr(H1);
test1 = b - F1;

%Write a MATLAB function to compute A(z)B(z) using yor representation for
%such matrices
 
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


