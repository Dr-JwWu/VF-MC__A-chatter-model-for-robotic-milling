function [Y,f,P1] = FourTran(X,Fs)

%{
Function:
    Fast Fourier Transform (FFT)

Input:
    X:	Column vector of signals
    Fs: Sampling frequency

Output:
    angle(Y):	Column vector of phases
    f:          Column vector of frequencies
    P1:         Column vector of amplitudes 

Remarks:
    The frequency band is 0 ~ Fs/2 and the resolution is Fs/L.
%}

%%
if mod(length(X),2) ~= 0 % Change the length L of X to an even number.
    X(end,:) = [];
end
L = length(X); 
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:L/2)'/L; 
end