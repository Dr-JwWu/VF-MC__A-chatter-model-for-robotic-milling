function [X,Fi] = FRF(wn,damp,K,wfre)

%{
Function:
    Calculate the amplitude and phase of the frequency response function (FRF)

Input:
    wn:     Modal angular frequency
    damp:   Damping ratio
    K:      Modal stiffness
    wfre:   Angular frequency of excitation

Output:
    X:  Amplitude of the FRF (i.e., dynamic compliance)
    Fi: Phase of the FRF
%}

%%
X = 1./ K ./ sqrt( (1-(wfre./wn).^2).^2 + (2.*damp.*wfre./wn).^2 );
Fi = -atan( (2.*damp.*wfre./wn) ./ (1-(wfre./wn).^2) ); 

Fi(Fi>0) = Fi(Fi>0) - pi; % Adjust the range of Fi to (-pi, 0]

end