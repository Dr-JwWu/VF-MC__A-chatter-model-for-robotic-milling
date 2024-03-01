function [DA,DP] = TimeDelay(omiga,To)
%{
Function:
    Calculate the influence coefficients of the time delay term

Input:
    omiga:  Angular frequency
    To:     Tooth passing period

Output:
    DA: Influence of the time delay term on the amplitude
    DP: Influence of the time delay term on the phase
%}

%%
TD = 1 - cos(omiga*To) + 1j*sin(omiga*To);
DA = norm(TD);
DP = angle(TD);
end