function outEEG = WPLI(inEEG)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage:
%
%   >>  outEEG = WPLI(inEEG);
%
% Input(s):
%           inEEG   - input EEG dataset
%    
% Output(s):
%           outEEG  - output EEG dataset
%
% Info:
%           Computes the phase lag index proposed by Martin Vinck et al., 2011 for 
%           each possible pair of EEG channels and for every band as well.
%
% Mathematical background:
%           According to Martin Vinck et al 2011 WPLI is defined in a complete 
%           mathematical formula as:
%
%                   PLI = mean(|sign([sin(Df(t))])|)        (参考)
%                 WPLI = |mean(|[sin(Df(t))]|*sign([sin(Df(t))]))|/mean(|[sin(Df(t))]|)        (1)
%
%           where:
%
%                   Df(t) = phase1(t)-phase2(t)), is the phase difference
%                   of the two signals at time t, 
%
%                   phase1(t) is the phase of the 1st signal and is equal
%                   to arctan(x1_H(t)/x1(t)) where x1_H(t) is the Hilbert
%                   transformed version of signal x1(t),
%
%                   phase2(t) is the phase of the 2nd signal and is equal
%                   to arctan(x2_H(t)/x2(t)) where x2_H(t) is the Hilbert
%                   transformed version of signal x2(t).
%
% Fundamental basis:
%               WPLI is a weighted PLI
%               The PLI ranges between 0 and 1. 
%               A PLI of zero indicates either no coupling or coupling with 
%               a phase difference centered around 0 mod p. A PLI of 1 
%               indicates perfect phase locking at a value of Df different 
%               from 0 mod p. The stronger this nonzero phase locking is, 
%               the larger PLI will be.
%           
% Important notes:
%               the sinus in (1) is used in order to convert the phases on 
%               the interval (0, 2pi] instead of (-pi, pi],
%
%               the mean value is used as an alternative way in order to
%               highlight the fact that when the distibution of Df is 
%               symmetric then the two signals originate from the same
%               source. The same happens when the mean Df is equal to or 
%               centered around 0 mod pi.
%
% Reference(s): Stam C. J., Nolte G., and Daffertshofer A (2007). Phase lag 
% index: assessment of functional connectivity from multi channel EEG and 
% MEG with diminished bias from common sources. Hum. Brain Mapp. 28(3),
% 1178-1193.            
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mf = size(inEEG.freband, 1);
[m, ~, o] = size(inEEG.data);
outEEG = [];

disp('>> FCLAB: PLI is being computed...');

for i = 1:mf
    testEEG = inEEG;
    temp_adj = zeros(m, m);
    freq_range = inEEG.freband{i,1};
    testEEG = eegfiltL(inEEG.data,inEEG.srate,freq_range(1),freq_range(2));
    
    if (o > 1)
        X = mean(testEEG, 3); %events data
    else
        X = testEEG;
    end
        
    for j = 1:m-1
        for k = j+1:m
            hilbert1 = hilbert(X(j, :));
            hilbert2 = hilbert(X(k, :));
            df = angle(hilbert1) - angle(hilbert2);
% % %             temp_adj(j, k) = abs(mean((sign(sin(df)))));
            temp_adj(j, k) = abs(mean(abs(sin(df)).*sign(sin(df))))./mean(abs(sin(df)));
        end
    end
    temp_adj = temp_adj + triu(temp_adj)';
    outEEG(i,:,:) = temp_adj;
% % % % %     eval(['outEEG.PLI.' strrep(inEEG.FC.parameters.bands{i,2},' ','_') '.adj_matrix = temp_adj;']);
end
