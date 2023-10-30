function qlp()
%% USER SETTINGS
setPulses = 1000;
useParamRefSig = false;

%CHOOSE VISIBLE PLOT
isMFO = false;
isRTI = true;
isRDM = false;
isSPEC = false;
isCFAR = false;
isCFARO = false;
isRange = false;

%CFAR
guardCellsCFAR = 200; % num each side of CUT
refCellsCFAR = 20; % num each side of CUT
pfa = 1e-6; % probability of false alarm
%SPEC
windowSizeSPEC = 20;
overlapSPEC = 0.95;
rangeBinSPEC = 1200;

%DATA FILES
setupParams = "e10_10_16_1800_40_P1_1_130000_S0_1_2047_node3.m";
dataNode = "e10_10_16_1800_40_P1_1_130000_S0_1_2047_node3.bin";
refSig = "..\reference_signals\refSigN3N2Pl1000.txt";

%SUM PULSES
isSum = false;


%% LOAD PARAMS
% Load parameters from .m file
run(setupParams);

%% LOAD DATA

% Define the naming pattern you want to match
patternA = 'n\dnumOfSamplesPerPulse';
patternB = 'n\dnumOfPulses';
patternC = 'n\drefsig';
patternD = 'n\dBandwidth';
patternE = 'n\dPRF';
patternF = 'n\dCapturetime';


% Get a list of variables in the workspace
workspaceVariables = whos;

% Loop through the variables and find the ones that match the pattern
for i = 1:numel(workspaceVariables)
    variableName = workspaceVariables(i).name;
    if ~isempty(regexp(variableName, patternA, 'once'))
        % Use 'eval' to assign the variable to a new variable in your workspace
        numSamples = eval(variableName);
    end
    if ~isempty(regexp(variableName, patternB, 'once'))
        % Use 'eval' to assign the variable to a new variable in your workspace
        numPulses = eval(variableName);
    end
    if ~isempty(regexp(variableName, patternC, 'once'))
        ref = eval(variableName);
    end
    if ~isempty(regexp(variableName, patternD, 'once'))
        % Use 'eval' to assign the variable to a new variable in your workspace
        bw = eval(variableName);
    end
    if ~isempty(regexp(variableName, patternE, 'once'))
        % Use 'eval' to assign the variable to a new variable in your workspace
        prf = eval(variableName);
    end
    if ~isempty(regexp(variableName, patternF, 'once'))
        % Use 'eval' to assign the variable to a new variable in your workspace
        ct = eval(variableName);
    end
end

numPulses = setPulses;
chunk = 0.01 * numPulses;
if isCFAR
    if ((refCellsCFAR + guardCellsCFAR + 1)* 2 > chunk)
        chunk = (refCellsCFAR + guardCellsCFAR + 1) * 2;
    end
end
pulsesPerLook = 0.2 * numPulses;
numPulses2Sum = 0.005 * numPulses;


% DATA
% Load data from .bin file
fnode = fopen(dataNode, "r");
binData = fread(fnode, numPulses * numSamples, "uint16"); % Assuming data stored as unsigned 16-bit integer
fclose(fnode);

% data offset
% Calculate the mean
meanValue = mean(binData);

% Subtract the mean from each element in binData
binData = binData - meanValue;

% convert to complex and reshape
complexData = hilbert(binData);

shapedData = reshape(complexData, numSamples, []).';

% REFERENCE

if ~useParamRefSig
% Load reference signal from .txt file
binRef = readmatrix(refSig);

% normalise & window
refSignal = binRef;
for i = 1:length(binRef)
    refSignal(i) = (binRef(i))/ max(binRef) * (0.5 * (1 - cos((2*pi*i) / (length(binRef) - 1))));
end

% convert to complex
refSignal = hilbert(refSignal);

else
    refSignal = ref;
end

% Zero pad ref signal
paddedRefSignal = zeros(1, numSamples) + 1i * 0;
paddedRefSignal(numSamples - length(refSignal) + 1: numSamples) = refSignal;

%% MATCHED FILTER & PLOTTING

% Initialize the result matrix
refSignal = conj(fliplr(paddedRefSignal));
matchedFilterOutput = zeros(size(shapedData));

isFirst = false;
% Perform convolution along rows
for row = 1:numPulses
    matchedFilter = ifft(fft(shapedData(row, :)) .* fft(refSignal)); % convolution
    matchedFilterOutput(row, :) = matchedFilter(1:length(refSignal));

    if ~isSum
        % PLOT RTI
        if (isRTI && mod(row, chunk) == 0)
            if row == chunk
                isFirst = true;
            end
            plotRTI(matchedFilterOutput(row - chunk + 1 : row, :), bw, prf, isFirst);
            isFirst = false;
        end
        % PLOT Range-Doppler Map
        if (isRDM && mod(row, pulsesPerLook) == 0)
            if row == pulsesPerLook
                isFirst = true;
            end
            plotRDM(matchedFilterOutput(row - pulsesPerLook + 1:row, :), bw, isCFAR, guardCellsCFAR, refCellsCFAR, pfa, isCFARO, isRange, isFirst, prf);
            isFirst = false;
        end
        % PLOT Spectrogram
        if (isSPEC && mod(row, pulsesPerLook) == 0)
            if row == pulsesPerLook
                isFirst = true;
            end
            plotSPEC(matchedFilterOutput(row - pulsesPerLook + 1:row , rangeBinSPEC/(299792458 / (2 * bw * 1e6))), windowSizeSPEC, overlapSPEC, prf, isFirst);
            isFirst = false;
        end
    end
end

if isSum
    matchedFilterOutput = sumPulses(matchedFilterOutput, numPulses2Sum);
    numPulses = size(matchedFilterOutput, 1);

    for row = 1:numPulses
        % PLOT RTI
        if (isRTI && mod(row, chunk) == 0)
            if row == chunk
                isFirst = true;
            end
            plotRTI(matchedFilterOutput(row - chunk + 1 : row, :), bw, prf, isFirst);
            isFirst = false;
        end
        % PLOT Range-Doppler Map
        if (isRDM && mod(row, pulsesPerLook) == 0)
            if row == pulsesPerLook
                isFirst = true;
            end
            plotRDM(matchedFilterOutput(row - pulsesPerLook + 1:row, :), bw, isCFAR, guardCellsCFAR, refCellsCFAR, pfa, isCFARO, isRange, isFirst);
            isFirst = false;
        end
        % PLOT Spectrogram
        if (isSPEC && mod(row, pulsesPerLook) == 0)
            if row == pulsesPerLook
                isFirst = true;
            end
            plotSPEC(matchedFilterOutput(row - pulsesPerLook + 1:row , rangeBinSPEC/(299792458 / (2 * bw * 1e6))), windowSizeSPEC, overlapSPEC, prf, isFirst);
            isFirst = false;
        end
    end
end


% Plot Matched Filter Output
if isMFO
    figure();
    plot(20*log10(abs(matchedFilterOutput(5000, :))));
    xlabel('Range Bin Number');
    ylabel('Magnitude [dB]')
    title('Output of Matched Filter');
end
end