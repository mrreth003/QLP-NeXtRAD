function plotSPEC(data, windowSize, overlap, prf, isFirst)
persistent himage;
% Calculate num windows required
numWindows = floor((length(data) - overlap * windowSize) / (windowSize - overlap * windowSize));

% Initialise spec matrix
spec = zeros(windowSize * 11, numWindows);

% Axes
tAxis = (0:size(data, 1) - 1) / (prf * 1e3);
fAxis = linspace(- (prf * 1e3) / 2, (prf * 1e3) / 2, size(data, 1));

% Compute spec
for i = 1:numWindows
    % Extract window
    startIdx = (i - 1) * floor(windowSize - overlap * windowSize) + 1;
    endIdx = startIdx + windowSize - 1;
    rangeWindow = data(startIdx:endIdx);


    % Apply hann window
    win = zeros(windowSize, 1);
    for n = 1:windowSize
        win(n) = 0.5 * (1 - cos(2 * pi * n / (windowSize - 1)));
    end
    rangeWindow = rangeWindow .* win;

    % Zero pad window
    rangeWindow = [rangeWindow.', zeros(length(rangeWindow) * 10, 1).']';

    % Store in spec matrix
    spec(:, i) = fftshift(fft(rangeWindow));
end

% PLOT
if isFirst
    figure(3);
    himage = imagesc(tAxis, fAxis, 20 * log10(abs(spec)));
    axis xy; colormap('jet'); colorbar;
    xlabel('Time (s)');
    ylabel('Doppler Frequency (Hz)');
    title('NeXtRAD Spectrogram Quick-Look');
else
    figure(3);
    set(himage, 'CData', 20 * log10(abs(spec)));

    % Adjust the axis limits
    axis([0, max(tAxis), 0, max(fAxis)]);

    % Automatically set the color axis limits based on the new data
    clim auto;

    drawnow;
end
end