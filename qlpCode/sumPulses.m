function result = sumPulses(data, numPulses)
    [numRows, numCols] = size(data);
    
    % Determine the number of loops and remaining pulses
    numLoops = floor(numRows / numPulses);
    remPulses = mod(numRows, numPulses);
    
    % Create an array to store the summed results
    result = zeros(numLoops + (remPulses > 0), numCols);
    
    % Sum the rows using vectorized operations
    for i = 1:numLoops
        rowsToSum = (i - 1) * numPulses + 1 : i * numPulses;
        result(i, :) = sum(data(rowsToSum, :));
    end
    
    % If there are remaining rows to sum
    if remPulses > 0
        rowsToSum = numLoops * numPulses + 1 : numRows;
        result(end, :) = sum(data(rowsToSum, :));
    end
end