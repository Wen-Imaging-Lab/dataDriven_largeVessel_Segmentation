function correlation_coefficient = corrcoefNew(matrixA, matrixB)
    % Check if the matrices have the same dimensions
    if ~isequal(size(matrixA,2), size(matrixB,2))
        error('Input matrices must have the same dimensions.');
    end
    % Calculate the coefficient
    coef = 1/(length(matrixA(1,:))-1);
    % Calculate mean of each row for matrixA and matrixB
    meanA = mean(matrixA, 2);
    meanB = mean(matrixB, 2);
    
    % Calculate the deviations from the means
    deviationA = matrixA - meanA;
    deviationB = matrixB - meanB;
    
    % Calculate the standard deviation of each rows
    stdA = std(matrixA,0,2);
    stdB = std(matrixB,0,2);
    % Calculate the sum of 
    product = sum((deviationA./stdA).*(deviationB./stdB),2);
    % Calculate the correlation coefficient for each row
    correlation_coefficient = product .* coef;
end
