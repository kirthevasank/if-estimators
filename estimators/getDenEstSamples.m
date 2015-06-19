function [Xden, Xest] = getDenEstSamples(X, numPartitions, currIdx)
% X is the entire dataset. numPartitions is the total number of partitions and
% currIdx is the current partition Idx

  if numPartitions == 1
    Xden = X;
    Xest = X;
  else
    n = size(X, 1);
    startIdx = round( (currIdx-1)*n/numPartitions + 1);
    endIdx = round( currIdx*n/numPartitions );
    estIdxs = startIdx:endIdx;
    denIdxs = [1:(startIdx-1), (endIdx+1):n];
    Xden = X(denIdxs, :);
    Xest = X(estIdxs, :);
  end

end

