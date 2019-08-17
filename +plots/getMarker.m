function [ marker ] = getMarker(index, markers)
%GETMARKER returns a marker for a index
    marker = markers{mod(index, length(markers))+1};
end

