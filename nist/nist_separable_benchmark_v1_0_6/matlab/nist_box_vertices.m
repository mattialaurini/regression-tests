function vertices = nist_box_vertices(lowerBound,upperBound)
%NIST_BOX_VERTICES Return the 2^d vertices of a box, one vertex per column.
lowerBound = double(lowerBound(:));
upperBound = double(upperBound(:));
if numel(lowerBound) ~= numel(upperBound) || any(lowerBound>=upperBound)
    error("Invalid box bounds.");
end
d = numel(lowerBound);
bits = dec2bin(0:2^d-1,d)-'0';
vertices = lowerBound + (upperBound-lowerBound).*bits.';
end
