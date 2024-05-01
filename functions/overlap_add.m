function y = overlap_add(b, data)
% loading parameters
nfft = b.nfft;
ndata = b.ndata;
L = b.window;
G = b.overlap;

B = fftshift(b.filter);

if length(B) == 1
    B = B(:); % make sure that B is a column
end
if size(B,2) == 1
    B = B(:, ones(1, size(data,2))); % replicate the column B
end
if size(data,2) == 1
    data = data(:, ones(1, size(b,2))); % replicate the column data
end

y = zeros(size(data));
istart = 1;
while istart <= ndata
    iend = min(istart+L-1, ndata);
    if (iend - istart) <= 0
        break
    else
        X = fft([zeros(1,G)'; data(istart:iend,:)], nfft);
    end
    Y = ifft(X.*B);
    yend = min(ndata, istart+nfft-1);
    y(istart:yend,:) = y(istart:yend,:) + Y(1:(yend-istart+1),:);
    istart = istart + L;
end

end