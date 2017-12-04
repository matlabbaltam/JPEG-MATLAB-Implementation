%% Jpeg implementation without coding
function jpeg_code(path)
Norm_Mat=[16 11 10 16 24 40 51 61       % Normalization matrix (8 X 8) used to Normalize the DCT Matrix.
          12 12 14 19 26 58 60 55
          14 13 16 24 40 57 69 56
          14 17 22 29 51 87 80 62
          18 22 37 56 68 109 103 77
          24 35 55 64 81 104 113 92
          49 64 78 87 103 121 120 101
          72 92 95 98 112 100 103 99];
order = [1 9 2 3 10 17 25 18 11 4 5 12 19 26 33  ...
         41 34 27 20 13 6 7 14 21 28 35 42 49 57 50 ...
         43 36 29 22 15 8 16 23 30 37 44 51 58 59 52 ...
         45 38 31 24 32 39 46 53 60 61 54 47 40 48 55 ...
         62 63 56 64];
image_data = imread(path);
imsize = size(image_data);
image_temp = image_data;
image_temp = double(image_temp)-128;
num_pad = [0,0];

%%zeor_padding
if mod(imsize(1),8)~=0,
    num_pad(1) = ceil(imsize(1)/8)*8-imsize(1);
end
if mod(imsize(2),8)~=0,
    num_pad(2) = ceil(imsize(2)/8)*8-imsize(2);
end
image_temp = padarray(image_temp,num_pad);
imsize_pad = size(image_temp);
eob = 100000;  %create eob symbol
r = zeros(numel(image_temp)+imsize_pad(1)*imsize_pad(2),1);

%% DCT of each sub-block and keep N coefficients
count = 0;
rec_pad = zeros(imsize_pad);
for i = 1:ceil(imsize(1)/8),
    for j = 1:ceil(imsize(2)/8),
        sub_block = image_temp(8*(i-1)+1:8*i,8*(j-1)+1:8*j);
        block_DCT = dct2(sub_block);
        block_quan = round(block_DCT./Norm_Mat);
        vector_quan = block_quan(:);
        vector_reord = vector_quan(order,:);
        ii = find(vector_reord, 1, 'last');   % find last non-zero element
        if isempty(ii),                   % check if there are no non-zero values
            ii = 0; 
        end 
        p = count + 1;
        q = p + ii;
        r(p:q)  = [vector_reord(1:ii); eob];     % truncate trailing zeros, add eob
        count = count + ii + 1;          % and add to output vector
    end
end
r((count + 1):end) = [];    
symbols = unique(r);
counts = hist(r,symbols);
probability = counts/sum(counts);
dict = huffmandict(symbols,probability);
comp = huffmanenco(r,dict);
binaryComp = de2bi(comp);
% encodedLen = numel(binaryComp);

fid = fopen('out.myJPEG','wt');
fwrite(fid,binaryComp);
fclose(fid);

% save('out.myJPEG','binaryComp')
save('dict.mat','dict');
