%% Reconstruction
function myJPEG(path)
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
dict = load('dict.mat','dict');
dict = getfield(dict,'dict');
[filepath,name,ext] = fileparts(path);
if strcmp(ext,'.myJPEG')==0,
    disp("only .myJPEG is allowed")
end
imsize = [128 128];
fid = fopen('out.myJPEG');
binaryComp = fread(fid);
fclose(fid);
eob = 100000;
comp = bi2de(binaryComp);
r = huffmandeco(comp,dict);
Recon = zeros(128);
i = 1;
j = 1;
DCT_vector = zeros(64,1);
k=0;
for jj = 1:length(r),
    if r(jj)~=eob,
        k=k+1;
        DCT_vector(order(k)) = r(jj);
    else %block finishes
        k = 0;
        DCT_block = reshape(DCT_vector,[8 8]);
        block = idct2(DCT_block.*Norm_Mat);
        if j<imsize(2)/8,
            Recon(8*(i-1)+1:8*i,8*(j-1)+1:8*j) = block;
            j = j+1;
        else,
            Recon(8*(i-1)+1:8*i,8*(j-1)+1:8*j) = block;
            i = i+1;
            j = 1;
        end
        DCT_vector = zeros(64,1);
    end
end
Recon = Recon+128;
Recon = (Recon-min(min(Recon)))/max(max(Recon));
imshow(Recon);
title('myJPEG Reconstruction')