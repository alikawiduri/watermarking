function x1=stirmark_imagebenchmark(watermarked_img,jenis)

% jenisserangan=[ {'no_attack'},{'jpeg_compression'}, {'rotation_attack'}, {'horizontal_flip'}, {'random_cropping'}, {'uniform_scaling'}, ...
%     {'non_uniform_scaling'}, {'random_row_col_removal'}, {'combined_transform'}, {'random_geometric_distortion'}, ...
%     {'geometric_distortion_with_jpeg'}, {'lpf_attack'}, {'sharpen_attack'}, {'histogram_modification'}, ...
%     {'gamma_correction'}, {'color_quantization'}, {'noise_addition'}, {'statistical_averaging'}, ...
%     {'overmarking'}, {'oracle_attack'}];
% Initialize parameters
attacks = {@no_attack,@jpeg_compression, @rotation_attack, @horizontal_flip, @random_cropping, @uniform_scaling, ...
    @non_uniform_scaling, @random_row_col_removal, @combined_transform, @random_geometric_distortion, ...
    @geometric_distortion_with_jpeg, @lpf_attack, @sharpen_attack, @histogram_modification, ...
    @gamma_correction, @color_quantization, @noise_addition, @statistical_averaging, ...
    @overmarking, @oracle_attack};
attack_names = {'Tanpa Serangan','JPEG Compression', 'Rotation', 'Horizontal Flip', 'Random Cropping', 'Uniform Scaling', ...
    'Non-Uniform Scaling', 'Row/Col Removal', 'Combined Transform', 'Random Geometric Distortion', ...
    'Geometric Distortion with JPEG', 'Low Pass Filter', 'Sharpening', 'Histogram Modification', ...
    'Gamma Correction', 'Color Quantization', 'Noise Addition', 'Statistical Averaging', ...
    'Overmarking', 'Oracle Attack'};

% attacked_img = attacks{[char(jenisserangan(jenis(1)))]}(watermarked_img,jenis(2));

% =====================================================
% FIX FINAL: Handle numeric & cell attack parameters
% =====================================================
if iscell(jenis)
    jenis1 = jenis{1};   % attack index (numeric)
    jenis2 = jenis{2};   % attack parameter (numeric / vector)
else
    jenis1 = jenis(1);
    jenis2 = jenis(2:end);
end

attacked_img = attacks{jenis1}(watermarked_img, jenis2);

x1=double(attacked_img);
% x2=reshape(x1,[size(x1,1)*size(x1,2),1]);

% if jenis(1)==1
%     serangan=@jpeg_compression;
% elseif jenis(1)==2
%     serangan=@rotation_attack;
% elseif jenis(1)==2
%     serangan=@horizontal_flip;
% elseif jenis(1)==2
%     serangan=@random_cropping;
% elseif jenis(1)==2
%     serangan=@jpeg_compression;
% elseif jenis(1)==2
%     serangan=@jpeg_compression;
% elseif jenis(1)==2
%     serangan=@jpeg_compression;
% elseif jenis(1)==2
%     serangan=@jpeg_compression;
% elseif jenis(1)==2
%     serangan=@jpeg_compression;

% Loop through each attack and evaluate
% for i = 1:length(attacks)
%     fprintf('Applying %s attack...\n', attack_names{i});
%     attacked_img = attacks{i}(watermarked_img);
%
%     % Calculate Bit Error Rate (BER)
%     ber = calculate_ber(original_img, attacked_img);
%     fprintf('%s Attack BER: %.4f\n', attack_names{i}, ber);
%
%     % Show the attacked image (optional)
%     figure;
%     imshow(attacked_img);
%     title([attack_names{i} ' Attack']);
% end
end


%% Attack functions 1
function attacked_img=no_attack(img,x)
attacked_img=img;
end
%2
function attacked_img = jpeg_compression(img,qf)
imwrite(uint8(img), 'jpeg_temp.jpg', 'jpg', 'Quality',qf);
% imwrite(img, 'jpeg_temp.jpg', 'JPEG 2000', 'CompressionRatio',1,'Mode','lossy');
attacked_img = imread('jpeg_temp.jpg');
% psnr(double(attacked_img),img)
% figure(1),clf,
% subplot(121),imshow(uint8(img)),title('Asli')
% subplot(122),imshow(uint8(attacked_img)),title('Attacked')

end
%3
function attacked_img = rotation_attack(img,angle)
% angle = 5; % Rotation angle in degrees
attacked_img = imrotate(img, angle, 'bilinear', 'crop');
end
%4
function attacked_img = horizontal_flip(img,col)
attacked_img = flip(img, col);
end
%5
function attacked_img = random_cropping(img, persen)

% --- pastikan ukuran target hanya [H W]
target_size = size(img);
target_size = target_size(1:2);

% --- hitung ukuran crop
crop_size = floor(target_size * persen/100);

% --- safety
crop_size = max(crop_size, [1 1]);

x = randi([1, target_size(1) - crop_size(1) + 1]);
y = randi([1, target_size(2) - crop_size(2) + 1]);

if ndims(img) == 3
    cropped = img(x:x+crop_size(1)-1, y:y+crop_size(2)-1, :);
else
    cropped = img(x:x+crop_size(1)-1, y:y+crop_size(2)-1);
end

% --- resize KEMBALI ke ukuran asli (H W)
attacked_img = imresize(cropped, target_size);

end

%6
function attacked_img = uniform_scaling(img,scale)
% scale = 0.9; % Scaling factor
attacked_img = imresize(img, scale);
attacked_img = imresize(attacked_img, size(img));
end

%7
function attacked_img = non_uniform_scaling(img,scale)
% scale_x = 0.8; % Horizontal scaling
% scale_y = 1.2; % Vertical scaling
% Kalau masih cell, ambil isinya
if iscell(scale)
    scale = scale{1};
end

% Kalau scalar → jadikan uniform scaling
if numel(scale) == 1
    scale_x = double(scale);
    scale_y = double(scale);
else
    scale_x = double(scale(1));
    scale_y = double(scale(2));
end

% Safety clamp (hindari ukuran nol / negatif)
scale_x = max(scale_x, 0.1);
scale_y = max(scale_y, 0.1);

attacked_img = imresize(img, ...
    [round(size(img,1)*scale_y), round(size(img,2)*scale_x)]);

attacked_img = imresize(attacked_img, size(img));
end

%8
% function attacked_img = random_row_col_removal(img,num)
% img(randi(size(img, 1), 1), :) = 0; % Remove random row
% img(:, randi(size(img, 2), 1)) = 0; % Remove random column
% attacked_img = img;
% end
function attacked_img = random_row_col_removal(img, num)
% Fungsi untuk menghapus num baris dan num kolom secara acak dari gambar

% Dapatkan ukuran gambar
[rows, cols] = size(img);

% Pastikan num tidak melebihi jumlah baris/kolom
num_rows = min(num, rows);
num_cols = min(num, cols);

% Pilih baris dan kolom acak untuk dihapus
row_indices = randperm(rows, num_rows);
col_indices = randperm(cols, num_cols);

% Hapus baris
img(row_indices, :) = 0;

% Hapus kolom
img(:, col_indices) = 0;

% Output gambar yang telah diserang
attacked_img = img;
end

%9
function attacked_img = combined_transform(img, val)
% =====================================================
% FINAL FIX: robust combined transform
% =====================================================

% Kalau masih cell, ambil isinya
if iscell(val)
    val = val{1};
end

% Jika scalar → default angle + scale
if numel(val) == 1
    sudut  = double(val) * 10;   % rotasi kecil
    ukuran = double(val);        % scaling
else
    sudut  = double(val(1));
    ukuran = double(val(2));
end

% Safety clamp
ukuran = max(ukuran, 0.1);

img2 = imrotate(img, sudut, 'bilinear', 'crop');
attacked_img = imresize(img2, ukuran);
attacked_img = imresize(attacked_img, size(img));
end

%10
function attacked_img = random_geometric_distortion(img,x)
tform = affine2d([1 0.1 0; 0.1 1 0; 0 0 1]);
output_view = imref2d(size(img));
attacked_img = imwarp(img, tform, 'OutputView', output_view);
end

%11
function attacked_img = geometric_distortion_with_jpeg(img,qf)
img = random_geometric_distortion(img);
attacked_img = jpeg_compression(img,qf);
end

%12
function attacked_img = lpf_attack(img,x)
h = fspecial('gaussian', [x, x], 1);
attacked_img = imfilter(img, h, 'replicate');
end

%13
function attacked_img = sharpen_attack(img,x)
attacked_img = imsharpen(img);
end

%14
function attacked_img = histogram_modification(img,x)
attacked_img = histeq(img);
end

%15
function attacked_img = gamma_correction(img,gamma)
% gamma = 0.8;
attacked_img = imadjust(img, [], [], gamma);
end

%16
function attacked_img = color_quantization(img,levels)
% levels = 16;
attacked_img = round(img / levels) * levels;
end

%17
function attacked_img = noise_addition(img,noise_intensity)
% noise_intensity = 0.02;
attacked_img = imnoise(img, 'gaussian', 0, noise_intensity);
end

%18
function attacked_img = statistical_averaging(img,var)
img2 = imnoise(img, 'gaussian', 0, var);
attacked_img = (double(img) + double(img2)) / 2;
end

%19
function attacked_img = overmarking(img,x)
watermark = randi([0, 1], size(img));
attacked_img = bitxor(img, watermark);
end

%20
function attacked_img = oracle_attack(img,x)
attacked_img = img;
% threshold = 0.1;
threshold = x;
while watermark_detected(attacked_img)
    attacked_img = attacked_img + threshold * randn(size(attacked_img));
end
end

%% Dummy watermark detection function
function detected = watermark_detected(img)
% Placeholder for actual watermark detection logic
detected = rand > 0.5;
end


%% BER Calculation
% function ber = calculate_ber(original_img, attacked_img)
%     original_bits = imbinarize(original_img);
%     attacked_bits = imbinarize(attacked_img);
%
%     % Ensure both images are of the same size
%     if size(original_bits) ~= size(attacked_bits)
%         error('Original and attacked images must have the same dimensions.');
%     end
%
%     % Calculate Bit Error Rate
%     error_bits = sum(original_bits(:) ~= attacked_bits(:));
%     total_bits = numel(original_bits);
%     ber = error_bits / total_bits;
% end
