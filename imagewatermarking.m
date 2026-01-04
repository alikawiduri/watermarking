function imagewatermarking(alfa)
% SWT-DST-DCT-SVD-SS + CPT

clc; close all;
rng(1); % supaya PN tidak random

%%%%%%%%%%%%%%%
%%%%%Jangan diubah-ubah
%%%%%%%%%%%%%%%
M=16;N=16;
Ldek=1;
mode=1;                  % SS
sub_eksis=[1 1 1 1 1];    % SWT-DST-DCT-SVD-CPT
resolusihost=256;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%BOLEH DIUBAH%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub=2;
blok=8;
alfass=180;
alfa=100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L=blok^2;
folderhost=[pwd '/host_image/'];
folderwatermark=[pwd '/watermark/'];

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% PREPARATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
x=imread([folderhost 'cameraman.jpg']);
x=im2gray(x);
x1=double(imresize(x,[resolusihost resolusihost]));

logo=imread([folderwatermark 'logo tel-u.png'],'png');
logor=imresize(logo,[size(x1,1)/blok size(x1,2)/blok]);
logogray=rgb2gray(logor);
logobw=~imbinarize(logogray,0.0000001);
logobw1d=logobw(:);

host1d=reshape(x1,[],1);
xw=zeros(size(host1d));

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% PN CODE %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
pn0=2*randi([0 1],blok,blok)-1;
pn1=2*randi([0 1],blok,blok)-1;

[U0,~,V0]=svd(pn0);
[U1,~,V1]=svd(pn1);

r_ref = zeros(length(logobw1d),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% EMBEDDING %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(logobw1d)

    seg=alfa*host1d(L*i-L+1:L*i);

    % SWT
    Sw=swt(seg,Ldek,'haar');

    % DST
    sdst=dst(Sw(sub,:));

    % DCT
    Swdct=dct(sdst);

    % 1D → 2D
    Swdct2=reshape(Swdct,[blok blok]);

    % SVD
    [U,S,V]=svd(Swdct2);

    % ===== CPT =====
    [r,theta_rad,~]=ctp(S(1,1),S(2,2));
    r_ref(i) = r;


    % ===== SS + CPT =====
    if logobw1d(i)==1
        r_emb=r+alfass;
    else
        r_emb=r-alfass;
    end

    % inverse CPT
    [S(1,1),S(2,2)]=ptc(r_emb,theta_rad);

    % ISVD
    Sisvd=U*S*V';

    % 2D → 1D
    Sisvd1=reshape(Sisvd,[blok^2 1]);

    % IDCT
    s_idct=idct(Sisvd1);

    % IDST
    Sw(sub,:)=idst(s_idct);

    % ISWT
    st=iswt(Sw,'haar');

    xw(L*i-L+1:L*i)=st/alfa;
end

xwb=(xw+1)*255/2;
xwb=min(max(xwb,0),255);
xwr=reshape(uint8(xwb),size(x1));

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% ATTACK %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
attacked_img=xwr;   % TANPA SERANGAN dulu

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% EXTRACTION %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
xwn=reshape(double(attacked_img),[],1);
wt=zeros(length(logobw1d),1);

for i=1:length(logobw1d)

    seg=alfa*(2*xwn(L*i-L+1:L*i)/255-1);

    % SWT
    Sw=swt(seg,Ldek,'haar');

    % DST
    sdst=dst(Sw(sub,:));

    % DCT
    Swdct=dct(sdst);

    % 1D → 2D
    Swdct2=reshape(Swdct,[blok blok]);

    % SVD
    [~,S,~]=svd(Swdct2);

    % ===== CPT extraction =====
    [r_hat,~,~]=ctp(S(1,1),S(2,2));

    % ===== SS decision =====
    if r_hat > r_ref(i)
        wt(i)=1;
    else
        wt(i)=0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% EVALUASI %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
ber=sum(wt~=logobw1d)/length(logobw1d);
mse=mean((x1(:)-double(xwr(:))).^2);
psnr=10*log10(255^2/mse);

hasil.ber=ber;
hasil.psnr=psnr;
hasil.payload=1/L;
hasil

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% DISPLAY %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
Nw = numel(logor);
logobw2d = reshape(wt(1:Nw),size(logor));
figure(1),clf
subplot(121),imshow(logobw),title('Watermark Asli')
subplot(122),imshow(logobw2d),title('Hasil Ekstraksi')
