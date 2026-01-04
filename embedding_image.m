function [xwr,So,U0,U1,V0,V1,pn0,pn1,wt]=embedding_image(host2d,host1d,logonrz,alfa,alfass,L,sub_eksis,Ldek,subband,nb,jenis,B,mode)
if mode==0 %QIM
    sub=subband;
    logobw1d=logonrz;
    if sub_eksis(1)==2
    blok=B/2;
    else
        blok=B;
    end

    %normalisasi
    host1dn=2*host1d/255-1;

    % host2d=host1d;
    for i=1:length(logobw1d)
        % segmen 1: x(1:L,1)
        % segmen 2: x(L+1:2L,1)
        %Segmentasi
        seg=alfa*host1dn(L*i-L+1:L*i,1);%segmentasi
        if sub_eksis(1)==1 %SWT
            Sw=swt(seg,Ldek,'haar');
        elseif sub_eksis(1)==0
            Sw=seg.';
        elseif sub_eksis(1)==2
            [A,D]=dwt(seg','haar');
            [LL,LH]=dwt(A,'haar');
            [HL,HH]=dwt(D,'haar');
            Sw=[LL;LH;HL;HH];
        end

        if sub_eksis(2)==0
            sdst=Sw(sub,:);
        elseif sub_eksis(2)==1 %DST
            sdst=dst(Sw(sub,:));
        end
        if sub_eksis(3)==0
            Swdct=sdst;
        elseif sub_eksis(3)==1 %DCT
            Swdct=dct(sdst);
        end
        Swdct2=reshape(Swdct,[blok blok]);
        if sub_eksis(4)==1 %SVD
            [U,S,V]=svd(Swdct2);
        elseif sub_eksis(4)==2
            [U,S]=qr(Swdct2);
        end

        %CPT dengan berbagai parameter
        if sub_eksis(5)==0 %Tanpa CPT
            swcpt(i)=S(1);
        elseif sub_eksis(4)==1 && sub_eksis(5)==1 %disisipkan pada r (SVD)
            %SVD-CPT
            [r, theta, theta_deg] = ctp(S(1),S(2));
            swcpt(i)=r;
        elseif sub_eksis(4)==1 && sub_eksis(5)==2 %disisipkan pada theta_deg (SVD)
            [r, theta_rad, theta_deg] = ctp(S(1),S(2));
            swcpt(i)=theta_deg;
        elseif sub_eksis(4)==1 && sub_eksis(5)==3 %disisipkan pada theta_rad (SVD)
            [r, theta_rad, theta_deg] = ctp(S(1),S(2));
            swcpt(i)=theta_rad;
        elseif sub_eksis(4)==2 && sub_eksis(5)==1 %disisipkan pada r (QR)
            [r, theta, phi] = ctp_3d(S(1,1), S(2,2), S(1,2));
            swcpt(i)=r;
        elseif sub_eksis(4)==2 && sub_eksis(5)==2 %disisipkan pada theta_deg (QR)
            [r, theta, phi] = ctp_3d(S(1,1), S(2,2), S(1,2));
            swcpt(i)=rad2deg(theta);
        elseif sub_eksis(4)==2 && sub_eksis(5)==3 %disisipkan pada theta_rad (QR)
            [r, theta, phi] = ctp_3d(S(1,1), S(2,2), S(1,2));
            swcpt(i)=(theta);
        elseif sub_eksis(4)==2 && sub_eksis(5)==4 %disisipkan pada phi_deg (QR)
            [r, theta, phi] = ctp_3d(S(1,1), S(2,2), S(1,2));
            swcpt(i)=rad2deg(phi);
        elseif sub_eksis(4)==2 && sub_eksis(5)==5 %disisipkan pada phi_rad (QR)
            [r, theta, phi] = ctp_3d(S(1,1), S(2,2), S(1,2));
            swcpt(i)=(phi);
        end



        %QIM embedding
        % r_emb(i)=QIM_emb(r,logobw1d(i),nb,jenis);
        % theta_emb(i)=QIM_emb(theta_deg,logobw1d(i),nb,jenis);
        qim_emb(i)=QIM_emb(swcpt(i),logobw1d(i),nb,jenis);
        % wt(i,1)=QIM_ext(qim_emb(i),nb,jenis);

        Sdgab=S;
        if sub_eksis(5)==0 %Tanpa CPT
            Sdgab(1)=qim_emb(i);


        elseif sub_eksis(4)==1 && sub_eksis(5)==1 %disisipkan pada r (SVD)
            %SVD-CPT
            r=qim_emb(i);
            [Semb1(i), Semb2(i)] = ptc(r, theta);
            Sdgab(1)=Semb1(i);
            Sdgab(2)=Semb2(i);

        elseif sub_eksis(4)==1 && sub_eksis(5)==2 %disisipkan pada theta_deg (SVD)
            theta=qim_emb(i);
            [Semb1(i), Semb2(i)] = ptc(r,theta);
            Sdgab(1)=Semb1(i);
            Sdgab(2)=Semb2(i);
        elseif sub_eksis(4)==1 && sub_eksis(5)==3 %disisipkan pada theta_rad (SVD)
            theta=qim_emb(i);
            [Semb1(i), Semb2(i)] = ptc(r,rad2deg(theta));
            Sdgab(1)=Semb1(i);
            Sdgab(2)=Semb2(i);
        elseif sub_eksis(4)==2 && sub_eksis(5)==1 %disisipkan pada r (QR)
            r=qim_emb(i);
            [Semb1(i), Semb2(i),Semb3(i)] = ptc_3d(r, theta, phi);
            Sdgab(1,1)=Semb1(i);
            Sdgab(2,2)=Semb2(i);
            Sdgab(1,2)=Semb3(i);
        elseif sub_eksis(4)==2 && sub_eksis(5)==2 %disisipkan pada theta_deg (QR)
            theta=qim_emb(i);
            [Semb1(i), Semb2(i),Semb3(i)] = ptc_3d(r, deg2rad(theta), phi);
            Sdgab(1,1)=Semb1(i);
            Sdgab(2,2)=Semb2(i);
            Sdgab(1,2)=Semb3(i);
        elseif sub_eksis(4)==2 && sub_eksis(5)==3 %disisipkan pada theta_rad (QR)
            theta=qim_emb(i);
            [Semb1(i), Semb2(i),Semb3(i)] = ptc_3d(r, theta, phi);
            Sdgab(1,1)=Semb1(i);
            Sdgab(2,2)=Semb2(i);
            Sdgab(1,2)=Semb3(i);
        elseif sub_eksis(4)==2 && sub_eksis(5)==4 %disisipkan pada phi_deg (QR)
            phi=qim_emb(i);
            [Semb1(i), Semb2(i),Semb3(i)] = ptc_3d(r, theta, deg2rad(phi));
            Sdgab(1,1)=Semb1(i);
            Sdgab(2,2)=Semb2(i);
            Sdgab(1,2)=Semb3(i);
        elseif sub_eksis(4)==2 && sub_eksis(5)==5 %disisipkan pada phi_rad (QR)
            phi=qim_emb(i);
            [Semb1(i), Semb2(i),Semb3(i)] = ptc_3d(r, theta, phi);
            Sdgab(1,1)=Semb1(i);
            Sdgab(2,2)=Semb2(i);
            Sdgab(1,2)=Semb3(i);
        end


        % %QIM extraction
        % wt(i,1)=QIM_ext(Semb(i),nb,jenis);
        %Penggabungan nilai Singular yg disisipkan dg yg tdk
        if sub_eksis(4)==1 %ISVD
            Sisvd=U*Sdgab*V.';
        elseif sub_eksis(4)==2
            Sisvd=U*Sdgab;
        end

        Sisvd1=reshape(Sisvd,[blok^2 1]);


        %ISVD
        %Menggabungkan subband yg disisipkan dg yg tdk disisipkan
        Swgab=Sw;

        if sub_eksis(3)==0
            s_idst=Sisvd1;
        elseif sub_eksis(3)==1
            s_idst=dct(Sisvd1);
        end

        if sub_eksis(2)==0
            %IDST
            Swgab(sub,:)=s_idst;

        elseif sub_eksis(2)==1
            %IDST
            Swgab(sub,:)=idst(s_idst);
        end

        if sub_eksis(1)==1
            %ISWT
            st=iswt(Swgab,'haar');
        elseif sub_eksis(1)==0
            st=Swgab;
        elseif sub_eksis(1)==2
            cA=idwt(Swgab(1,:),Swgab(2,:),'haar');
            cD=idwt(Swgab(3,:),Swgab(4,:),'haar');
            st=idwt(cA,cD,'haar');
        end

        %Penggabungan Segmen
        xw(L*i-L+1:L*i,1)=st/alfa;
    end
        xwb=(xw+1)*255/2;

    So=[];U0=[];U1=[];V0=[];V1=[];pn0=[];pn1=[];wt=[];
    xwr=reshape((xwb),[size(host2d,1),size(host2d,2)]);
%     mse=mean(mean((host1d-xwr(:)).^2));
%     psnr=10*log10(255^2/mse);
% psnr
    % xwr=reshape(round(xw),[size(host2d,1),size(host2d,2)]);

elseif mode==1 %SS

    if sub_eksis(1)==1 || sub_eksis(1)==0 %SWT
        %Generate pncode
        pn0=2*randi([0 1],B,B)-1;
        pn1=2*randi([0 1],B,B)-1;
    elseif sub_eksis(1)==2 %DWT
        %Generate pncode
        pn0=2*randi([0 1],B/2,B/2)-1;
        pn1=2*randi([0 1],B/2,B/2)-1;
    end

    if sub_eksis(4)==1 %SVD
        % pn0=randi([0 1],B,B);
        % pn1=randi([0 1],B,B);
        [U0,S0,V0]=svd(pn0);
        [U1,S1,V1]=svd(pn1);
    elseif sub_eksis(4)==2

        [U0,S0]=qr(pn0);
        [U1,S1]=qr(pn1);
        V0=[];V1=[];
        % elseif sub_eksis(4)==0
        %     S=Swdct2;
    end

    % % Embedding dengan Metode SS
    % for i=1:length(logonrz)
    %     xw(L*i-L+1:L*i,1)=host1d(L*i-L+1:L*i)+alfa*logonrz(i)*pn;
    % end
    % mean(mean((u*s*v.').*pn))
    % mean(mean((u*s*v.').*pn1))
    % mean(mean((u1*s1*v1.').*pn))
    % mean(mean((u1*s1*v1.').*pn1))
    % mean(mean((u*s*v.').*pn))
    % mean(mean((u*s1*v.').*pn1))

    %%%%%%%%%%%%%%EMBEDDING%%%%%%%%%%%%%%
    for i=1:size(host1d,1)/L
        %Segmentasi
        seg=alfa*host1d(L*i-L+1:L*i,1);%segmentasi
        %Wavelet Transform
        if sub_eksis(1)==1 %SWT
            Sw=swt(seg,Ldek,'haar');
        elseif sub_eksis(1)==0
            Sw=seg.';subband=1;
        elseif sub_eksis(1)==2
            [A,D]=dwt(seg','haar');
            [LL,LH]=dwt(A,'haar');
            [HL,HH]=dwt(D,'haar');
            Sw=[LL;LH;HL;HH];
        end
        %DST
        if sub_eksis(2)==0
            sdst=Sw(subband,:);
        elseif sub_eksis(2)==1 %DST
            sdst=dst(Sw(subband,:));
        end
        %DCT
        if sub_eksis(3)==0
            Swdct=sdst;
        elseif sub_eksis(3)==1 %DCT
            Swdct=dct(sdst);
        end
        %1D ke 2D
        if sub_eksis(1)==1 || sub_eksis(1)==0 %SWT

            Swdct2=reshape(Swdct,[B B]);
        elseif sub_eksis(1)==2 %DWT
            Swdct2=reshape(Swdct,[B/2 B/2]);

        end
        %SVD atau QR
        if sub_eksis(4)==1 %SVD
            [U,S,V]=svd(Swdct2);
        elseif sub_eksis(4)==2
            [U,S]=qr(Swdct2);
        elseif sub_eksis(4)==0
            S=Swdct2;
        end
        if mode==0
            %untuk QIM
            %CPT dengan berbagai parameter
            if sub_eksis(5)==0 %Tanpa CPT
                swcpt(i)=S(1);
            elseif sub_eksis(4)==1 && sub_eksis(5)==1 %disisipkan pada r (SVD)
                %SVD-CPT
                [r, theta, theta_deg] = ctp(S(1),S(2));
                swcpt(i)=r;
            elseif sub_eksis(4)==1 && sub_eksis(5)==2 %disisipkan pada theta_deg (SVD)
                [r, theta_rad, theta_deg] = ctp(S(1),S(2));
                swcpt(i)=theta_deg;
            elseif sub_eksis(4)==1 && sub_eksis(5)==3 %disisipkan pada theta_rad (SVD)
                [r, theta_rad, theta_deg] = ctp(S(1),S(2));
                swcpt(i)=theta_rad;
            elseif sub_eksis(4)==2 && sub_eksis(5)==1 %disisipkan pada r (QR)
                [r, theta, phi] = ctp_3d(S(1,1), S(2,2), S(1,2));
                swcpt(i)=r;
            elseif sub_eksis(4)==2 && sub_eksis(5)==2 %disisipkan pada theta_deg (QR)
                [r, theta, phi] = ctp_3d(S(1,1), S(2,2), S(1,2));
                swcpt(i)=rad2deg(theta);
            elseif sub_eksis(4)==2 && sub_eksis(5)==3 %disisipkan pada theta_rad (QR)
                [r, theta, phi] = ctp_3d(S(1,1), S(2,2), S(1,2));
                swcpt(i)=(theta);
            elseif sub_eksis(4)==2 && sub_eksis(5)==4 %disisipkan pada phi_deg (QR)
                [r, theta, phi] = ctp_3d(S(1,1), S(2,2), S(1,2));
                swcpt(i)=rad2deg(phi);
            elseif sub_eksis(4)==2 && sub_eksis(5)==5 %disisipkan pada phi_rad (QR)
                [r, theta, phi] = ctp_3d(S(1,1), S(2,2), S(1,2));
                swcpt(i)=(phi);
            end

            %QIM embedding
            % r_emb(i)=QIM_emb(r,logobw1d(i),nb,jenis);
            % theta_emb(i)=QIM_emb(theta_deg,logobw1d(i),nb,jenis);
            qim_emb(i)=QIM_emb(swcpt(i),logobw1d(i),nb,jenis);
            % qim_emb(i)=swcpt(i);
            wt(i,1)=QIM_ext(qim_emb(i),nb,jenis);
            Sdgab=S;
            if sub_eksis(5)==0 %Tanpa CPT
                Sdgab(1)=qim_emb(i);


            elseif sub_eksis(4)==1 && sub_eksis(5)==1 %disisipkan pada r (SVD)
                %SVD-CPT
                r=qim_emb(i);
                [Semb1(i), Semb2(i)] = ptc(r, theta);
                Sdgab(1)=Semb1(i);
                Sdgab(2)=Semb2(i);

            elseif sub_eksis(4)==1 && sub_eksis(5)==2 %disisipkan pada theta_deg (SVD)
                theta=qim_emb(i);
                [Semb1(i), Semb2(i)] = ptc(r,theta);
                Sdgab(1)=Semb1(i);
                Sdgab(2)=Semb2(i);
            elseif sub_eksis(4)==1 && sub_eksis(5)==3 %disisipkan pada theta_rad (SVD)
                theta=qim_emb(i);
                [Semb1(i), Semb2(i)] = ptc(r,rad2deg(theta));
                Sdgab(1)=Semb1(i);
                Sdgab(2)=Semb2(i);
            elseif sub_eksis(4)==2 && sub_eksis(5)==1 %disisipkan pada r (QR)
                r=qim_emb(i);
                [Semb1(i), Semb2(i),Semb3(i)] = ptc_3d(r, theta, phi);
                Sdgab(1,1)=Semb1(i);
                Sdgab(2,2)=Semb2(i);
                Sdgab(1,2)=Semb3(i);
            elseif sub_eksis(4)==2 && sub_eksis(5)==2 %disisipkan pada theta_deg (QR)
                theta=qim_emb(i);
                [Semb1(i), Semb2(i),Semb3(i)] = ptc_3d(r, deg2rad(theta), phi);
                Sdgab(1,1)=Semb1(i);
                Sdgab(2,2)=Semb2(i);
                Sdgab(1,2)=Semb3(i);
            elseif sub_eksis(4)==2 && sub_eksis(5)==3 %disisipkan pada theta_rad (QR)
                theta=qim_emb(i);
                [Semb1(i), Semb2(i),Semb3(i)] = ptc_3d(r, theta, phi);
                Sdgab(1,1)=Semb1(i);
                Sdgab(2,2)=Semb2(i);
                Sdgab(1,2)=Semb3(i);
            elseif sub_eksis(4)==2 && sub_eksis(5)==4 %disisipkan pada phi_deg (QR)
                phi=qim_emb(i);
                [Semb1(i), Semb2(i),Semb3(i)] = ptc_3d(r, theta, deg2rad(phi));
                Sdgab(1,1)=Semb1(i);
                Sdgab(2,2)=Semb2(i);
                Sdgab(1,2)=Semb3(i);
            elseif sub_eksis(4)==2 && sub_eksis(5)==5 %disisipkan pada phi_rad (QR)
                phi=qim_emb(i);
                [Semb1(i), Semb2(i),Semb3(i)] = ptc_3d(r, theta, phi);
                Sdgab(1,1)=Semb1(i);
                Sdgab(2,2)=Semb2(i);
                Sdgab(1,2)=Semb3(i);
            end

        elseif mode==1
            %Spread Spectrum
            if sub_eksis(4)==1 %SVD
                So(i).data={S};
                if logonrz(i)==1
                    Sdgab=S+alfass*S1;
                else
                    Sdgab=S+alfass*S0;
                end
                X00=(mean(mean(abs((U0*(Sdgab-cell2mat(So(i).data))/alfass*V0.').*pn0))));
                X11=(mean(mean(abs((U1*(Sdgab-cell2mat(So(i).data))/alfass*V1.').*pn1))));
                if X00>=X11
                    wt(i,1)=0;
                else
                    wt(i,1)=1;
                end
            elseif sub_eksis(4)==2 %QR
                So(i).data={S};
                if logonrz(i)==1
                    Sdgab=S+alfass*S1;
                else
                    Sdgab=S+alfass*S0;
                end
                X00=(mean(mean(abs((U0*(Sdgab-cell2mat(So(i).data))/alfass).*pn0))));
                X11=(mean(mean(abs((U1*(Sdgab-cell2mat(So(i).data))/alfass).*pn1))));
                if X00>=X11
                    wt(i,1)=0;
                else
                    wt(i,1)=1;
                end

            elseif sub_eksis(4)==0 %

            end

        end


        % %QIM extraction
        % wt(i,1)=QIM_ext(Semb(i),nb,jenis);
        %Penggabungan nilai Singular yg disisipkan dg yg tdk
        if sub_eksis(4)==1 %ISVD
            Sisvd=U*Sdgab*V.';
        elseif sub_eksis(4)==2%QR Rec.
            Sisvd=U*Sdgab;
        elseif sub_eksis(4)==0
            Sisvd=Sdgab;

        end
        %2D ke 1D
        if sub_eksis(1)==1 || sub_eksis(1)==0 %SWT

            Sisvd1=reshape(Sisvd,[(B)^2 1]);
        elseif sub_eksis(1)==2  %DWT
            Sisvd1=reshape(Sisvd,[(B/2)^2 1]);
        end
        %IDCT
        %Menggabungkan subband yg disisipkan dg yg tdk disisipkan
        Swgab=Sw;
        if sub_eksis(3)==0
            s_idst=Sisvd1;
        elseif sub_eksis(3)==1
            s_idst=idct(Sisvd1);
        end

        %IDST
        if sub_eksis(2)==0
            %IDST
            Swgab(subband,:)=s_idst;

        elseif sub_eksis(2)==1
            %IDST
            Swgab(subband,:)=idst(s_idst);
        end

        %Wavelet Reconstruction
        if sub_eksis(1)==1
            %ISWT
            st=iswt(Swgab,'haar');
        elseif sub_eksis(1)==0
            st=Swgab;
        elseif sub_eksis(1)==2
            cA=idwt(Swgab(1,:),Swgab(2,:),'haar');
            cD=idwt(Swgab(3,:),Swgab(4,:),'haar');
            st=idwt(cA,cD,'haar');
        end

        %Penggabungan Segmen
        xw(L*i-L+1:L*i,1)=st/alfa;
    end


xwr=reshape(round(xw),[size(host2d,1),size(host2d,2)]);

end