function wt=extraction_image(attacked_img,logobw1d,L,alfa,alfass,sub_eksis,Ldek,subband,B,nb,jenis,mode,So,U0,U1,V0,V1,pn0,pn1)
if mode==0
    sub=subband;
    % logobw1d=logonrz;
      if sub_eksis(1)==2
    blok=B/2;
    else
        blok=B;
    end
    xw2=reshape(attacked_img,[size(attacked_img,1)*size(attacked_img,2) 1]);
    
        %normalisasi
    xw3=2*xw2/255-1;

    for i=1:length(logobw1d)

        %Segmentasi
        seg=alfa*xw3(L*i-L+1:L*i,1);
        %SWT
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
        Swdct2=reshape(Swdct,[blok blok]);

        %SVD
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
        %QIM extraction
        % wt(i,1)=QIM_ext(rext(i),nb,jenis);
        % wt(i,1)=QIM_ext(theta_deg(i),nb,jenis);
        wt(i,1)=QIM_ext(swcpt(i),nb,jenis);
    end
elseif mode==1
    % xwn=attacked_img;
    xwn=reshape(attacked_img,[size(attacked_img,1)*size(attacked_img,2),1]);

    for i=1:size(xwn,1)/L
        % for i=1:1

        %Segmentasi
        seg=alfa*xwn(L*i-L+1:L*i,1);
        %SWT
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
        elseif sub_eksis(1)==2  %DWT
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
        if mode==0 %QIM
            %CPT dengan berbagai parameter
            if sub_eksis(5)==0 %Tanpa CPT
                swcptr(i)=S(1);
            elseif sub_eksis(4)==1 && sub_eksis(5)==1 %disisipkan pada r (SVD)
                %SVD-CPT
                [r, theta, theta_deg] = ctp(S(1),S(2));
                swcptr(i)=r;
            elseif sub_eksis(4)==1 && sub_eksis(5)==2 %disisipkan pada theta_deg (SVD)
                [r, theta_rad, theta_deg] = ctp(S(1),S(2));
                swcptr(i)=theta_deg;
            elseif sub_eksis(4)==1 && sub_eksis(5)==3 %disisipkan pada theta_rad (SVD)
                [r, theta_rad, theta_deg] = ctp(S(1),S(2));
                swcptr(i)=theta_rad;
            elseif sub_eksis(4)==2 && sub_eksis(5)==1 %disisipkan pada r (QR)
                [r, theta, phi] = ctp_3d(S(1,1), S(2,2), S(1,2));
                swcptr(i)=r;
            elseif sub_eksis(4)==2 && sub_eksis(5)==2 %disisipkan pada theta_deg (QR)
                [r, theta, phi] = ctp_3d(S(1,1), S(2,2), S(1,2));
                swcptr(i)=rad2deg(theta);
            elseif sub_eksis(4)==2 && sub_eksis(5)==3 %disisipkan pada theta_rad (QR)
                [r, theta, phi] = ctp_3d(S(1,1), S(2,2), S(1,2));
                swcptr(i)=(theta);
            elseif sub_eksis(4)==2 && sub_eksis(5)==4 %disisipkan pada phi_deg (QR)
                [r, theta, phi] = ctp_3d(S(1,1), S(2,2), S(1,2));
                swcptr(i)=rad2deg(phi);
            elseif sub_eksis(4)==2 && sub_eksis(5)==5 %disisipkan pada phi_rad (QR)
                [r, theta, phi] = ctp_3d(S(1,1), S(2,2), S(1,2));
                swcptr(i)=(phi);
            end
            %QIM extraction
            % wt(i,1)=QIM_ext(rext(i),nb,jenis);
            % wt(i,1)=QIM_ext(theta_deg(i),nb,jenis);
            wt(i,1)=QIM_ext(swcptr(i),nb,jenis);
        elseif mode==1 %Spread Spectrum
            if sub_eksis(4)==1 %SVD
                X0=(mean(mean(abs((U0*(S-cell2mat(So(i).data))/alfass*V0.').*pn0))));
                X1=(mean(mean(abs((U1*(S-cell2mat(So(i).data))/alfass*V1.').*pn1))));
                if X0>=X1
                    wt(i,1)=0;
                else
                    wt(i,1)=1;
                end
            elseif sub_eksis(4)==2 %QR
                X0=(mean(mean(abs((U0*(S-cell2mat(So(i).data))/alfass).*pn0))));
                X1=(mean(mean(abs((U1*(S-cell2mat(So(i).data))/alfass).*pn1))));
                if X0>=X1
                    wt(i,1)=0;
                else
                    wt(i,1)=1;
                end
            elseif sub_eksis(4)==0 %

            end


        end
    end

end