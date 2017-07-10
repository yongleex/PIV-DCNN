%% understand the turbulent energy spectrum

function tes = turb_energy_spectrum(u,v)
    [sx,sy]=size(u);
  %  Energy calculation
    Fu=fftshift(fft2(u)); Fv=fftshift(fft2(v));
    phy11=conj(Fu).*Fu;phy22=conj(Fv).*Fv;
    phy=phy11+phy22;  
  %  mesh(log(phy+1));
  
  % wave number estimation
    cx=(sx+1)./2;cy=(sy+1)./2;
    for i=1:sx
        for j=1:sy
            A(i,j)=round(sqrt((i-cx).^2+(j-cy).^2));
        end
    end
    A(A==0)=1;
    
  % add the energy at the same wave number
    for i=min(A(:)):max(A(:))
        esp(i) = sum(phy(A==i));
    end
%    // plot(log(esp+1));
    tes = esp;
end
    
    