function Moverlapped = foverlap(Mgamma,Mphi)
max_i=size(Mgamma,1);
max_j=size(Mphi,1);
Nframes= min(size(Mgamma,2),size(Mphi,2));
Moverlapped=zeros(max_i,max_j);
for i=1:max_i
    for j=1:max_j
        for k=1:Nframes
            if (Mgamma(i,k)~=0)&&(Mphi(j,k)~=0)
                Moverlapped(i,j)=Moverlapped(i,j)+1;
            end
        end
    end
end
end