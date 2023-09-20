using PyCall,PyPlot,LsqFit,Images,LinearAlgebra,StatsBase,FITSIO,Statistics,ImageFiltering
using LazCore,LazType
@pyimport numpy as np

function imfilter_gaussian(d::Mat,p::Number)
    im=imfilter(d,Kernel.gaussian(p),NA())
    return im
end

function de_resolution(b::Mat,dn)
    nx,ny=size(b)
    bx=zeros(div(nx,dn),div(ny,dn));
   for i in 1:div(nx,dn),j in 1:div(ny,dn)
        a=b[i*dn-dn+1:i*dn,j*dn-dn+1:j*dn]
        bx[i,j]=mean(a)
   end
      return bx
end

function avB2dx(Ax::Mat,Ay::Mat,dn)
	nx,ny=size(Ax)
	Ana=zeros(div(nx,dn),div(ny,dn));
	for  j in 1:div(ny,dn),i in 1:div(nx,dn)
	    is=(i-1)*dn+1;
	    ie=i*dn;
	    js=(j-1)*dn+1;
	    je=j*dn;
	    Axx=Ax[is:ie,js:je];
	    Ayy=Ay[is:ie,js:je];
	    Ana[i,j]=atan.(mean(Ayy[.~isnan.(Ayy)])/mean(Axx[.~isnan.(Ayy)]));
	end
	return Ana
end

function sban2d_mv(Ax::Mat,Ay::Mat,dn)
	nx,ny=size(Ax)
	Ana=zeros(nx,ny);
	Ans=zeros(nx,ny);
    for  j in dn/2:ny-dn/2,i in dn/2:nx-dn/2
        i=Int(i);
        j=Int(j);
        is=Int(i-dn/2+1);
        ie=Int(i+dn/2);
        js=Int(j-dn/2+1);
        je=Int(j+dn/2);
        Axx=Ax[is:ie,js:je];
        Ayy=Ay[is:ie,js:je];
        #Axx=Axx[.~isnan.(Axx)];
        #Ayy=Ayy[.~isnan.(Ayy)];
        if(length(Axx[.~isnan.(Axx)])>100)
            binsize=100;
            Apeak,Adisp=fit_gaussian_2d(Axx,Ayy,binsize);
            Ana[Int(i),Int(j)]=Apeak;
            Ans[Int(i),Int(j)]=Adisp;
        else
            Ana[Int(i),Int(j)]=NaN;
            Ans[Int(i),Int(j)]=NaN;
        #println("j="*string(j)*";i="*string(i))
        end
    end
    Ana[1:Int(dn/2-1),:].=NaN;
    Ana[:,1:Int(dn/2-1)].=NaN;
    Ana[Int(nx-dn/2+1):nx,:].=NaN;
    Ana[:,Int(ny-dn/2+1):ny].=NaN;
    Ans[isnan.(Ana)].=NaN; 
	return Ana,Ans
end

# read data
f=FITS("d:/ElGordo_hres_sync_blanked.fits");
d=read(f[1]);
header=read_header(f[1]);
d=d[:,:,1,1];

# calculate coordinates
CRPIX1b=read_header(f[1])["CRPIX1"];
CRVAL1b=read_header(f[1])["CRVAL1"];
CDELT1b=read_header(f[1])["CDELT1"];
CRPIX2b=read_header(f[1])["CRPIX2"];
CRVAL2b=read_header(f[1])["CRVAL2"];
CDELT2b=read_header(f[1])["CDELT2"];

nx,ny=size(d)
ra=zeros(nx);
dec=zeros(ny);
for i in 1:nx
	ra[i]=CRVAL1b+(i-CRPIX1b)*CDELT1b
end
for j in 1:ny
    dec[j]=CRVAL2b+(j-CRPIX2b)*CDELT2b
end

#calculate gradient
Ii=imfilter_gaussian(d,0);

dn=20;
ix,iy=sobel_conv_2d(Ii);
# ina: averaged gadient, ins: uncertainty
ina,ins=sban2d_mv(ix,iy,dn);

# constructing Stokes Q, U parameters
# calculating gradient angle psi
Qa=Ii.*cos.(2.0.*ina);
Ua=Ii.*sin.(2.0.*ina);
Qb=imfilter_gaussian(Qa,4);
Ub=imfilter_gaussian(Ua,4);
psi=0.5*atan.(Ub,Qb);
psi.+=pi/2;
psi[isnan.(Ii)].=NaN;

#calculting uncertainty psie
xx=div.(ins,pi);
for i in 1:size(xx)[1],j in 1:size(xx)[2]
	ins[i,j]-=xx[i,j]*pi;
end

noise=1.5e-6; # noise levels
ce=abs.(2*sin.(2.0.*ina).*ins)
se=abs.(2*cos.(2.0.*ina).*ins)
Qe=abs.(Qa).*sqrt.((noise./Ii).^2.0.+(ce./cos.(2.0.*ina)).^2);
Ue=abs.(Ua).*sqrt.((noise./Ii).^2.0.+(se./sin.(2.0.*ina)).^2);

UQ=Ub./Qb;
UQe=abs.(UQ).*sqrt.((Qe./Qb).^2.0.+(Ue./Ub).^2.0);
psie=0.5.*UQe./(1.0.+UQ.^2.0);
psie[isnan.(psi)].=NaN;

xx=div.(psie,pi);
for i in 1:size(xx)[1],j in 1:size(xx)[2]
    psie[i,j]-=xx[i,j]*pi;
end

# removing pixels, where uncertainty is larger than 30 degrees.
AA=copy(psi);
AA[.~isnan.(psi)].=1;
AA[isnan.(psi)].=0;
BB=imfill(Bool.(AA),(0,1000));
psi[Int.(BB).==0].=NaN;
psi[psie.*180/pi.>30].=NaN;

# averaging the gradient for visulization purpose
dnn=6;
pna=avB2dx(cos.(psi[26:290,26:220]),sin.(psi[26:290,26:220]),dnn);
dI=de_resolution(Ii[26:290,26:220],dnn);
pna[isnan.(dI)].=NaN;

nx,ny=size(Ii[26:290,26:220])
Xd,Yd=np.meshgrid(div(dnn,2)+1:dnn:div(dnn,2)+ny,div(dnn,2)+1:dnn:div(dnn,2)+nx);

Ia=log10.(Ii[26:290,26:220]);
Ia[isnan.(Ia)].=-1000;

figure(figsize=[10,5],tight_layout="true");
imshow(Ia',origin="lower", cmap="inferno");
cb=colorbar(pad=0,shrink=1);
cb.ax.tick_params(labelsize=10,direction="in");
cb.set_label(label="log\$_{10}\$(Intensity) [Jy/beam]",size=10);
clim(-5.5,-2.0);
cpp=contour(log10.(Ia)',levels=[-1000,-4.9,-4.5,-4,-3.5,-3,-2.5,-2,-1.5],alpha=0.3,cmap="Greys");
clabel(cpp, inline=2, fontsize=10);
xlabel("R.A.(J2000) [degree]",size=10);
ylabel("Dec.(J2000) [degree]",size=10);
tick_params(direction="in");
xticks([0,div(nx,4)*1,div(nx,4)*2,div(nx,4)*3,div(nx,4)*4-1],[round(ra[1+25];digits=3),round(ra[div(nx,4)+25];digits=3),round(ra[div(nx,4)*2+25];digits=3),round(ra[div(nx,4)*3+25];digits=3),round(ra[div(nx,4)*4+25];digits=3)],size=10);
yticks([0,div(ny,4)*1,div(ny,4)*2,div(ny,4)*3,div(ny,4)*4-1],[round(dec[1+25];digits=3),round(dec[div(ny,4)+25];digits=3),round(dec[div(ny,4)*2+25];digits=3),round(dec[div(ny,4)*3+25];digits=3),round(dec[div(ny,4)*4+25];digits=3)],size=10);
quiver(Yd[1:end-1,1:end-1],Xd[1:end-1,1:end-1],cos.(pna),sin.(pna),headwidth=0,scale=8,color="w",scale_units="inches")
quiver(Yd[1:end-1,1:end-1],Xd[1:end-1,1:end-1],-cos.(pna),-sin.(pna),headwidth=0,scale=8,color="w",scale_units="inches")


