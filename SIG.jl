using PyCall,PyPlot,Images,LsqFit,LinearAlgebra,StatsBase,FFTW,FITSIO,HDF5,Statistics,ImageFiltering
using LazCore,LazCyvecd,LazInstaller,LazIO,LazType,LazHY,Lazsyn,ath2h5
@pyimport numpy as np
@pyimport glob as glob

function getMF(d::Cube,jb::Cube,kb::Cube)
    nx,ny,nz=size(d);
    By=zeros(ny,nz);
    Bx=zeros(ny,nz)
    for k in 1:nz, j in 1:ny, i in 1:nx
        Bx[j,k]=d[i,j,k]*jb[i,j,k];
        By[j,k]=d[i,j,k]*kb[i,j,k]
    end
    return Bx,By
end

dbs=glob.glob("g:/run_new/*.007")

for n in 1:10

    n=10
    f=h5open(dbs[n]);
    ib=read(f,"i_mag_field");
    jb=read(f,"j_mag_field");
    kb=read(f,"k_mag_field");
    d=read(f,"gas_density");

    nu=Inf;
    Is,Qs,Us,Ps=syn_emissionX(d,ib,jb,kb,3,nu)
    P=abs.(Ps);
    phi=0.5*atan.(Us,Qs);
    #Ii,Qi,Ui = getdIQU(d, jb, kb);
    #phi=0.5*atan.(Ui,Qi);

    dn=22;
    ix,iy=sobel_conv_2d(Is);
    #pna,pns=sban2d(cos.(phi),sin.(phi),dn);
    pna=avB2dx(cos.(phi),sin.(phi),dn);
    ina,ins=sban2d(ix,iy,dn);
    ina.+=pi/2;
    pna.+=pi/2;
    AM(ina,pna)

    nx,ny=size(Is);
    Xd,Yd=np.meshgrid(div(dn,2)+1:dn:div(dn,2)+ny,div(dn,2)+1:dn:div(dn,2)+nx);

    figure(tight_layout="true")
    imshow(log10.(Is'./mean(Is)),origin="lower",cmap="coolwarm");
    tick_params(direction="in",labelsize=0);
    #imshow(log10.(Proj(d,3))',origin="lower", cmap="jet");
    cb=colorbar(pad=0,shrink=1)
    cb.ax.tick_params(labelsize=14,direction="in")
    cb.set_label(label="\$I/\\langle I\\rangle\$",size=16)
    #clim(0,1.6)
    quiver(Yd,Xd,cos.(pna),sin.(pna),headwidth=0,scale=35,color="b", label="Polarization")
    quiver(Yd,Xd,-cos.(pna),-sin.(pna),headwidth=0,scale=35,color="b")
    quiver(Yd,Xd,cos.(ina),sin.(ina),headwidth=0,scale=35,color="r")
    quiver(Yd,Xd,-cos.(ina),-sin.(ina),headwidth=0,scale=35,color="r", label="SIGs")

end

