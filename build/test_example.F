!$acc parallel loop async
      do j=0,jj
        do l=1,isp(j)
c
        do i=max(0,ifp(j,l)),min(ii,ilp(j,l))
          util1(i,j)=0.
          util2(i,j)=0.
        enddo
c
        do k=1,kk
          kn=k+nn
          do i=max(0,ifp(j,l)),min(ii,ilp(j,l))
            pbotl=max(p(i,j,k+1),p(i,j,kk+1)-thkbop)
            ptopl=max(p(i,j,k  ),p(i,j,kk+1)-thkbop)
            util1(i,j)=util1(i,j)+(u(i,j,kn)+u(i+1,j,kn))*(pbotl-ptopl)
            util2(i,j)=util2(i,j)+(v(i,j,kn)+v(i,j+1,kn))*(pbotl-ptopl)
          enddo
        enddo
c
        do i=max(0,ifp(j,l)),min(ii,ilp(j,l))
          ubot=(ubflxs_p(i  ,j,n)
     .          /max(epsil,pbu(i  ,j,n)*scuy(i  ,j))
     .         +ubflxs_p(i+1,j,n)
     .          /max(epsil,pbu(i+1,j,n)*scuy(i+1,j)))*tsfac
     .        +util1(i,j)/thkbop
          vbot=(vbflxs_p(i,j  ,n)
     .          /max(epsil,pbv(i,j  ,n)*scvx(i,j  ))
     .         +vbflxs_p(i,j+1,n)
     .          /max(epsil,pbv(i,j+1,n)*scvx(i,j+1)))*tsfac
     .        +util2(i,j)/thkbop
          ubbl=.5*sqrt(ubot*ubot+vbot*vbot)
          q=cb*(ubbl+cbar)
          drag(i,j)=q*g/(alpha0*thkbop)
          ustarb(i,j)=sqrt(q*ubbl)
        enddo
c
        enddo
      enddo
!$acc end parallel loop
c$OMP END PARALLEL DO
c
c --- store r.h.s. of barotropic u/v eqn. in -ubrhs,vbrhs-
c --- store wind forcing in -stresx,stresy-
c
c$OMP PARALLEL DO PRIVATE(l,i)
!$acc parallel loop async
      do j=1,jj
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
          ubrhs(i,j)=ubcors_p(i,j)*tsfac
        enddo
        enddo
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          vbrhs(i,j)=vbcors_p(i,j)*tsfac
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
!$acc end parallel loop
c
c$OMP PARALLEL DO PRIVATE(i)
!$acc parallel loop async collapse(2)
      do j=0,jj+1
        do i=0,ii+1
          dl2u(i,j)=0.
          dl2v(i,j)=0.
        enddo
      enddo
c$OMP END PARALLEL DO
!$acc end parallel loop
c
      do k=1,kk
!$acc parallel loop
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=-1,jj+2
        km=k+mm
          do l=1,isu(j)
          do i=max(-1,ifu(j,l)),min(ii+2,ilu(j,l))
            pu(i,j,k+1)=pu(i,j,k)+dpu(i,j,km)
          enddo
          enddo
          do l=1,isv(j)
          do i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
            pv(i,j,k+1)=pv(i,j,k)+dpv(i,j,km)
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
!$acc end parallel loop
      enddo
c
!$acc wait
