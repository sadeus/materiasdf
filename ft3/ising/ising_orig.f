ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   simulacion de monte carlo de modeldo de ising bidimensional
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
       double precision dseed

        common /datain/ nx,ny,xj,b
        common /param/  nterm,ngrupo,nfrec,nsize
        common /seed/   dseed
        dimension ratio(-4:4,-1:1)
c
        real                 xjj,xjjjj
c...........................................parametros de las
c
c        nnnn = 400
c        nx=nnnn
c        ny=nnnn
c        xjj=0.001
c        xjjj=xjj
c             do while(xjj.gt.0.0)
             iflag=1
             call ising_1(xjj)
            xj=xjj
                 if(xjj.gt.0.0) then
                 call ising00(iflag,ratio)   ! genera vector para salto
                 call ising01(ratio)
                 xjjj=xjj
                 endif
              xjjj=xjj
c              enddo
c
        close(20)
        open(21,file='randomdt.dat',status='old') ! data seed
	rewind(21)
        write(21,*) dseed
        close(21)
        close(22)
        close(23)
c
        stop
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        data input
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ising_1(xjj)
c
c   nx      =   dimx lattice,
c   ny      =   dimy lattice,
c   nterm   =   numero de pasos de termalizacion
c   ngrupo  =   cantidad de subcadenas
c   nfrec   =   separacion entre mediciones
c   nsize   =   cantidad de mediciones en un grupo
c   xj=fuerza del potencial
c   b=campo externo
c............................................................
c
c       double precision dseed,random
        double precision dseed
        real                 xjj,xjn
        character*12        aaa
c............................................................
        common /datain/ nx,ny,xj,b
        common /param/  nterm,ngrupo,nfrec,nsize
        common /seed/   dseed
        dimension ratio(-4:4,-1:1)
c............................................................
        open(20,file='isingdat.dat',status='old') ! data 
        read(20,*) nx,ny
        read(20,*) nterm,ngrupo,nfrec
        read(20,*) iflag,nsize
        read(20,*) xjj
        read(20,*) b
        read(20,*) aaa
        close(20)
200     format(2i8)
201     format(3i8)
c        
        open(21,file='randomdt.dat',status='old') ! data seed
        read(21,*) dseed
        close(21)
c..............         data para test
c       nx=10
c       ny=10
c       nterm=40    ! termalizacion 
c       ngrupo=10   ! subcadenas
c       nfrec=10    ! muestras por subcadema
c       nsize=10    ! pasos intermedios
c       xjj=0.2
c       b=0.0
c       dseed=123457.d0
c       aaa="output01.dat"
c       write(6,*) 'aaa ',aaa
c
c.....................................
c         open(22,file='confprev.dat',status='unknown')
c.............          data output
        open(23,file=aaa,status='unknown')
        write(23,*) '# dimension de la red :',nx
        write(23,*) '# acoplamiento spines :',xjj
        write(23,*) '# campo externo       :',b
        write(23,*) '#                     '
        write(23,*) '#paso,  <mag.total>,<mag.cadena>,energ,calor esp'
        write(23,*) '#                     '
c
        open(21,file='randomdt.dat',status='old') ! data seed
        rewind(21)
        write(21,*) dseed
        close(21)
c
        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   iniciales
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ising00(iflag,ratio)
c
        double precision dseed,random
        parameter(idimx =500,idimy=500)
c        common /dataux/ratio(-4:4,-1:1)
        common /config/ispin(idimx,idimy)
        common /datain/ nx,ny,xj,b
        common /param/  nterm,ngrupo,nfrec,nsize
        dimension ratio(-4:4,-1:1)
c........................................
c........................................matrix para aceptacion
            do if=-4,4,2
            ratio(if,-1)=exp(2*(xj*if+b))
            ratio(if,1) =1./ratio(if,-1)
c            write(6,*) if,ratio(if,-1),ratio(if,1)
            enddo
c        pause
c........................................inicial random, si iflag=1
            if(iflag.eq.1) then
                do ix=1,nx
                    do iy=1,ny
                        if (random(dseed).gt..5) then
                        ispin(ix,iy)=1
                        else
                        ispin (ix,iy)=-1
                        endif
c                    write(6,*) ix,iy,ispin(ix,iy)
                    enddo
                enddo
            endif
c.........................................fin
        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   calculamos : energia magnetizacion suceptibilidad calor especifico
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ising01(ratio)
        parameter(idimx =500,idimy=500)
c........................................termalizacion
        common /config/ ispin(idimx,idimy)
        common /fismag/ xmag(3,3),xene(3,3),xce(3,3),xchi(3,3)
        common /param/  nterm,ngrupo,nfrec,nsize
        common /datain/ nx,ny,xj,b
        common/dataout/ prou,prod,e1,c1,xms1,es1,cs1
        common /screenn/ halfx,halfy,xwidth,yheight,cols,rows
        dimension ratio(-4:4,-1:1)
        data isw,igr,itot/1,2,3/
        data ival,isq,isi/1,2,3/
c........................................termalizacion
            DO iterm=1,nterm
            call mmc(rate,ratio)
            ENDDO
        call poncero(itot)                          ! ponemos a cero los totales
        more=ngrupo
        kk = 0
            DO igrp=ngrupo-more+1,ngrupo            ! loop sobre grupos
c 		write(6,*) 'igrp ',igrp
            call poncero(igr)                       ! cero promedios grupo
                DO iter=1,nfrec*nsize
                call mmc(rate,ratio)                ! sobre todas las part
                    IF (mod(iter,nfrec).eq.0) then  ! una salida
                    isweep=iter/nfrec               ! which sweep is it
                    call poncero(isw)               ! zero sweep totals
                    call sum                        ! sweep totals, add to group
                    END IF
                ENDDO
            kk=kk+1
            call promedio(igrp,xm1,xm2,e1,c1,xms1,es1,cs1) !CALC TOTAL AVERAGES
c                        prou = -xm1
c                        prod = -xm2
                         prou = xm1
                         prod = xm2
                        kk = igrp
            ENDDO
c
             xmore=0.
c        more = 0
c            IF (more .gt. 0) then
c            ngrupo=ngrupo+more
c            goto 15
c            END IF
c
      RETURN
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   hacemos un loop completo sobrela red
c   nnsum puede valer -4,-2,0,2,4
cccccccccccccccccccccccccccccccc
      subroutine mmc(rate,ratio)
c      
        double precision dseed,random,ran
        parameter(idimx =500,idimy=500)
        common /config/ispin(idimx,idimy)
        common /datain/nx,ny,xj,b
        common /param/  nterm,ngrupo,nfrec,nsize
        common /seed/dseed
        dimension ratio(-4:4,-1:1)
c
        nspin = nx*ny
        nacep=0
            do ix=1,nx                      ! loop sobre elementos (x)
            ixp1=ix+1                       ! vecino
            if (ix .eq. nx) ixp1=1          ! condiciones periodicas
            ixm1=ix-1
            if (ix .eq. 1) ixm1=nx          ! condiciones periodicas
                do iy=1,ny                  ! loop sobrew vecinos (y)
                iyp1=iy+1                   ! vecino
                if (iy .eq. ny) iyp1=1      ! cond.periodicas
                iym1=iy-1
                if (iy .eq. 1) iym1=ny
                nnsum=ispin(ix,iyp1)+ispin(ix,iym1)+ispin(ixp1,iy)
     &                  +ispin(ixm1,iy)     ! interaccion vecinos
                ran=random(dseed)           ! random
c
                    if (ran.lt.ratio(nnsum,ispin(ix,iy))) then   ! ya calc.
                    ispin(ix,iy)=-ispin(ix,iy) !cambiamos spin - exito!
                    nacep=nacep+1                   !otra aceptacion
                    end if
                enddo
            enddo
        rate=float(nacep)/float(nspin)      ! rate
c        write(6,*) 'rate : ', rate, nacep,nspin
c        pause
C
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   sumamnos    j - i - k - l       con lo que se suma sinrepeticion
c                   |   |   |
c                   j   i   k
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine sum()
c        
        parameter(idimx =500,idimy=500)
        common /config/ ispin(idimx,idimy)
        common /fismag/ xmag(3,3),xene(3,3),xce(3,3),xchi(3,3)
        common /datain/nx,ny,xj,b
        common /param/  nterm,ngrupo,nfrec,nsize
        data isw,igr,itot/1,2,3/
        data ival,isq,isi/1,2,3/
c............................................................
        pairs=0
            do iy=1,ny                  ! loop y
            iym1=iy-1
            if (iy .eq. 1) iym1=ny
                do ix=1,nx              ! loop x
                ixm1=ix-1
                if (ix .eq. 1) ixm1=nx
                pairs=pairs+ispin(ix,iy)*(ispin(ix,iym1)+ispin(ixm1,iy))
C.....................   magnet. = suma de spins sobre todos
                xmag(isw,ival)=xmag(isw,ival)+ispin(ix,iy)
                enddo
            enddo
c
      xmag(isw,isq)=xmag(isw,ival)**2
      xene(isw,ival)=-xj*pairs-b*xmag(isw,ival)   !eQ 8.18
      xene(isw,isq)=xene(isw,ival)**2
c
c....................     sumamos contributions a grupos
      xmag(igr,ival)=xmag(igr,ival)+xmag(isw,ival)
      xmag(igr,isq)=xmag(igr,isq)+xmag(isw,isq)
      xene(igr,ival)=xene(igr,ival)+xene(isw,ival)
      xene(igr,isq)=xene(igr,isq)+xene(isw,isq)
C
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc ccccccccccc
c se calculan los promedios de grupos y se suma a totales;
c se calculan los sigmas y se muestran resultados
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc ccccccccccc
c
        subroutine promedio(iigr,xm1,xm2,e1,c1,xms1,es1,cs1)
c
c       double precision dseed
c
        common /fismag/ xmag(3,3),xene(3,3),xce(3,3),xchi(3,3)
        common /param/  nterm,ngrupo,nfrec,nsize
        common /datain/ nx,ny,xj,b
c
        data isw,igr,itot/1,2,3/
        data ival,isq,isi/1,2,3/
c
        xnspin = float(nx)*float(ny)
C....................calcula prom grupos et  incert. de sumas grupos
            do iquant=ival,isq
            xmag(igr,iquant)=xmag(igr,iquant)/nsize
            xene(igr,iquant)=xene(igr,iquant)/nsize
            enddo
        xchi(igr,ival)=xmag(igr,isq)-xmag(igr,ival)**2
        xmag(igr,isi)=xchi(igr,ival)/nsize
        if (xmag(igr,isi) .lt. 0.) xmag(igr,isi)=0.
        xce(igr,ival)=xene(igr,isq)-xene(igr,ival)**2
        xene(igr,isi)=xce(igr,ival)/nsize
        if (xene(igr,isi) .lt. 0.) xene(igr,isi)=0.
        xchi(igr,isq)=xchi(igr,ival)**2
        xce(igr,isq)=xce(igr,ival)**2
C.......................................
            do iquant=ival,isi
            xmag(itot,iquant)=xmag(itot,iquant)+xmag(igr,iquant)
            xene(itot,iquant)=xene(itot,iquant)+xene(igr,iquant)
            xchi(itot,iquant)=xchi(itot,iquant)+xchi(igr,iquant)
            xce(itot,iquant)=xce(itot,iquant)+xce(igr,iquant)
            enddo
c......................................
C     promedios
        xm=xmag(itot,ival)/iigr
        xmsig1=(xmag(itot,isq)/iigr-xm**2)/iigr/nsize
        if (xmsig1 .lt. 0) xmsig1=0.
        xmsig1=sqrt(xmsig1)
        xmsig2=sqrt(xmag(itot,isi))/iigr
C........................................
        e=xene(itot,ival)/iigr
        esig1=(xene(itot,isq)/iigr-e**2)/iigr/nsize
        if (esig1 .lt. 0) esig1=0.
        esig1=sqrt(esig1)
        esig2=sqrt(xene(itot,isi))/iigr
C........................................
        sus=xchi(itot,ival)/iigr
        sussig=(xchi(itot,isq)/iigr-sus**2)/iigr
        if (sussig .lt. 0.) sussig=0.
        sussig=sqrt(sussig)
C.........................................
        c=xce(itot,ival)/iigr
        csig=(xce(itot,isq)/iigr-c**2)/iigr
        if (csig .lt. 0.) csig=0.
        csig=sqrt(csig)
c.........................................
        xm1 = xm/xnspin
        xm2 =xmag(igr,ival)/xnspin
        xms1=xmsig1/xnspin
        e1=e/xnspin
        es1=esig1/xnspin
        c1=c/xnspin
        cs1=csig/xnspin
c.........................................sera el res para graficos
        write(23,1111)iigr,xm1,xm2,e1,c1
1111    format(2x,i6,4(2x,f14.6))
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine poncero(ilevel)
        common/fismag/xmag(3,3),xene(3,3),xce(3,3),xchi(3,3)
            do iq=1,3
            xmag(ilevel,iq)=0.
            xene(ilevel,iq)=0.
            xchi(ilevel,iq)=0.
            xce(ilevel,iq)=0.
            enddo
        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C returns a uniformly distributed random number between 0 and 1
c
c      real function random(dseed)
      double precision function random(dseed)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision dseed,d2p31m,d2p31
      data d2p31m/2147483647.d0/
      data d2p31 /2147483711.d0/
c..........
      dseed = mod(16807.d0*dseed,d2p31m)
      random= dseed / d2p31
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
