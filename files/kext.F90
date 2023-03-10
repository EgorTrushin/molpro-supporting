subroutine kext(icaa,ikaa,iq1,iq2,nt1,lent1,ipaa,ipee,nums,numt, &
     & ica,npkx,lenn,llcp,iadmax)
  use common_zahl
  use common_tapes
  use common_cvec
  use common_cthr
  use common_cprint
  use common_cpar
  use common_corb
  use common_clseg
  use common_clocal
  use common_cdirect
  use common_cbas
  use common_big
  use common_cmpp
  use cpfil
  use common_cjadr, nlgx => nlgc
  implicit double precision (a-h,o-z)
  logical sparse
  !
  dimension ipaa(8),ipee(8),ica(8,8,2),nums(8),numt(8)
  ipe1=ipee(1)-ipaa(1)+1
  dmax=0d0
  thr=0d0
  nzero=0
  sparse=.false.
  if(thrkcp.gt.0d0) then
    do i=1,llcp
      if(abs(q(icaa+i-1)).gt.1.d10) write(6,*) 'i=',i,q(icaa+i-1)
      dmax=max(dmax,abs(q(icaa+i-1)))
    end do
    thr=dmax*thrkcp
    do i=1,llcp
      if(abs(q(icaa-1+i)).lt.thr) then
        nzero=nzero+1
        q(icaa-1+i)=0d0
      end if
    end do
    spars=dble(nzero)/dble(llcp)
    sparse=spars.gt.sparfac
    if(sparse) write(iout,96) dmax,thr,nzero,llcp,spars,sparse
96  format(' KEXT: DMAX=',d12.2,'  THR=',d12.2,'  NZERO=',i8, &
         & '  LCP=',i8,'  SPARSE=',f10.3,' (',l1,')')
  end if
  !
  if(direct) then
    iadmax=max(iadmax,icorr(0))
    call dkext(q(icaa),q(ikaa),ipaa,ipee,ispa,ica,npkx,lenn)
    if (nprocs.gt.1) then
      call global_sync(1)
      mpp_state=1
    end if
    return
  end if
  ! .....n1,n2,n3,n4: space for integral matrices and intermediate quantities
  nmax=0
  do i=1,nsk
    nmax=max(nmax,ntb(i))
  end do
  n1=iq1
  n2=icorr(ipe1*nmax)
  n3=icorr(ipe1*nmax)
  n4=icorr(ipe1*nmax)
  iadmax=max(iadmax,n4+ipee(1)*nmax)
  if(nprocs.gt.1 .and. iprint(39).ge.1) then
    call global_sync(100)
    call add_mpptim('KEXT(TOT)',0,0)
    call add_mpptim('KEXT(SYNC)',0,0)
    call add_mpptim('KEXT(UNSYNC)',0,0)
  end if
  numnodes=nprocs
  inode=iprocs
  call fzero(q(ikaa),lenn)
  ikee=ikaa+lenn-1
  !
  ! basis integrals of type (aa/aa)
  !
  call add_mpptim('KEXTA',0,0)
  if(nsk.gt.1) call add_mpptim('KEXTB',0,-1)
  if(nsk.ge.4) call add_mpptim('KEXTC',0,-1)
  if (nprocs.gt.1) mpp_state=2
  do iska=1,nsk
    inta=nt(iska)
    if(inta.eq.0) cycle
    ip1a=ipaa(1)
    ip1e=ipee(1)
    nlga=inta*(inta+1)/2
    call kexta(q(icaa),q(ikaa),q(n1),q(n2),q(n2),q(n3),q(n4), &
         & q(nt1),lent1,ica,.false.)
  end do
  call add_mpptim('KEXTA',0,1)
  call corlsr(n2)
  !
  ! basis integrals of type (ab/ab)
  !
  if(nsk.lt.2) go to 500
  call add_mpptim('KEXTB',0,0)
  n2=iq2
  nlb=icorrm()
  intyp=nsk
  do iska=2,nsk
    iske=iska-1
    do iskb=1,iske
      intyp=intyp+1
      isy=mult(iska,iskb)
      inta=nt(iska)
      intb=nt(iskb)
      if(inta.eq.0.or.intb.eq.0) cycle
      nd2=intb*(intb+1)/2
      n5=icorr(nd2)
      n6=icorr(nd2)
      n7=icorr(nd2)
      n8=icorr(nd2)
      intbb=ntb(iskb)
      ip1a=ipaa(isy)
      ip1e=ipee(isy)
      ip2a=ipaa(1)
      ip2e=ipee(1)
      nlb=icorrm()
      ippp=max(nums(1),numt(1))
      mblk=min(inta,nlb/(2*nd2+ippp+1))
      if(mblk.eq.0.and.inta.gt.0) call error('address error mblka', &
           & 'kext')
      indx=icori(mblk)
      n3=icorr(mblk*nd2)
      n4=icorr(mblk*nd2)
      n9=icorr(mblk*ippp)
      iadmax=max(iadmax,n9+mblk*ippp)
      if(n9.le.ikee+1) call error('CHAOS B','kext')
      call kextb(q(icaa),q(ikaa),q(n1),q(n2),q(n3),q(n4),q(nt1),lent1, &
           & q(n5),q(n6),q(n7),q(n8),q(n9),q(n9), &
           & isy,nd2,mblk,ica,iq(indx))
      ! start debug
      ! ;      call outvec(q(itmpk),lenn,'global sum kaa before kextb')
      ! ;      call outvec(q(ikaa),lenn,'kaa after kextb')
      ! ;      itmpk1=icorr(lenn)
      ! ;      call fmove(q(ikaa),q(itmpk1),lenn)
      ! ;      call global_sum(q(itmpk1),lenn,'+')
      ! ;      call outvec(q(itmpk),lenn,'global sum kaa before kextb')
      ! ;      call outvec(q(itmpk1),lenn,'global sum kaa after kextb')
      ! ;      do k=1,lenn
      ! ;       q(itmpk1-1+k)=q(itmpk1-1+k)-q(itmpk-1+k)
      ! ;      end do
      ! ;      call outvec(q(itmpk1),lenn,'contribution of kextb')
      ! ;      call corlsr(itmpk)
      ! end
      call corlsr(n5)
    end do
  end do
  call add_mpptim('KEXTB',0,1)
  !
  ! basis integrals of type (s,xy/x,y) etc.
  !
  if(nsk.lt.4) go to 500
  call add_mpptim('KEXTC',0,0)
  ! .....n7,n8,n9: buffers for integral matrix blocks
  call corlsr(n5)
  nlb=n6-ikee-1
  do iska=4,nsk
    do iskb=3,iska-1
      iskab=mult(iska,iskb)
      do iskc=2,iskb-1
        iskd=mult(iskab,iskc)
        if(iskd.ge.iskc) cycle
        intyp=intyp+1
        if(nt(iska).eq.0) cycle
        if(nt(iskb).eq.0) cycle
        if(nt(iskc).eq.0) cycle
        if(nt(iskd).eq.0) cycle
        iskac=mult(iska,iskc)
        iskad=mult(iska,iskd)
        ip1a=ipaa(iskab)
        ip1e=ipee(iskab)
        ip2a=ipaa(iskac)
        ip2e=ipee(iskac)
        ip3a=ipaa(iskad)
        ip3e=ipee(iskad)
        nxs=nt(iskc)*nt(iskd)
        n5=icorr(nxs)
        n6=icorr(nxs)
        nlb=icorrm()
        npmax=ip1e-ip1a+1
        mblk=min(nt(iskb),nlb/(npmax+2*nxs+1))
        if(mblk.le.0.and.nt(iskb).gt.0) call error('address error mblkb', &
             & 'kext')
        indx=icori(mblk)
        npblk=icorr(mblk*npmax)
        n3=icorr(mblk*nxs)
        n4=icorr(mblk*nxs)
        iadmax=max(iadmax,n4+mblk*nxs)
        call kextc(q(icaa),q(ikaa),q(n1),q(n2),q(n3),q(n4),q(nt1),lent1, &
             & q(n5),q(n6),ica,nxs,mblk,iq(indx),q(npblk),npmax)
        call corlsr(n5)
        ! start debug
        ! ;      call outvec(q(ikaa),lenn,'kaa after kextc')
        ! ;      call outvec(q(itmpk),lenn,'global sum kaa before kextc')
        ! ;      call outvec(q(ikaa),lenn,'kaa after kextc')
        ! ;      itmpk1=icorr(lenn)
        ! ;      call fmove(q(ikaa),q(itmpk1),lenn)
        ! ;      call global_sum(q(itmpk1),lenn,'+')
        ! ;      call outvec(q(itmpk),lenn,'global sum kaa before kextc')
        ! ;      call outvec(q(itmpk1),lenn,'global sum kaa after kextc')
        ! ;      do k=1,lenn
        ! ;       q(itmpk1-1+k)=q(itmpk1-1+k)-q(itmpk-1+k)
        ! ;      end do
        ! ;      call outvec(q(itmpk1),lenn,'contribution of kextc')
        ! ;      call corlsr(itmpk)
        ! end
      end do
    end do
  end do
  call add_mpptim('KEXTC',0,1)
  !
  ! .....transform operators back to mo basis; store in final order
  !
  ! .....read orbitals again since overwritten by integrals
500 continue
  if (nprocs.gt.1) then
    if(iprint(39).ge.1) then
      call add_mpptim('KEXT(UNSYNC)',0,1)
      call global_sync(500)
      call add_mpptim('KEXT(SYNC)',0,1)
    end if
    call global_sum(q(ikaa),lenn,'+')
    if(iprint(39).ge.1) call add_mpptim('KEXT(TOT)',0,1)
    if (nprocs.gt.1) mpp_state=1
  end if
end subroutine kext


subroutine kexta(a,b,tint,amn,bmn,ail,bil,t,lent,ica,sparse)
  use common_tapes
  use common_openmp
  use common_cthr
  use common_cprint
  use common_clseg
  use common_cmpp
  use common_cjadr
  implicit double precision (a-h,o-z)
  !
  ! .....amn and bmn may overlap
  !
  logical nosing,notrip,sparse,para
  dimension a(*),b(*),tint(*),t(*),amn(*),bmn(*),ail(*),bil(*)
  dimension ica(8,8,2)
  !
  ! evaluates operators b = k(a) for integral type (aa/aa)
  !
  nskip=0
  ipa=ip1a-1
  ipe=ip1e-ipa
  if(ipe.le.0) return
  if(iprint(39).gt.1) t1=second_x()
  ncpux=1
  ! uncomment next line when the multithread stuff below actually works!
  ! s/cc!\$/!\$/
  ! c!$    ncpux=omp_get_max_threads()
  iaa=ica(iska,1,1)
  ia=iaa-1
  nop=0
  op=0d0
  ita=ispa(1)-ipa
  ise=ita-1
  notrip=ita.gt.ipe
  nosing=ise.le.0
  call fzero(tint,inta**2)
  call get_inpf('KEXTA','THRESH',thr)
  do i = 1,inta
    ma = iaa
    ll=(i-1)*i
    nchp=(ipe-1)/ncpux+1
    para=ncpux.gt.1.and.4*nchp*i*i.gt.minbr1
    ncpu=1
    if(para) ncpu=ncpux
    nchp=(ipe-1)/ncpu+1
    !
    ! c!$omp parallel do if(para),LOCAL(icpu,ip1,ip2,np,ila,ip,n)
    do icpu=1,ncpu
      ip1=(icpu-1)*nchp+1
      ip2=min(ipe,ip1+nchp-1)
      np=(ip1-1)*i
      ila=ia+(ip1-1)*nlga
      do ip=ip1,ip2
        do n=1,i
          bil(np+n)=0d0
          ail(np+n)=a(ila+n)
        end do
        ila=ila+nlga
        np=np+i
      end do
    end do
    !
    do m = 1,i
      !
      ! .....expand integral block
      !
      call ao_integral_matrix_get(i,m,iska,iska,iska,t,lent,ii)
      if (ii.gt.0) then
        if(thr.gt.0d0) then
          len=i*(i-1)/2+m
          do ln=0,len-1
            if(abs(t(ii+ln)).gt.thr) go to 110
          end do
          nskip=nskip+1
          go to 200
        end if
110     ln=0
        nl=1-i
        do l = 1,i
          tint(i*l)=0d0
          tint(ll+l)=0d0
          ne = l
          if(i.eq.l) ne = m
          !
          !
          do n = 1,ne
            tint(ln+n)=t(ii)
            tint(nl+n*i)=t(ii)
            ii=ii+1
          end do
          !
          tint(ln+l)=tint(ln+l)*2d0
          nop=nop+2*ne
          ln=ln+i
          nl=nl+1
        end do
        !
        ! c!$omp parallel do  LOCAL(icpu,ip1,ip2,nch,np,mna,imn,ip,n,mn) shared(i)
        do icpu=1,ncpu
          ip1=(icpu-1)*nchp+1
          ip2=min(ipe,ip1+nchp-1)
          nch=ip2-ip1+1
          np=(ip1-1)*i
          imn=np+1
          mna=ma-1+(ip1-1)*nlga
          do ip=ip1,ip2
            do n=1,m
              amn(np+n)=a(mna+n)
            end do
            mn=mna+2*m
            if(ip.le.ise) then
              do n=m+1,i
                amn(np+n)=a(mn)
                mn=mn+n
              end do
            else
              do n=m+1,i
                amn(np+n)=-a(mn)
                mn=mn+n
              end do
            end if
            mna=mna+nlga
            np=np+i
          end do
          !
          if(sparse.and.i.gt.50) then
            call mxmbs(tint,1,i,amn(imn),1,i,bil(imn),1,i,i,i,nch)
            call mxmas(tint,1,i,ail(imn),1,i,bmn(imn),1,i,i,i,nch)
          else
            call mxmb(tint,i,1,amn(imn),1,i,bil(imn),1,i,i,i,nch)
            call mxma(tint,i,1,ail(imn),1,i,bmn(imn),1,i,i,i,nch)
          end if
          ! start debug
          ! ;      write (6,*) 'imn=',imn
          ! ;      call outsqr (ail(imn),i,i,nch,'ail(imn)')
          ! ;      call outsqr (amn(imn),i,i,nch,'amn(imn)')
          ! ;      call outsqr (bil(imn),i,i,nch,'bil(imn)')
          ! ;      call outsqr (bmn(imn),i,i,nch,'bmn(imn)')
          ! end
        end do
        nop=nop+(4*i*i+2*i)*ipe
        if(nop.gt.100000000) then
          op=op+dble(nop)
          nop=0
        end if
        !
        ! c!$omp parallel do  LOCAL(icpu,ip1,ip2,np,mna,ip,n,mn)
        do icpu=1,ncpu
          ip1=(icpu-1)*nchp+1
          ip2=min(ipe,ip1+nchp-1)
          np=(ip1-1)*i
          mna=ma-1+(ip1-1)*nlga
          do ip=ip1,ip2
            do n=1,m
              b(mna+n)=b(mna+n)+bmn(np+n)
            end do
            mn=mna+2*m
            if(ip.le.ise) then
              do n=m+1,i
                b(mn)=b(mn)+bmn(np+n)
                mn=mn+n
              end do
            else
              do n=m+1,i
                b(mn)=b(mn)-bmn(np+n)
                mn=mn+n
              end do
            end if
            mna=mna+nlga
            np=np+i
          end do
        end do
        !
      end if ! ii
200   ma = ma + m
    end do
    !
    ! c!$omp parallel do LOCAL(icpu,ip1,ip2,np,ila,ip,n)
    do icpu=1,ncpu
      ip1=(icpu-1)*nchp+1
      ip2=min(ipe,ip1+nchp-1)
      np=(ip1-1)*i
      ila=ia+(ip1-1)*nlga
      do ip=ip1,ip2
        do n=1,i
          b(ila+n)=b(ila+n)+bil(np+n)
        end do
        ila=ila+nlga
        np=np+i
      end do
    end do
    nop=nop+2*ipe*i
    !
    ia = ia + i
  end do
  if(nskip.gt.0) write(iout,*) 'thr=',thr,'  nskip=',nskip
  if(iprint(39).gt.1) then
    t1=second_x()-t1
    if(t1.gt.0d0) then
      op=op+dble(nop)
      op=1.d-6*op
      flop=op/t1
      write(6,111) t1,op,flop,inta,ipe
111   format(' CPU FOR KEXTA: ',t20,f9.2, &
           & ' SEC (',f8.2,' MOPS ',f8.2,' MFLOPS), INTA=',i4,'  IPE=',i4)
      call flush6
    end if
  end if
  if(notrip) go to 1000
  i1 = ica(iska,1,2)-1
  do ip = ita,ipe
    ii = i1
    do i = 1,inta
      ii = ii + i
      b(ii) = 0.0
    end do
    i1 = i1 + nlga
  end do
1000 return
end subroutine kexta


subroutine kextb(a,b,tint,tint1,tint2,tint3,t,lent,t1,t2,t3,t4, &
     & aim,bim,isy,nd2,mblk,ica,indx)
!  use common_openmp
  use common_cthr
  use common_cprint
  use common_clseg
  use common_cmpp
  use common_cjadr
  implicit double precision (a-h,o-z)
  !
  ! ...  aim, bim may overlap
  !
  logical nt01,nosng1,nosng2,notrp1,notrp2
  logical skipped, holed
  dimension a(*),b(*),t(*),tint(*),tint1(nd2), &
       & tint2(nd2,mblk),tint3(nd2,mblk),t1(nd2),t2(nd2), &
       & t3(nd2),t4(nd2),ica(8,8,2),aim(*),bim(*)
  dimension indx(mblk)
  !
  ! evaluates operators b=k(a) for integral type (ab/ab)
  !
  minbr=2147483647
  ! $    minbr = minbr1 * omp_get_max_threads()
  thr=thrint
  ipa1=ip1a-1
  ipe1=ip1e-ipa1
  ita1=ispa(isy)-ipa1
  ise1=ita1-1
  nt01=ipe1.le.0
  nosng1=ise1.le.0
  notrp1=ita1.gt.ipe1
  ipa2=ip2a-1
  ipe2=ip2e-ipa2
  if(nt01.and.ipe2.le.0) return
  op=0d0
  nop=0
  if(iprint(39).gt.1) cpu=second_x()
  nd1=inta*intb
  nd3=inta*(inta+1)/2
  ita2=ispa(1)-ipa2
  ise2=ita2-1
  ite1=ipe1-ita1+1
  ite2=ipe2-ita2+1
  nd23=3*nd2
  nosng2=ise2.le.0
  notrp2=ite2.le.0
  ims=ica(iska,1,1)
  imt=ica(iska,1,2)
  lns=ica(iskb,1,1)
  lnt=ica(iskb,1,2)
  ins=ica(iska,isy,1)
  int=ica(iska,isy,2)
  ims0 = ims-1
  imt0 = imt-1
  do i = 1,inta
    mls = ica(iska,isy,1)
    mlt = ica(iska,isy,2)
    !
    !
    ! .....tint1(ln) = (im/ln)
    ! .....tint2(ln,m) = (il/mn)+(in/lm)
    ! .....tint3(ln,m) = (il/mn)-(in/lm)
    !
    mm=0
    skipped = .false.
    holed = nprocs.ne.1
    do m=1,i
      ln=nd2
      call ao_integral_matrix_get(i,m,iska,iska,iskb,t,lent,ii)
      skipped = skipped .or. ii.le.0
      if (ii.le.0) go to 200
      holed = holed .or. skipped
      mm=mm+1
      indx(mm)=m
      call dcopy_X(nd2,t(ii),3,tint1(1),1)
      call dcopy_X(nd2,t(ii+1),3,tint2(1,mm),1)
      call dcopy_X(nd2,t(ii+2),3,tint3(1,mm),1)
      !
      if(nt01) go to 200
      nop=nop+intb*intb*ipe1
      if(notrp1) then
        do ln=1,nd2
          t1(ln)=tint1(ln)+0.5d0*(tint2(ln,mm)+tint3(ln,mm))
        end do
        do ln=1,nd2
          ! 2/06 BUG FIXED PJK was tint3(ln,m)
          t2(ln)=t1(ln)-tint3(ln,mm)
        end do
        nop=nop+4*nd2
        call expdi(t1,t2,tint,intb,intbb)
      else if(nosng1) then
        do ln=1,nd2
          t3(ln)=tint1(ln)-0.5d0*(tint2(ln,mm)+tint3(ln,mm))
        end do
        do ln=1,nd2
          t4(ln)=t3(ln)+tint3(ln,mm)
        end do
        nop=nop+4*nd2
        call expdi(t3,t4,tint1,intb,intbb)
      else
        do ln=1,nd2
          t3ln=0.5d0*(tint2(ln,mm)+tint3(ln,mm))
          t4ln=t3ln-tint3(ln,mm)
          t1(ln)=tint1(ln)+t3ln
          t3(ln)=tint1(ln)-t3ln
          t2(ln)=tint1(ln)+t4ln
          t4(ln)=tint1(ln)-t4ln
        end do
        nop=nop+7*nd2
        call expdi(t1,t2,tint,intb,intbb)
        call expdi(t3,t4,tint1,intb,intbb)
      end if
      !
      if(nosng1) go to 490
      ! $    if(intb*intb*ise1.gt.minbr) then
      ! $    call mxmbm (tint,1,intbb,a(mls),1,intb,b(ins),1,intb,
      ! $   1           intb,intb,ise1)
      ! $    call mxmbm (tint,intbb,1,a(ins),1,intb,b(mls),1,intb,
      ! $   1           intb,intb,ise1)
      ! $    else
      call mxmb (tint,1,intbb,a(mls),1,intb,b(ins),1,intb, &
           & intb,intb,ise1)
      call mxmb (tint,intbb,1,a(ins),1,intb,b(mls),1,intb, &
           & intb,intb,ise1)
      ! $    endif
490   continue
      !
      if(notrp1) go to 200
      !
      ! $    if(intb*intb*ite1.gt.minbr) then
      ! $    call mxmbm(tint1,1,intbb,a(mlt),1,intb,b(int),1,intb,
      ! $   1          intb,intb,ite1)
      ! $    call mxmbm(tint1,intbb,1,a(int),1,intb,b(mlt),1,intb,
      ! $   1          intb,intb,ite1)
      ! $    else
      call mxmb(tint1,1,intbb,a(mlt),1,intb,b(int),1,intb, &
           & intb,intb,ite1)
      call mxmb(tint1,intbb,1,a(int),1,intb,b(mlt),1,intb, &
           & intb,intb,ite1)
      ! $    endif
      !
200   continue
      mls=mls+intb*ise1
      mlt=mlt+intb*ite1
      if (m.eq.i .or. mm.eq.mblk) then
        if (mm.eq.0) go to 210
        !
        ! .....process block of integrals
        !
        nop=nop+nd2*mm*ipe2
        if(nop.gt.100000000) then
          op=op+dble(nop)
          nop=0
        end if
        ! start debug
        ! ;      write (6,*) 'm block: i,m,mm,mblk ',i,m,mm,mblk
        ! ;      call outive (indx,mm,'indx')
        ! ;      call outsqr (tint2,nd2,nd2,mm,'tint2')
        ! ;      call outsqr (tint3,nd3,nd3,mm,'tint3')
        ! end
        !
        if(nosng2) go to 265
        maa=0
        ima=ims0
        do ip=1,ise2
          do mmm=1,mm
            aim(maa+mmm)=a(ima+indx(mmm))
          end do
          maa=maa+mm
          ima=ima+nd3
        end do
        ! start debug
        ! ;      call outsqr (aim,mm,mm,ise2,'aim')
        ! end
        call mxmb(tint2,1,nd2,aim,1,mm,b(lns),1,nd2,nd2,mm,ise2)
        call mxma(tint2,nd2,1,a(lns),1,nd2,bim,1,mm,mm,nd2,ise2)
        ! start debug
        ! ;      call outsqr (bim,mm,mm,ise2,'bim')
        ! end
        maa=0
        ima=ims0
        do ip=1,ise2
          do mmm=1,mm
            b(ima+indx(mmm))=b(ima+indx(mmm))+bim(maa+mmm)
          end do
          maa=maa+mm
          ima=ima+nd3
        end do
265     if(notrp2) go to 210
        maa=0
        ima=imt0
        do ip=1,ite2
          do mmm=1,mm
            aim(maa+mmm)=a(ima+indx(mmm))
          end do
          maa=maa+mm
          ima=ima+nd3
        end do
        call mxmb(tint3,1,nd2,aim,1,mm,b(lnt),1,nd2,nd2,mm,ite2)
        call mxma(tint3,nd2,1,a(lnt),1,nd2,bim,1,mm,mm,nd2,ite2)
        maa=0
        ima=imt0
        do ip=1,ite2
          do mmm=1,mm
            b(ima+indx(mmm))=b(ima+indx(mmm))+bim(maa+mmm)
          end do
          maa=maa+mm
          ima=ima+nd3
        end do
        !
210     continue
        ! ... reset the block
        mm = 0
        skipped = .false.
        holed = nprocs.ne.1
      end if
    end do !m
    ims0 = ims0 + i
    imt0 = imt0 + i
    !
    ins=ins+intb*ise1
    int=int+intb*ite1
  end do
  if(iprint(39).gt.1) then
    cpu=second_x()-cpu
    if(cpu.gt.0d0) then
      op=op+dble(nop)
      op=(4.d-6)*op
      flop=op/cpu
      write(6,111) cpu,op,flop,inta,intb,mblk,ipe1,ipe2
111   format(' CPU FOR KEXTB: ',t20,f9.2, &
           & ' SEC (',f8.2,' MOPS ',f8.2,' MFLOPS), INTA=',i4,' INTB=',i4, &
           & ' MBLK=',i4,' IPE1=',i4,'  IPE2=',i4)
      call flush6
    end if
  end if
end subroutine kextb


subroutine kextc(a,b,tint,tint1,tint2,tint3,t,lent,t2,t3, &
     & ica,nxs,mblk,indx,pblk,npblk)
  use common_tapes
!  use common_openmp
  use common_cthr
  use common_cprint
  use common_corb
  use common_clseg
  use common_cbas
  use common_cmpp
  use common_cjadr
  implicit double precision (a-h,o-z)
  ! .....real version
  logical nt01,nt02,nt03,nosng1,nosng2,nosng3,notrp1,notrp2,notrp3
  logical nt023
  logical skipped, holed
  dimension a(*),b(*),tint(*),tint1(*),tint2(nxs,mblk), &
       & tint3(nxs,mblk),t(*),ica(8,8,2),t2(*),t3(*), &
       & indx(mblk),pblk(mblk,npblk)
  !
  ! evaluates operators b=k(a) for integral type (ab/cd)
  !
  minbr=2147483647
  ! $    minbr = minbr1 * omp_get_max_threads()
  thr=thrint
  nt01=ip1a.gt.ip1e
  nt02=ip2a.gt.ip2e
  nt03=ip3a.gt.ip3e
  nt023=nt02.and.nt03
  if(nt01.and.nt023) return
  if(iprint(39).gt.1) cpu=second_x()
  nop=0
  op=0d0
  ipa1=ip1a-1
  ipe1=ip1e-ipa1
  ita1=ispa(iskab)-ipa1
  ise1=ita1-1
  ipa2=ip2a-1
  ipe2=ip2e-ipa2
  ita2=ispa(iskac)-ipa2
  ise2=ita2-1
  ipa3=ip3a-1
  ipe3=ip3e-ipa3
  ita3=ispa(iskad)-ipa3
  ise3=ita3-1
  ite1=ipe1-ita1+1
  ite2=ipe2-ita2+1
  ite3=ipe3-ita3+1
  notrp1=ite1.le.0
  notrp2=ite2.le.0
  notrp3=ite3.le.0
  nosng1=ise1.le.0
  nosng2=ise2.le.0
  nosng3=ise3.le.0
  lns=ica(iskc,iskab,1)
  lnt=ica(iskc,iskab,2)
  ils=ica(iska,iskac,1)
  ilt=ica(iska,iskac,2)
  ins=ica(iska,iskad,1)
  int=ica(iska,iskad,2)
  ns=nt(iskd)
  nx=nt(iskc)
  ny=nt(iskb)
  nz=nt(iska)
  ! nxs=nx*ns
  nxs3=3*nxs
  if (ise1.gt.npblk .or. ite1.gt.npblk) then
    write (iout,*) 'programming error, npblk, kextc '
    call fehler
  end if
  do i=1,nz
    mls=ica(iskb,iskad,1)
    mlt=ica(iskb,iskad,2)
    mns=ica(iskb,iskac,1)
    mnt=ica(iskb,iskac,2)
    ! start debug
    ! ;      call outive (ica,128,'ica')
    ! ;      write (6,*) 'initial mlt,mns,mnt,mls ',m,mlt,mns,mnt,mls
    ! end
    ims=ica(iska,iskab,1)+(i-1)*ise1*ny
    imt=ica(iska,iskab,2)+(i-1)*ite1*ny
    ! me=0
    mm=0
    skipped = .false.
    holed = .false.
    ims0=ims-1
    imt0=imt-1
    do m=1,ny
      ! do 210 mb = 1,ny,mblk
      ! ma=me+1
      ! me=min(ny,me+mblk)
      ! mm=me-ma+1
      !
      ! .....tint1(ln) = (im/ln)
      ! .....tint2(ln,m) = (il/mn)+(in/lm)
      ! .....tint3(ln,m) = (il/mn)-(in/lm)
      ! call fzero (tint2,nxs*mblk)
      ! call fzero (tint3,nxs*mblk)
      !
      ! do m=1,mm
      call ao_integral_matrix_get(i,m,iska,iskb,iskc,t,lent,ii)
      ! start debug
      ! ;      write (6,*) 'kextc i,m,ii ',i,m,ii
      ! ;      call global_outsqr(t(ii),3,3,nxs,'integral matrix','+')
      ! end
      skipped = skipped .or. ii.le.0
      if (ii.gt.0) then
        holed = holed .or. skipped
        mm=mm+1
        indx(mm)=m
        call dcopy_X(nxs,t(ii),3,tint1(1),1)
        call dcopy_X(nxs,t(ii+1),3,tint2(1,mm),1)
        call dcopy_X(nxs,t(ii+2),3,tint3(1,mm),1)
        !
        if(nt023) go to 200
        do ln=1,nxs
          t2(ln)=0.5d0*(tint2(ln,mm)+tint3(ln,mm))
        end do
        do ln=1,nxs
          t3(ln)=t2(ln)-tint3(ln,mm)
        end do
        nop=nop+3*nxs
        !
        if(nt02) go to 500
        if(nosng2) go to 480
        !
        do ln=1,nxs
          tint(ln)=tint1(ln)+t3(ln)
        end do
        nop=nop+nxs*(ise2+1)
        ! $     if(ns*nx*ise2.gt.minbr) then
        ! $     call mxmbm(tint,ns,1,a(mns),1,ns,b(ils),1,nx,nx,ns,ise2)
        ! $     call mxmbm(tint,1,ns,a(ils),1,nx,b(mns),1,ns,ns,nx,ise2)
        ! $     else
        call mxmb (tint,ns,1,a(mns),1,ns,b(ils),1,nx,nx,ns,ise2)
        call mxmb (tint,1,ns,a(ils),1,nx,b(mns),1,ns,ns,nx,ise2)
        ! $     endif
480     continue
        !
        if(notrp2) go to 500
        !
        do ln=1,nxs
          tint(ln)=tint1(ln)-t3(ln)
        end do
        nop=nop+nxs*(ite2+1)
        ! $    if(ns*nx*ite2.gt.minbr) then
        ! $    call mxmbm(tint,ns,1,a(mnt),1,ns,b(ilt),1,nx,nx,ns,ite2)
        ! $    call mxmbm(tint,1,ns,a(ilt),1,nx,b(mnt),1,ns,ns,nx,ite2)
        ! $    else
        call mxmb (tint,ns,1,a(mnt),1,ns,b(ilt),1,nx,nx,ns,ite2)
        call mxmb (tint,1,ns,a(ilt),1,nx,b(mnt),1,ns,ns,nx,ite2)
        ! $    end if
500     continue
        !
        if(nt03) go to 200
        if(nosng3) go to 520
        !
        do ln=1,nxs
          tint(ln)=tint1(ln)+t2(ln)
        end do
        nop=nop+nxs*(ise3+1)
        ! $    if(ns*nx*ise3.gt.minbr) then
        ! $    call mxmbm(tint,ns,1,a(ins),1,ns,b(mls),1,nx,nx,ns,ise3)
        ! $    call mxmbm(tint,1,ns,a(mls),1,nx,b(ins),1,ns,ns,nx,ise3)
        ! $    else
        call mxmb (tint,ns,1,a(ins),1,ns,b(mls),1,nx,nx,ns,ise3)
        call mxmb (tint,1,ns,a(mls),1,nx,b(ins),1,ns,ns,nx,ise3)
        ! $    end if
520     continue
        !
        if(notrp3) go to 200
        !
        do ln=1,nxs
          tint(ln)=tint1(ln)-t2(ln)
        end do
        nop=nop+nxs*(ite3+1)
        ! $    if(ns*nx*ite3.gt.minbr) then
        ! $    call mxmbm(tint,ns,1,a(int),1,ns,b(mlt),1,nx,nx,ns,ite3)
        ! $    call mxmbm(tint,1,ns,a(mlt),1,nx,b(int),1,ns,ns,nx,ite3)
        ! $    else
        call mxmb (tint,ns,1,a(int),1,ns,b(mlt),1,nx,nx,ns,ite3)
        call mxmb (tint,1,ns,a(mlt),1,nx,b(int),1,ns,ns,nx,ite3)
        ! $    end if
      end if ! ii
200   mlt=mlt+nx*ite3
      mns=mns+ns*ise2
      mnt=mnt+ns*ite2
      mls=mls+nx*ise3
      ! start debug
      ! ;      write (6,*) 'end loop m,mlt,mns,mnt,mls ',m,mlt,mns,mnt,mls
      ! ;      call global_outsqr(b(mns),ns,ns,ise2,'b(mns)','+')
      ! ;      call global_outsqr(b(mls),nx,nx,ise3,'b(mls)','+')
      ! ;      call global_outsqr(b(mnt),ns,ns,ite2,'b(mnt)','+')
      ! ;      call global_outsqr(b(mlt),nx,nx,ite3,'b(mlt)','+')
      ! ;      call global_outsqr(tint2,nxs,nxs,mm,'tint2','+')
      ! ;      call global_outsqr(tint3,nxs,nxs,mm,'tint3','+')
      ! end
      if (m.eq.ny .or. mm.eq.mblk) then
        ! ...  finish m block
        if(mm.eq.0) go to 210
        nop=nop+nxs*mm*ipe1
        if(nop.gt.100000000) then
          op=op+dble(nop)
          nop=0
        end if
        if(nosng1) go to 430
        iaa=0
        lna=lns
        do ix=1,nx
          naa=0
          do mmm=1,mm
            do n=1,ns
              tint(naa+n)=tint2(iaa+n,mmm)
            end do
            naa=naa+ns
          end do
          iaa=iaa+ns
          !if (holed) then
          if (.true. .or. holed) then
            ! write (6,*) 'holed!',mm,ny
            do iise=1,ise1
              do mmm=1,mm
                pblk(mmm,iise) = a(ims0+(iise-1)*ny+indx(mmm))
              end do
            end do
            ! $     if(ns*mm*ise1.gt.minbr) then
            ! $      call mxmbm(tint,1,ns,pblk,1,mblk,b(lna),1,ns,ns,mm,ise1)
            ! $      call mxmam(tint,ns,1,a(lna),1,ns,pblk,1,mblk,mm,ns,ise1)
            ! $     else
            call mxmb (tint,1,ns,pblk,1,mblk,b(lna),1,ns,ns,mm,ise1)
            call mxma (tint,ns,1,a(lna),1,ns,pblk,1,mblk,mm,ns,ise1)
            ! $     end if
            do iise=1,ise1
              do mmm=1,mm
                b(ims0+(iise-1)*ny+indx(mmm)) = &
                     & b(ims0+(iise-1)*ny+indx(mmm)) + pblk(mmm,iise)
              end do
            end do
          else
            ims=ims0+indx(1)
            ! $     if(ns*mm*ise1.gt.minbr) then
            ! $      call mxmbm(tint,1,ns,a(ims),1,ny,b(lna),1,ns,ns,mm,ise1)
            ! $      call mxmbm(tint,ns,1,a(lna),1,ns,b(ims),1,ny,mm,ns,ise1)
            ! $     else
            call mxmb (tint,1,ns,a(ims),1,ny,b(lna),1,ns,ns,mm,ise1)
            call mxmb (tint,ns,1,a(lna),1,ns,b(ims),1,ny,mm,ns,ise1)
            ! $     end if
          end if
          lna=lna+ns*ise1
        end do
        !
430     if(notrp1) go to 210
        iaa=0
        lna=lnt
        do ix=1,nx
          naa=0
          do mmm=1,mm
            do n=1,ns
              tint(naa+n)=tint3(iaa+n,mmm)
            end do
            naa=naa+ns
          end do
          iaa=iaa+ns
          !if (holed) then
          if (.true. .or. holed) then
            !write (6,*) 'holed!',mm,ny
            do iite=1,ite1
              do mmm=1,mm
                pblk(mmm,iite) = a(imt0+(iite-1)*ny+indx(mmm))
              end do
            end do
            ! $     if(ns*mm*ite1.gt.minbr) then
            ! $      call mxmbm(tint,1,ns,pblk,1,mblk,b(lna),1,ns,ns,mm,ite1)
            ! $      call mxmam(tint,ns,1,a(lna),1,ns,pblk,1,mblk,mm,ns,ite1)
            ! $     else
            call mxmb (tint,1,ns,pblk,1,mblk,b(lna),1,ns,ns,mm,ite1)
            call mxma (tint,ns,1,a(lna),1,ns,pblk,1,mblk,mm,ns,ite1)
            ! $     end if
            do iite=1,ite1
              do mmm=1,mm
                b(imt0+(iite-1)*ny+indx(mmm)) = &
                     & b(imt0+(iite-1)*ny+indx(mmm)) + pblk(mmm,iite)
              end do
            end do
          else
            imt=imt0+indx(1)
            ! $     if(ns*mm*ite1.gt.minbr) then
            ! $      call mxmb (tint,1,ns,a(imt),1,ny,b(lna),1,ns,ns,mm,ite1)
            ! $      call mxmb (tint,ns,1,a(lna),1,ns,b(imt),1,ny,mm,ns,ite1)
            ! $     else
            call mxmb (tint,1,ns,a(imt),1,ny,b(lna),1,ns,ns,mm,ite1)
            call mxmb (tint,ns,1,a(lna),1,ns,b(imt),1,ny,mm,ns,ite1)
            ! $     end if
          end if
          lna=lna+ns*ite1
        end do
210     continue
        ! ... reset the block
        mm = 0
        skipped = .false.
        holed = .false.
        !
      end if
    end do ! m
    !
    ilt=ilt+nx*ite2
    ils=ils+nx*ise2
    int=int+ns*ite3
    ins=ins+ns*ise3
  end do
  if(iprint(39).gt.1) then
    cpu=second_x()-cpu
    if(cpu.gt.0d0) then
      op=op+dble(nop)
      op=(4.d-6)*op
      flop=op/cpu
      write(6,111) cpu,op,flop,ns,nx,ny,nz,mblk,ipe1,ipe2,ipe3
111   format(' CPU FOR KEXTC: ',t20,f9.2, &
           & ' SEC (',f8.2,' MOPS ',f8.2,' MFLOPS), INTA=',i3,' INTB=',i3, &
           & ' INTC=',i3,' INTD=',i3,' MBLK=',i4,' IPE1=',i3,'  IPE2=',i3, &
           & ' IPE3=',i3)
      call flush6
    end if
  end if
end subroutine kextc


subroutine expdi(a,b,c,n,nb)
  implicit double precision (a-h,o-z)
  dimension a(*),b(*),c(nb,n)
  ij=0
  do i=1,n
    do j=1,i
      c(j,i)=a(ij+j)
    end do
    do j=1,i
      c(i,j)=b(ij+j)
    end do
    ij=ij+i
    c(i,i)=a(ij)+b(ij)
  end do
end subroutine expdi


subroutine tranie(e,a,q,w,isym,np,moao)
  use common_openmp
  use common_cpar
  use common_corb
  use common_clseg
  use common_cbas
  ! subr transform int/ext parts of operator matrix a, symmetry isym,np
  ! subr by orbitals q
  ! subr w workspace 2*ntqg
  ! subr if moao.gt.0  e <- q(dag) * a * q   (ao->mo)
  ! subr if moao.le.0  a <- q * e * q(dag)   (mo->ao)
  ! subr routine uses actual dimensions in cbas
  ! subr routine calls excom(2) to force ao dimensions
  implicit double precision (a-h,o-z)
  dimension e(*),a(*),q(*),w(*)
  minbr=2147483647
  ! $    minbr = minbr1 * omp_get_max_threads()
  call excom(2)
  iop = 1 + ntqgsx(isym)
  if (moao.gt.0) then
    ! ...  ao->mo transformation
    ! ...  first expand or move the operator matrix
    if (np.ne.0) then
      call expan (a,w(iop),isym,np,0)
    else
      call fmove (a,w(iop),ntqgsx(isym))
    end if
    ie = 1
    do isyme=1,nskcp
      isymi = mult(isyme,isym)
      ni = ivac(isymi)
      ne = nty(isyme)
      if(ni.eq.0.or.ne.eq.0) go to 100
      nxi = ntx(isymi)
      nxe = ntx(isyme)
      nxib = ntbx(isymi)
      nxeb = ntbx(isyme)
      !
      ! The second dimension of matrices sigma is number of active if keepcl.eq.0
      ! or occupied (active+inactive) if keepcl.ne.0 in symmetry isymi.
      ! For keepcl.eq.0 we must point to the first active MO in symmetry isymi
      ! skipping therefore core and inactive in this symmetry.
      ! For keepcl.ne.0 we must point to the first inactive MO in symmetry isymi
      ! skipping only core.
      ! Variable ipoint counts the number of orbitals we must skip.
      !
      ipoint=icore(isymi)
      if(keepcl.eq.0) ipoint=icore(isymi)+icloss(isymi)
      !
      call fzero (w,ni*nxe)
      ! $    if(nxe*nxi*ni.gt.minbr) then
      ! $      call mxmbm(w(iop+ntqsx(isym,isyme)),1,nxeb,
      ! $   >  q(ntqsx(1,isymi)+nxib*ipoint+1),1,nxib,
      ! $   >  w,1,nxe,  nxe,nxi,ni)
      ! $    else
      call mxmb(w(iop+ntqsx(isym,isyme)),1,nxeb, &
           & q(ntqsx(1,isymi)+nxib*ipoint+1),1,nxib, &
           & w,1,nxe,  nxe,nxi,ni)
      ! $    end if
      call fzero (e(ie),ne*ni)
      ! $    if(ne*nxe*ni.gt.minbr) then
      ! $      call mxmbm(q(ntqsx(1,isyme)+nxeb*(ntx(isyme)-nty(isyme))+1),
      ! $   >   nxeb,1,w,1,nxe,  e(ie),1,ne,  ne,nxe,ni)
      ! $    else
      call mxmb (q(ntqsx(1,isyme)+nxeb*(ntx(isyme)-nty(isyme))+1), &
           & nxeb,1,w,1,nxe,  e(ie),1,ne,  ne,nxe,ni)
      ! $    end if
100   ie = ie + ne*ni
    end do
  else
    ! ...   mo->ao transformation
    ie = 1
    call fzero (w(iop),ntqgsx(isym))
    do isyme=1,nskcp
      isymi = mult(isyme,isym)
      ni = ivac(isymi)
      ne = nty(isyme)
      if(ni.eq.0.or.ne.eq.0) go to 200
      nxi = ntx(isymi)
      nxe = ntx(isyme)
      nxib = ntbx(isymi)
      nxeb = ntbx(isyme)
      !
      ! Explanation above ...
      !
      ipoint=icore(isymi)
      if(keepcl.eq.0) ipoint=icore(isymi)+icloss(isymi)
      !
      call fzero (w,ni*nxe)
      ! $    if(nxe*ne*ni.gt.minbr) then
      ! $      call mxmbm(q(ntqsx(1,isyme)+nxeb*(ntx(isyme)-nty(isyme))+1),
      ! $   >  1,nxeb,e(ie),1,ne,  w,1,nxe,  nxe,ne,ni)
      ! $    else
      call mxmb (q(ntqsx(1,isyme)+nxeb*(ntx(isyme)-nty(isyme))+1), &
           & 1,nxeb,e(ie),1,ne,  w,1,nxe,  nxe,ne,ni)
      ! $    end if
      ! $    if(nxe*ni*nxi.gt.minbr) then
      ! $      call mxmbm(w,1,nxe,
      ! $   >  q(ntqsx(1,isymi)+nxib*ipoint+1),nxib,1,
      ! $   >  w(iop+ntqsx(isym,isyme)),1,nxeb,  nxe,ni,nxi)
      ! $    else
      call mxmb (w,1,nxe, &
           & q(ntqsx(1,isymi)+nxib*ipoint+1),nxib,1, &
           & w(iop+ntqsx(isym,isyme)),1,nxeb,  nxe,ni,nxi)
      ! $    end if
200   ie = ie + ne*ni
    end do
    if (np.ne.0) then
      call reduc (w(iop),a,isym,np,1.0d0,1)
    else
      call fmove (w(iop),a,ntqgsx(isym))
    end if
  end if
end subroutine tranie


subroutine kciadr(npaa,npee,npaa1,npee1,icase,nstat1,nmat,ica, &
     & lenn,icpup)
  use common_corb
  use common_cdirect
  use common_cbas
  use cpair
  implicit double precision (a-h,o-z)
  dimension ica(8,8,2),icpup(*)
  !
  ! .....computes symmetry block addresses for cikext
  !
  dimension num(8,2)
  call izero(ica,128)
  call izero(num,16)
  if(icase.ne.0) then
    do ip=npaa,npee
      i=ipair(1,ip)
      j=ipair(2,ip)
      np=ipair(3,ip)
      isy=ipair(4,ip)
      isp=1
      if(np.lt.0) isp=2
      num(isy,isp)=num(isy,isp)+nstat1*nmat
      if(isy.eq.1) then
        do is=1,nsk
          ica(is,isy,isp)=ica(is,isy,isp)+nstat1*nmat*nt(is)*(nt(is)+1)/2
        end do
      else
        do is=1,nsk
          js=mult(is,isy)
          if(js.gt.is) cycle
          ica(is,isy,isp)=ica(is,isy,isp)+nstat1*nmat*nt(is)*nt(js)
        end do
      end if
    end do
  end if
  !
  if(icase.eq.0.or.icase.eq.3) then
    do ip=npaa1,npee1
      i=ipair1(1,ip)
      j=ipair1(2,ip)
      np=ipair1(3,ip)
      isy=ipair1(4,ip)
      isp=1
      if(np.lt.0) isp=2
      num(isy,isp)=num(isy,isp)+nstat1
      if(isy.eq.1) then
        do is=1,nsk
          ica(is,isy,isp)=ica(is,isy,isp)+nstat1*nt(is)*(nt(is)+1)/2
        end do
      else
        do is=1,nsk
          js=mult(is,isy)
          if(js.gt.is) cycle
          ica(is,isy,isp)=ica(is,isy,isp)+nstat1*nt(is)*nt(js)
        end do
      end if
    end do
  end if
  if(direct) then
    len=1
    do isy=1,nskcp
      do isp=1,2
        ica(1,isy,isp)=len
        do is=1,nsk
          js=mult(isy,is)
          if(is.ge.js) ica(is,isy,isp)=len+ntds(isy,is)
        end do
        len=len+num(isy,isp)*ntdgs(isy)
      end do
    end do
    lenn=len-1
    return
  end if
  !
  ia=1
  do isy=1,nskcp
    do is=1,nsk
      js=mult(is,isy)
      if(js.gt.is) cycle
      do isp=1,2
        idum=ica(is,isy,isp)
        ica(is,isy,isp)=ia
        ia=ia+idum
      end do
    end do
  end do
  lenn=ia-1
end subroutine kciadr


subroutine global_outsqr(z,iz,m,n,name,op)
  use common_big
  use common_cmpp
  implicit double precision (a-h,o-z)
  character(len=*) name,op
  character(len=1024) string
  dimension z(iz,n)
  ! call outsqr (z,iz,m,n,'Local '//
  ! &   name(1:len_trim(name)))
  ii=icorr(m*n)
  jj=ii-1
  do j=1,n
    do i=1,m
      q(jj+i)=z(i,j)
    end do
    jj=jj+iz
  end do
  if (nprocs.gt.1) call global_sum(q(ii),m*n,op)
  write(string,*) 'Global ',op(1:len_trim(op)),': ', &
       & name(1:len_trim(name))
  call outsqr (q(ii),m,m,n,string(1:len_trim(string)))
  call corlsr(ii)
end subroutine global_outsqr
