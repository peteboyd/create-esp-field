      program create_esp
      implicit none
      logical xfrac, yfrac, zfrac 
!-----General parameters
      integer n_dim, n_atoms_max, n_atoms_type_max, voxels_max_1d,
     & MAXK, KMAX, KSQMAX, order_type_max, type_list_max
       parameter(n_dim=3,n_atoms_max=5000,n_atoms_type_max=107,
     & voxels_max_1d=3000000,MAXK=5000,KMAX=7,KSQMAX=49,
     & order_type_max=50,type_list_max=80)

      integer n_grid(n_dim), n_atoms, atom_index(n_atoms_max),
     & NMAX(n_dim)
      common/int_general/n_grid, atom_index, NMAX

      double precision axis_zero(n_dim), axis_vector(n_dim,n_dim),
     & real_box_vector(n_dim,n_dim), 
     & Box_volume, recip_box_vector(n_dim,n_dim), KVEC(MAXK)
      common/dbl_general/axis_zero, axis_vector,
     & real_box_vector, Box_volume, recip_box_vector, KVEC  

!-----Mathematical constants
      double precision pi, TWOPI
      parameter(pi=3.14159265,TWOPI=6.2831853)
      double precision perm, k_esp
c     Units from : http://folk.uio.no/michalj/node72.html
      parameter(perm=0.07957747154594767)
c     in atomic units (hartree*bohr/e)
      parameter(k_esp=1.d0/(4.d0*pi*perm)) 
c     Daniele's chioce to add an aditional space outside of the vdw
c     radii to ignore ESP contribution (in Angstrom)
      double precision DO_FACTOR
      parameter(DO_FACTOR=1.0)
!-----Input section
      integer flag_cutoff, fit_RESP, symm_flag, const_flag, 
     & q_tot(4)
      common/int_input/flag_cutoff, fit_RESP, symm_flag, const_flag, 
     & q_tot
      double precision vdw_fact, R_cutoff
      common/dbl_input/vdw_fact, R_cutoff
      double precision vdw_radii(n_atoms_max)

      character*1, allocatable :: rootname(:)
      character(len=70) file_name
      common/char_input/file_name

!-----Process cube section
      integer counter_grid
      common/int_process/counter_grid

      double precision V_pot_good(voxels_max_1d),
     & grid_pos_good(voxels_max_1d,n_dim), sum_good
      common/dbl_process/V_pot_good, grid_pos_good, sum_good


!-----Process recp space section
      double precision alpha
      common/dbl_reciprocal/alpha

c     WARNING: Setting lenrec too high (eg 256) results in
c     undefined behaviour when performing string manipulation
c     in the cif file.
      integer, parameter :: lenrec = 120
      integer, parameter :: ncif = 205
      real(8), parameter :: angs2bohr = 1.889725989d0
     
      integer, allocatable :: atom_number(:)
      real(8), allocatable :: mass(:),eps(:),sig(:)
      character*2, allocatable :: atom(:)
      character*3, allocatable :: atmname(:)
      character(len=lenrec), allocatable :: loopatom(:)
      character*1 title(lenrec)
      integer i_atom, seed, i_dim, i, j, k, i_neigh_k, n_loop,
     & i_neigh_j, i_neigh_i, i_extra, i_loop, j_dim, loopcount
    
      integer idum, ind 
      double precision, allocatable :: q_part(:)
      double precision, allocatable, dimension(:,:) :: atom_pos
      double precision, allocatable, dimension(:,:) :: atom_pos_frac
      double precision q_part_max, tmpx,tmpy,tmpz,
     & L_box(n_dim), grid_pos(n_dim), atom_pos_tmp(n_dim),
     & delta_dist(n_dim), Phi_sum, Phi_real, Phi_recp, Phi_SIC,
     & Phi_EXCL, adf, bdf, cdf, xdf, ydf, zdf,
     & dist, q_sum, grid_spacing, delta_fdist(n_dim), 
     & grid_pos_frac(n_dim)
      real(8) alen,blen,clen,alph,beta,gamma
      character*1 record(lenrec)
      double precision, allocatable, dimension(:,:,:) :: V_coul
      integer, allocatable, dimension(:,:,:) :: V_flag
      real(8) time,ftime,etime
      character*4 p

      character(len=100) :: arg,ciffile
      save record,mass,atom,atmname,rootname
      save loopatom,q_part,atom_number, V_coul, V_flag
c$$$      seed = 159463
c$$$      q_part_max = 20.0
      R_cutoff = 20.
      alpha = (pi/R_cutoff)**2
    
      if (iargc().ne.2) then
          write(*,*)"Usage: create-esp-field [cif file] [cube spacing]"
          call exit(1)
      end if
c     start timing
      call timchk(0, time)
c     get filename from command line 
      call getarg(1, ciffile)

      call getarg(2, arg)
      grid_spacing=dblstr(arg, 100, idum)
      call scancif(ncif, ciffile, n_atoms, loopcount)
      allocate(atom_pos(n_atoms,n_dim),atom(n_atoms),
     &atom_number(n_atoms))
      allocate(atom_pos_frac(n_atoms,n_dim))
      allocate(q_part(n_atoms),mass(n_atoms),atmname(n_atoms))
      call readcif(ncif,ciffile,n_atoms,loopcount,
     & alen,blen,clen,alph,beta,gamma,title,xfrac,yfrac,zfrac)
      call getroot(ciffile,ind)
c      call process_cube(q_part)
c      call process_box
      call getcell(alen,blen,clen,alph,beta,gamma)
      call VDW_radii_array(atom_number,vdw_radii)
      Phi_SIC = 0.d0
      Phi_EXCL = 0.d0
      do i=1,n_atoms
        if((xfrac).and.(yfrac).and.(zfrac))then
            tmpx = atom_pos_frac(i,1)*real_box_vector(1,1) + 
     &             atom_pos_frac(i,2)*real_box_vector(2,1) +
     &             atom_pos_frac(i,3)*real_box_vector(3,1)
            tmpy = atom_pos_frac(i,1)*real_box_vector(1,2) + 
     &             atom_pos_frac(i,2)*real_box_vector(2,2) +
     &             atom_pos_frac(i,3)*real_box_vector(3,2)
            tmpz = atom_pos_frac(i,1)*real_box_vector(1,3) + 
     &             atom_pos_frac(i,2)*real_box_vector(2,3) +
     &             atom_pos_frac(i,3)*real_box_vector(3,3)
            atom_pos(i,1) = tmpx
            atom_pos(i,2) = tmpy
            atom_pos(i,3) = tmpz
        end if
        ! compute EWALD SIC terms
        !Phi_SIC = Phi_SIC + q_part(i)*q_part(i)
        ! compute EWALD Excluded atom terms
c        do j = 1, i-1
c          ! compute distance in fractional coordinates
c          adf = atom_pos_frac(i,1) - atom_pos_frac(j,1)
c          bdf = atom_pos_frac(i,2) - atom_pos_frac(j,2)
c          cdf = atom_pos_frac(i,3) - atom_pos_frac(j,3)
c
c          ! find minimum image distance
c          adf = adf - dble(nint(adf))
c          bdf = bdf - dble(nint(bdf))
c          cdf = cdf - dble(nint(cdf))
c
c          ! convert to cartesian
c          xdf = adf*real_box_vector(1,1) + 
c     &          bdf*real_box_vector(2,1) +
c     &          cdf*real_box_vector(3,1)
c          ydf = adf*real_box_vector(1,2) + 
c     &          bdf*real_box_vector(2,2) +
c     &          cdf*real_box_vector(3,2)
c          zdf = adf*real_box_vector(1,3) + 
c     &          bdf*real_box_vector(2,3) +
c     &          cdf*real_box_vector(3,3)
c          dist = sqrt(xdf*xdf + ydf*ydf + zdf*zdf)
c
c          call Ewald_excl_term(q_part(i),q_part(j),dist,Phi_EXCL)
c
c        end do
      end do
c      Phi_SIC = Phi_SIC*k_esp
c      write(*,'(a,2f12.6)')
c     &"Self Interaction and Exclusion terms from Ewald sum:",
c     &Phi_SIC,Phi_EXCL
      write(*,'(a,3f12.6)')"Exact spacing of grid points = ",
     &alen/n_grid(1),
     &blen/n_grid(2), 
     &clen/n_grid(3)

c      do i_atom = 1, n_atoms
c        write(*,*)atom_number(i_atom)
c      end do
c      call exit(1)

      call timchk(0,etime)
      do j_dim=1, n_dim
       dist = 0.d0
       do i_dim=1, n_dim
         dist = dist + real_box_vector(j_dim,i_dim)**2
       end do
       dist = sqrt(dist)
       NMAX(j_dim) = int(R_cutoff/dist) + 1
c       write(*,*) "# of cells in",j_dim, "direction: ", NMAX(j_dim) 
      end do

      allocate (V_coul(n_grid(1),n_grid(2),n_grid(3))) 
      allocate (V_flag(n_grid(1),n_grid(2),n_grid(3))) 
      do k=1, n_grid(3)
       do j=1, n_grid(2)
        do i=1, n_grid(1)
         V_coul(i,j,k) = 0.d0
         V_flag(i,j,k) = 1
         if(DO_FACTOR.ne.0)then
           ! flag if not to compute ESP at this grid point
           grid_pos_frac(1) = dble(i-1)/dble(n_grid(1))
           grid_pos_frac(2) = dble(j-1)/dble(n_grid(2))
           grid_pos_frac(3) = dble(k-1)/dble(n_grid(3)) 
           do i_atom=1, n_atoms
             delta_fdist(1) = grid_pos_frac(1)-atom_pos_frac(i_atom, 1)
             delta_fdist(2) = grid_pos_frac(2)-atom_pos_frac(i_atom, 2)
             delta_fdist(3) = grid_pos_frac(3)-atom_pos_frac(i_atom, 3)
             ! shift by pbc
             delta_fdist(1) = delta_fdist(1)-dble(nint(delta_fdist(1)))
             delta_fdist(2) = delta_fdist(2)-dble(nint(delta_fdist(2)))
             delta_fdist(3) = delta_fdist(3)-dble(nint(delta_fdist(3)))

             ! compute cartesian
             dist=0.d0
             do i_dim=1, 3
               delta_dist(i_dim) = 
     &           delta_fdist(1)*real_box_vector(1,i_dim)+
     &           delta_fdist(2)*real_box_vector(2,i_dim)+
     &           delta_fdist(3)*real_box_vector(3,i_dim)
             !  atom_pos_tmp(i_dim) = atom_pos(i_atom,i_dim)
             !  delta_dist(i_dim) = grid_pos(i_dim) - atom_pos_tmp(i_dim)
               dist = dist + delta_dist(i_dim)**2
             end do
             ! check for nearby atoms
             dist = sqrt(dist)
             if(dist.le.(vdw_radii(i_atom)+DO_FACTOR))then
               V_flag(i,j,k)=0
             endif
           end do 
         endif
        end do
       end do
      end do

      call k_vectors_coeff
      do k=1, n_grid(3)
       do j=1, n_grid(2)
        do i=1, n_grid(1)
          if(V_flag(i,j,k).eq.1)then
            do i_dim=1, 3
              grid_pos(i_dim) = ((i-1)*axis_vector(1,i_dim) + 
     &      (j-1)*axis_vector(2,i_dim) + (k-1)*axis_vector(3,i_dim))
            end do
            grid_pos_frac(1) = dble(i-1)/dble(n_grid(1))
            grid_pos_frac(2) = dble(j-1)/dble(n_grid(2))
            grid_pos_frac(3) = dble(k-1)/dble(n_grid(3)) 
            do i_atom=1, n_atoms
              delta_fdist(1) = grid_pos_frac(1)-atom_pos_frac(i_atom,1)
              delta_fdist(2) = grid_pos_frac(2)-atom_pos_frac(i_atom,2)
              delta_fdist(3) = grid_pos_frac(3)-atom_pos_frac(i_atom,3)
              ! shift by pbc
              delta_fdist(1) = delta_fdist(1)-dble(nint(delta_fdist(1)))
              delta_fdist(2) = delta_fdist(2)-dble(nint(delta_fdist(2)))
              delta_fdist(3) = delta_fdist(3)-dble(nint(delta_fdist(3)))

              ! compute cartesian
              dist=0.d0
              do i_dim=1, 3
               delta_dist(i_dim) = 
     &            delta_fdist(1)*real_box_vector(1,i_dim)+
     &            delta_fdist(2)*real_box_vector(2,i_dim)+
     &            delta_fdist(3)*real_box_vector(3,i_dim)
                  atom_pos_tmp(i_dim) = atom_pos(i_atom,i_dim)
              ! delta_dist(i_dim) = grid_pos(i_dim) - atom_pos_tmp(i_dim)
              end do
!--------   Creating a file with the Coulomb potential for testing purposes
              call Ewald_recp_sum(0,1,delta_dist,Phi_recp)
              call Ewald_real_sum(0,1,grid_pos,atom_pos_tmp,Phi_real)
              ! TODO(pboyd): add self interaction correction.
              ! TODO(pboyd): remove interactions between MOF atoms in
              ! RECP sum.

              Phi_sum = q_part(i_atom)*(Phi_real + Phi_recp)*k_esp
              V_coul(i,j,k) = V_coul(i,j,k) + Phi_sum
              !write(*,*)Phi_real,Phi_recp,V_coul(i,j,k)
            enddo
c            V_coul(i,j,k) = V_coul(i,j,k) - Phi_SIC - Phi_EXCL
          endif
        end do 
       end do
      end do
      call timchk(1,etime)

      call timchk(0,ftime)
      call writecube(ind)

      call timchk(1,ftime)
      call timchk(1,time)

      write(*,'(a,f20.5,a)')"ESP generating time :",etime," seconds."
      write(*,'(a,f20.5,a)')"File writing time   :",ftime," seconds."
      write(*,'(a,f20.5,a)')"Total wall time     :",time," seconds."
      contains 
      
      subroutine writecube(ind)
c***********************************************************************
c     
c     routine to write a cube file from the potential generated
c     in the main program
c
c***********************************************************************
      implicit none
      integer i_dim, i_atoms, ind, n_loop, i_extra
      integer i, j, k
      character(len=ind+5) cubefile
      character(len=ind) root

      do i=1, ind
        root(i:i) = rootname(i)
      end do
      cubefile=root//".cube"

      open(10,file=cubefile,status="new")
!----Write cube file
      write(10,*)"Cube file generated with create-esp-field"
      write(10,*)"*****************************************"

!-----Writing the number of atoms and origin
      write(10,'(i5,3f12.6)') n_atoms,
     &      (axis_zero(i_dim), i_dim=1,n_dim)

!-----Writing the voxels arrray
      write(10,'(i5,3f12.6)')n_grid(1), 
     &      (axis_vector(1,i_dim), i_dim=1,3)
      write(10,'(i5,3f12.6)')n_grid(2), 
     &      (axis_vector(2,i_dim), i_dim=1,3)
      write(10,'(i5,3f12.6)')n_grid(3), 
     &      (axis_vector(3,i_dim), i_dim=1,3)

      do i_atom=1, n_atoms
c        write(*,*)atom_number(i_atom),atom_index(i_atom),
c     &atom_pos(i_atom,1:3)
        write(10,'(2i5,3f12.6)')atom_number(i_atom), atom_index(i_atom),
     & (atom_pos(i_atom,i_dim), i_dim=1, 3) 
      end do
      n_loop = floor(real(n_grid(3))/6)
      i_extra = mod(n_grid(3),6)
      write(p,'(i1)') i_extra
c      write(*,*) p, i_extra, n_loop

      do i=1, n_grid(1)
       do j=1, n_grid(2)
        if(i_extra.ne.0) then    
         do i_loop=1, n_loop 
          write(10,'(6e13.5)') (V_coul(i,j,k+(i_loop-1)*6), k=1, 6) 
         end do          
         write(10,'('//trim(p)//'e13.5)') (V_coul(i,j,k+(i_loop-1)*6), 
     &   k=1, i_extra)         
        else
         do i_loop=1, n_loop 
          write(10,'(6e13.5)') (V_coul(i,j,k+(i_loop-1)*6), k=1, 6)
         end do
        end if          
       end do
      end do
      close(10)
      end subroutine writecube

      real(8) function atmmass(i) 

c***********************************************************************
c     
c     function to generate an atom label from an atomic 
c     number 
c
c***********************************************************************
      implicit none
      character*2 i
      if (i.eq.'H ')then
        atmmass=1.00794
      elseif(i.eq.'He')then
        atmmass=4.002602
      elseif(i.eq.'Li')then
        atmmass=6.941
      elseif(i.eq.'Be')then
        atmmass=9.012182
      elseif(i.eq.'B ')then
        atmmass=10.811
      elseif(i.eq.'C ')then
        atmmass=12.0107
      elseif(i.eq.'N ')then
        atmmass=14.00674
      elseif(i.eq.'O ')then
        atmmass=15.9994
      elseif(i.eq.'F ')then
        atmmass=18.9984032
      elseif(i.eq.'Ne')then
        atmmass=20.1797
      elseif(i.eq.'Na')then
        atmmass=22.989770
      elseif(i.eq.'Mg')then
        atmmass=24.3050
      elseif(i.eq.'Al')then
        atmmass=26.981538
      elseif(i.eq.'Si')then
        atmmass=28.0855
      elseif(i.eq.'P ')then
        atmmass=30.973761
      elseif(i.eq.'S ')then
        atmmass=32.066
      elseif(i.eq.'Cl')then
        atmmass=35.4527
      elseif(i.eq.'Ar')then
        atmmass=39.948
      elseif(i.eq.'K ')then
        atmmass=39.0983
      elseif(i.eq.'Ca')then
        atmmass=40.078
      elseif(i.eq.'Sc')then
        atmmass=44.955910
      elseif(i.eq.'Ti')then
        atmmass=47.867
      elseif(i.eq.'V ')then
        atmmass=50.9415
      elseif(i.eq.'Cr')then
        atmmass=51.9961
      elseif(i.eq.'Mn')then
        atmmass=54.938049
      elseif(i.eq.'Fe')then
        atmmass=55.845
      elseif(i.eq.'Co')then
        atmmass=58.9332
      elseif(i.eq.'Ni')then
        atmmass=58.6934
      elseif(i.eq.'Cu')then
        atmmass=63.546
      elseif(i.eq.'Zn')then
        atmmass=65.39
      elseif(i.eq.'Ga')then
        atmmass=69.723
      elseif(i.eq.'Ge')then
        atmmass=72.61
      elseif(i.eq.'As')then
        atmmass=74.9216
      elseif(i.eq.'Se')then
        atmmass=78.96
      elseif(i.eq.'Br')then
        atmmass=79.904
      elseif(i.eq.'Kr')then
        atmmass=83.80
      elseif(i.eq.'Y ')then
        atmmass=88.90585
      elseif(i.eq.'Zr')then
        atmmass=91.224
      elseif(i.eq.'Nb')then
        atmmass=92.90638
      elseif(i.eq.'Mo')then
        atmmass=95.94
      elseif(i.eq.'Ru')then
        atmmass=101.07
      elseif(i.eq.'Rh')then
        atmmass=102.90550
      elseif(i.eq.'Pd')then
        atmmass=106.42
      elseif(i.eq.'Ag')then
        atmmass=107.8682
      elseif(i.eq.'Cd')then
        atmmass=112.411
      elseif(i.eq.'In')then
        atmmass=114.818
      elseif(i.eq.'Sn')then
        atmmass=118.710
      elseif(i.eq.'Sb')then
        atmmass=121.760
      elseif(i.eq.'Te')then
        atmmass=127.760
      elseif(i.eq.'I ')then
        atmmass=126.90447
      elseif(i.eq.'Xe')then
        atmmass=131.29
      elseif(i.eq.'Ba')then
        atmmass=137.327
      elseif(i.eq.'Hf')then
        atmmass=178.49
      elseif(i.eq.'Ta')then
        atmmass=180.9479
      elseif(i.eq.'W ')then
        atmmass=183.84
      elseif(i.eq.'Re')then
        atmmass=186.207
      elseif(i.eq.'Os')then
        atmmass=190.23
      elseif(i.eq.'Ir')then
        atmmass=192.217
      elseif(i.eq.'Pt')then
        atmmass=195.078
      elseif(i.eq.'Au')then
        atmmass=196.96655
      elseif(i.eq.'Hg')then
        atmmass=200.59
      elseif(i.eq.'Tl')then
        atmmass=204.3833
      elseif(i.eq.'Pb')then
        atmmass=207.2
      endif
      return
      end function atmmass

      integer function intstr(word,len,lst)

c***********************************************************************
c     
c     function for extracting integers from a 
c     character string
c     
c     parameters:
c     word   - input character string
c     len    - working length of character string
c     lst    - location of space character at end of
c     integer string
c
c***********************************************************************
      
      implicit none

      logical flag,count,final
      character*1 n,word,ksn
      integer lst,len,j,isn

      dimension n(0:9),word(len)
      data n/'0','1','2','3','4','5','6','7','8','9'/

      isn=1
      lst=0
      ksn='+'
      intstr=0
      flag=.false.
      final=.false.
      count=.false.
      
      do while(lst.lt.len.and.(.not.final))

        lst=lst+1
        flag=.false.

        do j=0,9
          
          if(n(j).eq.word(lst))then
            
            intstr=10*intstr+j
            count=.true.
            flag=.true.
            
          endif
          
        enddo

        if(count.and.(.not.flag))final=.true.
        if(flag.and.ksn.eq.'-')isn=-1
        ksn=word(lst)

      enddo

      intstr=isn*intstr

      do j=lst,len
        word(j-lst+1)=word(j)
      enddo
      do j=len-lst+2,len
        word(j)=' '
      enddo

      return
      end function intstr

      real(8) function dblstr(word,len,lst)

c***********************************************************************
c     
c     Function for extracting double precisions from a 
c     character string. 
c     
c     parameters:
c     word   - input character string
c     len    - working length of character string
c     lst    - location of space character at end of
c     double precision string
c     
c***********************************************************************
      
      implicit none
      
      character*1 n,word,ksn,dot,d,e
      logical flag,ldot,start,final
      integer len,lst,iexp,idum,i,j,fail
      real(8) sn,ten,one
      dimension n(0:9),word(len)
      character*1, allocatable :: work(:)

      data n/'0','1','2','3','4','5','6','7','8','9'/
      data dot/'.'/
      data d/'d'/
      data e/'e'/
      
      allocate(work(len),stat=fail)

      lst=0
      sn=1.d0
      ksn='+'
      ten=10.d0
      one=1.d0
      
      dblstr=0.d0
      iexp=0
      idum=0
      start=.false.
      ldot=.false.
      final=.false.

      do while(lst.lt.len.and.(.not.final))
        
        lst=lst+1
        flag=.false.
        
        do j=0,9
          
          if(n(j).eq.word(lst))then
            
            dblstr=ten*dblstr+one*dble(j)
            flag=.true.
            start=.true.
            
          endif
          
        enddo
        
        if(dot.eq.word(lst))then
          
          flag=.true.
          ten=1.d0
          ldot=.true.
          start=.true.
          
        endif

        if(flag.and.ksn.eq.'-') sn=-1.d0
        if(ldot) one=one/10.d0
        ksn=word(lst)
        if(ksn.eq."D")ksn="d"
        if(ksn.eq."E")ksn="e"
        
        if(start)then
          
          if(d.eq.ksn.or.e.eq.ksn)then
            
            do i=1,len-lst
              work(i)=word(i+lst)
            enddo
            iexp=intstr(work,len-lst,idum)
            final=.true.

          endif

          if(.not.flag)final=.true.
          
        endif
        
      enddo
      
      dblstr=sn*dblstr*(10.d0**iexp)
      lst=lst+idum
      
      do j=lst,len
        word(j-lst+1)=word(j)
      enddo
      do j=len-lst+2,len
        word(j)=' '
      enddo

      deallocate(work,stat=idum)

      return
      end function dblstr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----Subroutine to compute the generalized reciprocal space vectors
      subroutine k_vectors_coeff
      implicit none      

      integer TOTK, KX, KY, KZ, KSQ

      double precision beta, RKX, RKY, RKZ, RKSQ
      
      beta = 1.0/(4*alpha)

      TOTK = 0
      DO KX = 0, KMAX             
       DO KY = -KMAX, KMAX
        DO KZ = -KMAX, KMAX          
         RKX = recip_box_vector(1,1)*REAL(KX) + 
     &   recip_box_vector(2,1)*REAL(KY) +         
     &   recip_box_vector(3,1)*REAL(KZ)

         RKY = recip_box_vector(1,2)*REAL(KX) + 
     &   recip_box_vector(2,2)*REAL(KY) +         
     &   recip_box_vector(3,2)*REAL(KZ)
         RKZ = recip_box_vector(1,3)*REAL(KX) + 
     &   recip_box_vector(2,3)*REAL(KY) +         
     &   recip_box_vector(3,3)*REAL(KZ)

         KSQ = KX * KX + KY * KY + KZ * KZ
         IF ((KSQ.LT.KSQMAX).AND.(KSQ.NE.0)) THEN
           TOTK = TOTK + 1
           IF (TOTK.GT.MAXK) STOP 'KVEC IS TOO SMALL'
           RKSQ = RKX * RKX + RKY * RKY + RKZ * RKZ
           KVEC(TOTK) = (2*TWOPI/Box_volume) * EXP(-beta*RKSQ) / RKSQ
           !PB - DEBUG
           !write(*,*)KVEC(TOTK),RKSQ,KSQ,KX,KY,KZ,RKX,RKY,RKZ
         ENDIF            
        END DO
       END DO
      END DO

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----Subroutine to compute the reciprocal space series within the Ewald
!-----technique
      subroutine Ewald_recp_sum(term_flag,input_dim,delta_r,Phi)
      implicit none
      integer TOTK, KX, KY, KZ, KSQ, term_flag, input_dim

      double precision RX, RY, RZ, FACTOR, Phi, delta_r(n_dim), RKvec

      COMPLEX EIKX(0:KMAX), EIKY(-KMAX:KMAX), EIKZ(-KMAX:KMAX), EIKR

      RX = delta_r(1)
      RY = delta_r(2)
      RZ = delta_r(3)
      EIKX(0) = (1.0, 0.0)
      EIKY(0) = (1.0, 0.0)
      EIKZ(0) = (1.0, 0.0)

      EIKX(1) = CMPLX (COS(recip_box_vector(1,1)*RX + 
     & recip_box_vector(1,2)*RY + recip_box_vector(1,3)*RZ),
     & SIN(recip_box_vector(1,1)*RX + 
     & recip_box_vector(1,2)*RY + recip_box_vector(1,3)*RZ))

      EIKY(1) = CMPLX (COS(recip_box_vector(2,1)*RX + 
     & recip_box_vector(2,2)*RY + recip_box_vector(2,3)*RZ),
     & SIN(recip_box_vector(2,1)*RX + 
     & recip_box_vector(2,2)*RY + recip_box_vector(2,3)*RZ))

      EIKZ(1) = CMPLX (COS(recip_box_vector(3,1)*RX + 
     & recip_box_vector(3,2)*RY + recip_box_vector(3,3)*RZ),
     & SIN(recip_box_vector(3,1)*RX + 
     & recip_box_vector(3,2)*RY + recip_box_vector(3,3)*RZ))
      EIKY(-1) = CONJG (EIKY(1))
      EIKZ(-1) = CONJG (EIKZ(1))

      DO KX = 2, KMAX     
         EIKX(KX) = EIKX(KX-1) * EIKX(1)
      END DO

      DO KY = 2, KMAX          
         EIKY(KY) = EIKY(KY-1) * EIKY(1)
         EIKY(-KY) = CONJG(EIKY(KY))
      END DO

      DO KZ = 2, KMAX
         EIKZ(KZ) = EIKZ(KZ-1) * EIKZ(1)
         EIKZ(-KZ) = CONJG ( EIKZ(KZ) )
      END DO

      Phi = 0.d0
      TOTK = 0
     
      DO KX = 0, KMAX
        IF (KX.EQ.0) THEN
          FACTOR = 1.0
        ELSE
          FACTOR = 2.0
        ENDIF
        DO KY = -KMAX, KMAX
           DO KZ = -KMAX, KMAX
              KSQ = KX * KX + KY * KY + KZ * KZ
              IF ((KSQ.LT.KSQMAX).AND.(KSQ.NE.0)) THEN
                 TOTK = TOTK + 1                                  
                 EIKR = EIKX(KX) * EIKY(KY) * EIKZ(KZ)
                 if(term_flag.eq.0) then                  
                  Phi = Phi + FACTOR * KVEC(TOTK) * real(EIKR)
                 else if(term_flag.eq.1) then
                  RKvec = recip_box_vector(1,input_dim)*REAL(KX) + 
     &            recip_box_vector(2,input_dim)*REAL(KY) +         
     &            recip_box_vector(3,input_dim)*REAL(KZ)
                  Phi = Phi + FACTOR * KVEC(TOTK) * 
     &            real((0.,1.)*EIKR) * RKvec
                 end if
              ENDIF
           END DO
        END DO
      END DO
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----Subroutine to compute the real space series within the Ewald
!-----technique
      subroutine Ewald_real_sum(term_flag,input_dim,pos_grid,
     & pos_charge,Phi)
      implicit none
      integer i_cellx, i_celly, i_cellz, i_dim, term_flag, input_dim 

      double precision pos_grid(n_dim), pos_charge(n_dim), Phi, 
     & delta_dist(n_dim)
      real dist

      Phi = 0.d0
      do i_cellx = -NMAX(1), NMAX(1)
       do i_celly = -NMAX(2), NMAX(2)
        do i_cellz = -NMAX(3), NMAX(3)
         dist = 0.d0
         do i_dim=1, 3
          delta_dist(i_dim) = pos_grid(i_dim) 
     &    - ( pos_charge(i_dim) + 
     &    i_cellx*real_box_vector(1,i_dim) +   
     &    i_celly*real_box_vector(2,i_dim)
     &    + i_cellz*real_box_vector(3,i_dim) )
          dist = dist + delta_dist(i_dim)**2
         end do          
         dist = sqrt(dist)
         if(dist.le.R_cutoff) then
          if(term_flag.eq.0) then
           Phi = Phi + erfc(sqrt(alpha)*dist)/dist
          else if(term_flag.eq.1) then
           Phi = Phi + (delta_dist(input_dim)/dist**2)*
     &    (2.*sqrt(alpha)*exp(-alpha*dist**2)/sqrt(pi) + 
     &     erfc(sqrt(alpha)*dist)/dist)
          end if
         end if
        end do
       end do
      end do
      end subroutine

      subroutine Ewald_excl_term(qi,qj,dist,Phi_EXCL)
      implicit none
      double precision qi,qj,dist,Phi_EXCL
      double precision alpr, alpr2, erfr, chgprd, exp1
      double precision a1, a2, a3, a4, a5, pp
      double precision rr3, r10, r42, r216, tt

      data a1,a2,a3/0.254829592d0,-0.284496736d0,1.421413741d0/
      data a4,a5,pp/-1.453152027d0,1.061405429d0,0.3275911d0/
      data rr3/0.333333333333d0/,r10/0.1d0/,r42/0.02380952381d0/
      data r216/4.62962962963d-3/
      alpr=dist*alpha
      alpr2 = alpr*alpr
      chgprd = qi*qj * k_esp
      if(alpr.lt.1.d-2)then

          erfr = 2.d0*chgprd*(alpha/pi/pi)*
     &(1.d0+alpr2*(-rr3+alpr2*(r10+alpr2*(-r42+alpr2*r216))))

      else
          tt=1.d0/(1.d0 + pp*alpr)
          exp1 = exp(-(alpr)**2)
          erfr = (1.d0-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1)*
     &             chgprd/dist
      endif
      Phi_EXCL = Phi_EXCL + erfr
      return  

      end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----Complementary error function
      double precision function erfc(X)

C    *******************************************************************
C    ** APPROXIMATION TO THE COMPLEMENTARY ERROR FUNCTION             **
C    **                                                               **
C    ** REFERENCE:                                                    **
C    **                                                               **
C    ** ABRAMOWITZ AND STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS,    **
C    **    NATIONAL BUREAU OF STANDARDS, FORMULA 7.1.26               **
C    *******************************************************************

      double precision  A1, A2, A3, A4, A5, P

      PARAMETER ( A1 = 0.254829592, A2 = -0.284496736 )
      PARAMETER ( A3 = 1.421413741, A4 = -1.453152027 )
      PARAMETER ( A5 = 1.061405429, P  =  0.3275911   )
 
      double precision  T, X, XSQ, TP

C    ****************************************************************

      T  = 1.0 / ( 1.0 + P * X )
      XSQ = X * X

      TP = T * ( A1 + T * ( A2 + T * ( A3 + T * ( A4 + T * A5 ) ) ) )

      erfc = TP * exp ( -XSQ )

      return
      end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----Subroutine to deal with the simulation box and atomic variables
!-----according to the PBCs    
      subroutine process_box
      implicit none
      integer i_dim, j_dim, i_atom

      double precision a(n_dim), b(n_dim), c(n_dim), ab(n_dim),
     & Vect_Dot, r(n_dim), cross_tmp(n_dim), signo, dot, dot_1,
     & cross_tmp_1(n_dim)

!-----Creating the real space box
      do j_dim=1, n_dim
       do i_dim=1, n_dim
        real_box_vector(j_dim,i_dim) = 
     &  n_grid(j_dim)*axis_vector(j_dim,i_dim)
       end do
      end do

      a(1) = real_box_vector(1,1)
      a(2) = real_box_vector(1,2)
      a(3) = real_box_vector(1,3)
      b(1) = real_box_vector(2,1)
      b(2) = real_box_vector(2,2)
      b(3) = real_box_vector(2,3)
      c(1) = real_box_vector(3,1)
      c(2) = real_box_vector(3,2)
      c(3) = real_box_vector(3,3)

      call Vect_Cross(b,c,ab)
      Box_Volume = Vect_Dot(a,ab)
      write(*,*) "------PRINTING SOME USEFUL INFO---------"
      write(*,*) "----------------------------------------"
      write(*,*) "----------------------------------------"
      write(*,*) "Real box volume ", Box_volume

!-----Creating the reciprocal space box
      call Vect_Cross(a,b,ab)      
      do i_dim=1, n_dim
       recip_box_vector(3,i_dim) = (2*pi/Box_volume)*ab(i_dim)
      end do
      call Vect_Cross(c,a,ab)      
      do i_dim=1, n_dim
       recip_box_vector(2,i_dim) = (2*pi/Box_volume)*ab(i_dim)
      end do      
      call Vect_Cross(b,c,ab)      
      do i_dim=1, n_dim
       recip_box_vector(1,i_dim) = (2*pi/Box_volume)*ab(i_dim)
      end do       

!-----Setting the origin of the real box to (0,0,0) 
      do i_atom=1, n_atoms
       do i_dim=1, n_dim      
        atom_pos(i_atom,i_dim) = 
     &  atom_pos(i_atom,i_dim) - axis_zero(i_dim) 
        r(i_dim) = atom_pos(i_atom,i_dim)      
       end do
!------Including the atoms within the box (applying pbc's)       
       call Vect_Cross(r,b,cross_tmp) 
       dot = Vect_Dot(c,cross_tmp)
       call Vect_Cross(a,b,cross_tmp_1) 
       dot_1 = Vect_Dot(c,cross_tmp_1)
       signo = dot/dot_1
       if(signo.lt.0) then
        do i_dim=1, n_dim
         atom_pos(i_atom,i_dim) = 
     &   atom_pos(i_atom,i_dim) + a(i_dim)
        end do
       else if(signo.gt.1) then
        do i_dim=1, n_dim
         atom_pos(i_atom,i_dim) = 
     &   atom_pos(i_atom,i_dim) - a(i_dim)
        end do
       end if

       call Vect_Cross(r,c,cross_tmp) 
       dot = Vect_Dot(a,cross_tmp)
       call Vect_Cross(b,c,cross_tmp_1) 
       dot_1 = Vect_Dot(a,cross_tmp_1)
       signo = dot/dot_1 
       if(signo.lt.0) then
        do i_dim=1, n_dim
         atom_pos(i_atom,i_dim) = 
     &   atom_pos(i_atom,i_dim) + b(i_dim)
        end do
       else if(signo.gt.1) then
        do i_dim=1, n_dim
         atom_pos(i_atom,i_dim) = 
     &   atom_pos(i_atom,i_dim) - b(i_dim)
        end do
       end if

       call Vect_Cross(r,a,cross_tmp) 
       dot = Vect_Dot(b,cross_tmp)
       call Vect_Cross(c,a,cross_tmp_1) 
       dot_1 = Vect_Dot(b,cross_tmp_1)
       signo = dot/dot_1
       if(signo.lt.0) then
        do i_dim=1, n_dim
         atom_pos(i_atom,i_dim) = 
     &   atom_pos(i_atom,i_dim) + c(i_dim)
        end do
       else if(signo.gt.1) then
        do i_dim=1, n_dim
         atom_pos(i_atom,i_dim) = 
     &   atom_pos(i_atom,i_dim) - c(i_dim)
        end do
       end if
 
      end do

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----Subroutine to compute the cross product of two vectors
      Subroutine Vect_Cross(v1,v2,v3)
      Implicit NONE

      double precision v1(3),v2(3),v3(3)

      v3(1)=v1(2)*v2(3)-v1(3)*v2(2)
      v3(2)=v1(3)*v2(1)-v1(1)*v2(3)
      v3(3)=v1(1)*v2(2)-v1(2)*v2(1)
      
      End subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----Function to compute the dot product of two vectors
      Function Vect_Dot(a,b)
      Implicit NONE

      double precision a(3),b(3), Vect_Dot

      Vect_Dot=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
      
      End function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----Subroutine to process cube file
      subroutine process_cube(charge)
      implicit none
      integer i, j, k, i_dim, i_atom, n_loop, i_extra, i_loop,
     & flag_status, i_neigh_i, i_neigh_j, i_neigh_k

      double precision junk, dist, 
     & grid_pos(n_dim), grid_pos_tmp(n_dim), sum_V_pot, 
     & delta_dist(n_dim), charge(n_atoms_max), sum_q

      double precision, allocatable, dimension(:,:,:) :: V_pot

      character*4 p

      sum_q = 0.d0
      write(*,*) "Enter file name"
      read(*,*) file_name

      open(10,file=file_name,status="old")

       read(10,*)
       read(10,*)
       read(10,*) n_atoms, (axis_zero(i_dim), i_dim=1,n_dim)       
       read(10,*) n_grid(1), (axis_vector(1,i_dim), i_dim=1,n_dim)
       read(10,*) n_grid(2), (axis_vector(2,i_dim), i_dim=1,n_dim)
       read(10,*) n_grid(3), (axis_vector(3,i_dim), i_dim=1,n_dim)
       do i_atom=1, n_atoms
        read(10,*) atom_index(i_atom), junk,
     & (atom_pos(i_atom,i_dim), i_dim=1, n_dim), charge(i_atom)
        sum_q = sum_q + charge(i_atom)
       end do
       n_loop = floor(real(n_grid(3))/6)
       i_extra = mod(n_grid(3),6)
       write(p,'(i1)') i_extra
       write(*,*) p, i_extra, n_loop
       close(10)
       write(*,*) "Total charge: ", sum_q

       end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
      subroutine scancif(ncif,ciffile,n_atoms, loopcount)
c*********************************************************************
c
c     routine to scan cif file to determine number of atoms
c
c*********************************************************************
      implicit none
      logical done,safe,loopchk,atmchk,entrychk,chgchk
      integer ncif,n_atoms
      integer i,idum, loop_count, iatm, loopcount
      character*100 ciffile
      character(len=lenrec), dimension(100) :: loops

      done=.false.
      loopchk=.false.
      atmchk=.false.
      chgchk=.false.
      entrychk = .false.
      open(ncif,file=ciffile,status='old')
      iatm=0
 
      do while(.not.done)
        call getrec(safe,ncif)
        call strip(record,lenrec)
        if(loopchk)then
            if(findstring('_',record,idum))then
                loop_count = loop_count + 1
c               check for atomic coordinate positions to flag
c               atom entries
                call copystring(record, loops(loop_count),lenrec)
                if((findstring('atom', record, idum)).and.
     &              (findstring('_x', record, idum)))then
                    atmchk=.true.
                else if ((findstring('_atom_site_charge',record,idum))
     &      .or.(findstring('_atom_type_charge',record,idum)).or.
     &      (findstring('_atom_type_partial_charge',record,idum)))then
                    chgchk=.true.
                end if
            else
                if(atmchk)then
                  allocate(loopatom(loop_count))
                  do i=1,loop_count
                    loopatom(i)=loops(i)
                  end do
                  loopcount=loop_count
                end if
                loopchk=.false.
                entrychk=.true.
            end if

        endif
        if(findstring('loop_',record,idum))then
            loop_count = 0
            entrychk=.false.
            atmchk=.false.
            loopchk=.true.
        end if
        if((.not.loopchk).and.(entrychk))then
            if((record(1).eq.'#').or.(record(1).eq.' '))then
                entrychk=.false.
                atmchk=.false.
                loopchk=.false.
            else
                if(atmchk)then
                    iatm = iatm + 1
c                    write(*,*)record
                end if
            end if
        end if
        if(.not.safe)done=.true.
      end do

      close(ncif)
      if(.not.chgchk)then
         write(*,*) "ERROR: there were no charges in the cif"
         call exit(1)
      end if
      n_atoms=iatm 
      return
      end subroutine scancif

      subroutine readcif(ncif,ciffile,n_atoms,loopcount,
     & alen,blen,clen,alph,beta,gamma,title,xfrac,yfrac,zfrac)
c*********************************************************************
c
c     routine to read cif files and store all the values in memory 
c
c*********************************************************************
      implicit none
      logical done,safe,loopchk,atmchk,entrychk
      logical xfrac, yfrac, zfrac
      integer ncif,n_atoms,iatm,idum,loop_count
      integer i,ijunk,rmkcnt,loopcount
      real(8) xcoord,ycoord,zcoord,charge, rjunk
      real(8) alen,blen,clen,alph,beta,gamma
      character*8 atm,atmn
      character(len=lenrec) cjunk
      character*1 title(80)
      character*100 ciffile

      done=.false.
      loopchk=.false.
      atmchk=.false.
      entrychk = .false.
      xfrac = .false.
      yfrac = .false.
      zfrac = .false.
      open(ncif,file=ciffile,status='old')
      iatm=0
 
      do while(.not.done)
        call getrec(safe,ncif)
c        call lowcase(record,7)
        call strip(record,lenrec)
        if(findstring('data_',record,idum))then
          rmkcnt=rmkcnt+1
          if(rmkcnt.eq.1)then
            call getword(cjunk,record,6,lenrec)
            call copystring(record,title(1),80)
          endif
        elseif(findstring('_cell_length_a',record,idum))then
          call getword(cjunk,record,14,lenrec)  
          alen=dblstr(record,lenrec,idum)
        elseif(findstring('_cell_length_b',record,idum))then
          call getword(cjunk,record,14,lenrec)  
          blen=dblstr(record,lenrec,idum)
        elseif(findstring('_cell_length_c',record,idum))then
          call getword(cjunk,record,14,lenrec)  
          clen=dblstr(record,lenrec,idum)
        elseif(findstring('_cell_angle_alpha',record,idum))then
          call getword(cjunk,record,17,lenrec)  
          alph=dblstr(record,lenrec,idum)
        elseif(findstring('_cell_angle_beta',record,idum))then
          call getword(cjunk,record,16,lenrec)  
          beta=dblstr(record,lenrec,idum)
        elseif(findstring('_cell_angle_gamma',record,idum))then
          call getword(cjunk,record,17,lenrec)  
          gamma=dblstr(record,lenrec,idum)
        end if
        if(loopchk)then
            if(findstring('_',record,idum))then
                loop_count = loop_count + 1
c               check for atomic coordinate positions to flag
c               atom entries
                if((findstring('atom', record, idum)).and.
     &              (findstring('_x', record, idum)))then
                    atmchk=.true.
                end if
            else
                loopchk=.false.
                entrychk=.true.
            end if

        endif
        if(findstring('loop_',record,idum))then
            loop_count = 0
            entrychk=.false.
            atmchk=.false.
            loopchk=.true.
        end if
        if((.not.loopchk).and.(entrychk))then
            if((record(1).eq.'#').or.(record(1).eq.' '))then
                entrychk=.false.
                atmchk=.false.
                loopchk=.false.
            else
                if(atmchk)then
                    iatm = iatm + 1
                    atom_index(iatm) = iatm
                    do i=1, loopcount
                      call getword(cjunk,record,lenrec,lenrec)
                      if ((findstring('atom', loopatom(i), idum)).and.
     &                   (findstring('_x', loopatom(i), idum)))then
                         if (findstring('frac', loopatom(i),idum))then
                            xfrac=.true.
                            atom_pos_frac(iatm,1) = 
     &dblstr(cjunk,lenrec,idum)
                         else
                            atom_pos(iatm,1) = dblstr(cjunk,lenrec,idum)
                         end if
                      else if ((findstring('atom', loopatom(i), idum))
     &                   .and.
     &                   (findstring('_y', loopatom(i), idum)))then
                         if (findstring('frac', loopatom(i),idum))then
                            yfrac=.true. 
                            atom_pos_frac(iatm,2) = 
     &dblstr(cjunk,lenrec,idum)
                         else
                            atom_pos(iatm,2) = dblstr(cjunk,lenrec,idum)
                         end if
                      else if ((findstring('atom', loopatom(i), idum))
     &                   .and.
     &                   (findstring('_z', loopatom(i), idum)))then
                         if (findstring('frac', loopatom(i),idum))then
                            zfrac=.true. 
                            atom_pos_frac(iatm,3) = 
     &dblstr(cjunk,lenrec,idum)
                         else
                            atom_pos(iatm,3) = dblstr(cjunk,lenrec,idum)
                         end if
                      else if ((findstring('atom', loopatom(i), idum))
     &                   .and.
     &                   (findstring('symbol', loopatom(i), idum)))then
                         call getword(atm,cjunk,2,lenrec)
                         atom(iatm) = atm
                         mass(iatm)=atmmass(atm)
                         atom_number(iatm) = atmnumber(mass(iatm))
                      else if ((findstring('atom', loopatom(i), idum))
     &                   .and.
     &                   (findstring('label', loopatom(i), idum)))then
                         call getword(atmn,cjunk,4,lenrec)
                         atmname(iatm) = atmn
                      else if ((findstring('atom', loopatom(i), idum))
     &                   .and.
     &                   (findstring('charge', loopatom(i), idum)))then
                         q_part(iatm) = dblstr(cjunk,lenrec,idum)
                      end if
                    end do
                end if
            end if
        end if
        if(.not.safe)done=.true.
      end do
      close(ncif)
      
      return
      end subroutine readcif

      logical function findstring(seek,string,here)

c***********************************************************************
c     
c     routine to find an explicit string in an input record
c     note: variable `seek' is a character string while variable
c    `string' is a character*1 array i.e. code is application specific
c     
c***********************************************************************

      implicit none

      integer i,n,m,here
      character*(*) seek
      character*1 string(lenrec)

      m=lenrec
      n=len(seek)
      findstring=.false.

      here=0
      do while(here.le.m-n.and.(.not.findstring))

        findstring=.true.

        do i=1,n
          if(seek(i:i).ne.string(here+i))findstring=.false.
        enddo

        here=here+1

      enddo

      return
      end function findstring
      
      subroutine copystring(oldstr,newstr,length)

c***********************************************************************
c     
c     routine to copy one string into another
c     
c***********************************************************************

      implicit none

      character*1 newstr(*),oldstr(*)
      integer i,length

      do i=1,length

        newstr(i)=oldstr(i)

      enddo

      return
      end subroutine copystring
      
      subroutine getword(word,string,len1,len2)

c***********************************************************************
c     
c     routine to fetch an 8 character word from a string
c     while ignoring leading blanks
c     
c***********************************************************************

      implicit none

      logical final
      character*8 word
      integer len1,len2,i,j,k
      character*1 wrdseq(len1),string(len2)
c      character*1 word(len1),string(len2)
      do i=1,len1
        wrdseq(i)=' '
c         word(i)=' '
      enddo

      i=0
      k=0
      final=.false.
      
      do while(.not.final.and.i.lt.len2)
        
        i=i+1
        
        if(string(1).eq.' ')then
          
          if(k.gt.0)final=.true.
          
        else
          
          k=k+1
          wrdseq(k)=string(1)
c          word(k)=string(1)
          if(k.eq.len1)final=.true.

        endif
        
        do j=1,len2-1
          
          string(j)=string(j+1)
          
        enddo
        
        string(len2)=' '
          
      enddo
      
      word=mkwd8(wrdseq)

      return
      end subroutine getword
      
      subroutine getrec(safe,ifile)

c*********************************************************************
c     
c      subroutine to read a character string on one node
c     
c*********************************************************************

      implicit none
      
      logical safe

      character*150 line
      integer ifile,i
      
      safe=.true.
      
      read(ifile,'(a150)',end=100)line
 
      do i=1,lenrec

        record(i)=line(i:i)
         
      enddo
             
      return
        
  100 safe=.false.
                  
      end subroutine getrec


      subroutine lowcase(string,length)

c***********************************************************************
c     
c     routine to lowercase a string of up to 255 characters.
c     Transportable to non-ASCII machines
c     
c***********************************************************************

      implicit none

      character*1 string(*)
      character*1 letter
      integer i,length

      do i=1,min(255,length)

        letter=string(i)

        if(letter.eq.'A')then
          letter='a'
        else if(letter.eq.'B')then
          letter='b'
        else if(letter.eq.'C')then
          letter='c'
        else if(letter.eq.'D')then
          letter='d'
        else if(letter.eq.'E')then
          letter='e'
        else if(letter.eq.'F')then
          letter='f'
        else if(letter.eq.'G')then
          letter='g'
        else if(letter.eq.'H')then
          letter='h'
        else if(letter.eq.'I')then
          letter='i'
        else if(letter.eq.'J')then
          letter='j'
        else if(letter.eq.'K')then
          letter='k'
        else if(letter.eq.'L')then
          letter='l'
        else if(letter.eq.'M')then
          letter='m'
        else if(letter.eq.'N')then
          letter='n'
        else if(letter.eq.'O')then
          letter='o'
        else if(letter.eq.'P')then
          letter='p'
        else if(letter.eq.'Q')then
          letter='q'
        else if(letter.eq.'R')then
          letter='r'
        else if(letter.eq.'S')then
          letter='s'
        else if(letter.eq.'T')then
          letter='t'
        else if(letter.eq.'U')then
          letter='u'
        else if(letter.eq.'V')then
          letter='v'
        else if(letter.eq.'W')then
          letter='w'
        else if(letter.eq.'X')then
          letter='x'
        else if(letter.eq.'Y')then
          letter='y'
        else if(letter.eq.'Z')then
          letter='z'
        endif

        string(i)=letter

      enddo

      return
      end subroutine lowcase
       
      subroutine getcell(alen,blen,clen,alph,beta,gam) 
c***********************************************************************
c     
c     routine to convert to cell vectors 
c     from alen,blen,clen
c     alph = angle between b and c
c     beta  = angle between a and c
c     gamma = angle between a and b
c
c***********************************************************************
      implicit none
      integer i
      real(8), parameter :: deg2rad=0.0174532925d0
      real(8), parameter :: rad2deg=57.2957795d0
      real(8), parameter :: pi=3.141592653589793d0
      real(8) alen,blen,clen,alph,beta,gam
      real(8) pr
      real(8), dimension(3) :: avec,bvec,cvec,v1,v2,v3
      double precision a(n_dim), b(n_dim), c(n_dim), ab(n_dim),
     & r(n_dim), cross_tmp(n_dim), signo, dot, dot_1,
     & cross_tmp_1(n_dim)

      alph=alph*deg2rad
      beta=beta*deg2rad
      gam=gam*deg2rad
c     the a axis is oriented to the xaxis
c     start with orthogonal basis
      v1=(/1.d0,0.d0,0.d0/)

      v2=(/0.d0,1.d0,0.d0/)

      v3=(/0.d0,0.d0,1.d0/)
      
      n_grid(1) = nint(alen/grid_spacing)
      n_grid(2) = nint(blen/grid_spacing)
      n_grid(3) = nint(clen/grid_spacing)

      avec=v1

      bvec=dsin(gam)*v2+dcos(gam)*v1

      pr=(dcos(alph)-dcos(gam)*dcos(beta))/dsin(gam)

      cvec=v3*sqrt(1-dcos(beta)**2-pr**2)
     &+dcos(beta)*v1+pr*v2
c     write to cell
      real_box_vector(1,1:3)=avec*alen*angs2bohr
      real_box_vector(2,1:3)=bvec*blen*angs2bohr
      real_box_vector(3,1:3)=cvec*clen*angs2bohr
c     fix for rounding
      do i=1,3
        real_box_vector(i,1:3)=float(
     &nint(real_box_vector(i,1:3)*1.d6))/1.d6
      enddo
c     create axis vectors
      axis_vector(1,1:3) = real_box_vector(1,1:3)/n_grid(1)
      axis_vector(2,1:3) = real_box_vector(2,1:3)/n_grid(2)
      axis_vector(3,1:3) = real_box_vector(3,1:3)/n_grid(3)

!-----Creating the reciprocal space box
      a(1) = real_box_vector(1,1)
      a(2) = real_box_vector(1,2)
      a(3) = real_box_vector(1,3)
      b(1) = real_box_vector(2,1)
      b(2) = real_box_vector(2,2)
      b(3) = real_box_vector(2,3)
      c(1) = real_box_vector(3,1)
      c(2) = real_box_vector(3,2)
      c(3) = real_box_vector(3,3)
      call Vect_Cross(b,c,ab)
      Box_Volume = Vect_Dot(a,ab)
!-----Creating the reciprocal space box
      call Vect_Cross(a,b,ab)      
      do i_dim=1, n_dim
       recip_box_vector(3,i_dim) = (2*pi/Box_volume)*ab(i_dim)
      end do
      call Vect_Cross(c,a,ab)
      do i_dim=1, n_dim
       recip_box_vector(2,i_dim) = (2*pi/Box_volume)*ab(i_dim)
      end do      
      call Vect_Cross(b,c,ab)      
      do i_dim=1, n_dim
       recip_box_vector(1,i_dim) = (2*pi/Box_volume)*ab(i_dim)
      end do
c     INTEL FORTRAN SOMTIMES CHANGES THE INDICES OF THE VECTORS
C     CREATING INCONSISTENCIES IN THE INVERSION.
      !write(*,*)real_box_vector 
      !write(*,*)recip_box_vector
      return
      end subroutine getcell

      subroutine strip(string,imax)

c***********************************************************************
c     
c     Routine to strip blanks from start of a string
c     maximum length is 255 characters
c     
c***********************************************************************

      implicit none

      integer i,imax,j
      character*1 string(imax)
      do i=1,imax
    
        if(string(1).eq.' ')then

          do j=1,imax-1

            string(j)=string(j+1)

          enddo

          string(imax)=' '

        endif

      enddo

      return
      end subroutine strip
      
      character*8 function mkwd8(string)

c***********************************************************************
c     
c     Routine to make an 8 character word from a string
c
c***********************************************************************

      implicit none

      integer i
      character*1 string(*)
      
      do i=1,8
         mkwd8(i:i)=string(i)
      enddo
      
      return
      end function mkwd8
      
      integer function atmnumber(i) 
c*******************************************************************
c     generate an atomic number from an atomic mass 
c     EDIT (pb 09/01/13): this function reads the mass reported
c     on a standard periodic table and assigns an atomic number.
c     You will run into problems if you are using atomic masses
c     of isotopes in the FIELD file.
c*******************************************************************
      implicit none
      real(8) i
      if ((i.ge.0.0).and.(i.le.1.5))then
        atmnumber=1
      elseif((i.ge.3.9).and.(i.le.4.5))then
        atmnumber=2
      elseif((i.ge.6.5).and.(i.le.7.1))then
        atmnumber=3
      elseif((i.ge.8.9).and.(i.le.9.5))then
        atmnumber=4
      elseif((i.ge.10.5).and.(i.le.11.1))then
        atmnumber=5
      elseif((i.ge.11.9).and.(i.le.12.5))then
        atmnumber=6
      elseif((i.ge.13.9).and.(i.le.14.5))then
        atmnumber=7
      elseif((i.ge.15.5).and.(i.le.16.1))then
        atmnumber=8
      elseif((i.ge.18.5).and.(i.le.19.1))then
        atmnumber=9
      elseif((i.ge.19.9).and.(i.le.20.5))then
        atmnumber=10
      elseif((i.ge.22.5).and.(i.le.23.1))then
        atmnumber=11
      elseif((i.ge.23.9).and.(i.le.24.5))then
        atmnumber=12
      elseif((i.ge.26.5).and.(i.le.27.1))then
        atmnumber=13
      elseif((i.ge.27.9).and.(i.le.28.5))then
        atmnumber=14
      elseif((i.ge.30.5).and.(i.le.31.1))then
        atmnumber=15
      elseif((i.ge.31.9).and.(i.le.32.5))then
        atmnumber=16
      elseif((i.ge.34.9).and.(i.le.36.1))then
        atmnumber=17
c     Ar (18) has mass range that overlaps with Ca (20). Be careful 
c     with mass rounding here!
      elseif((i.ge.39.5).and.(i.le.39.9999))then
        atmnumber=18
      elseif((i.ge.38.9).and.(i.le.39.4))then
        atmnumber=19
      elseif((i.ge.40.0).and.(i.le.40.5))then
        atmnumber=20
      elseif((i.ge.44.5).and.(i.le.45.1))then
        atmnumber=21
      elseif((i.ge.47.5).and.(i.le.48.1))then
        atmnumber=22
      elseif((i.ge.50.5).and.(i.le.51.1))then
        atmnumber=23
      elseif((i.ge.51.5).and.(i.le.52.1))then
        atmnumber=24
      elseif((i.ge.54.5).and.(i.le.55.1))then
        atmnumber=25
      elseif((i.ge.55.5).and.(i.le.56.1))then
        atmnumber=26
c     Co (27) and Ni (28) have very close mass ranges
      elseif((i.ge.58.76).and.(i.le.59.1))then
        atmnumber=27
      elseif((i.ge.58.5).and.(i.le.59.75))then
        atmnumber=28
      elseif((i.ge.62.9).and.(i.le.64.1))then
        atmnumber=29
      elseif((i.ge.64.9).and.(i.le.66.1))then
        atmnumber=30
      elseif((i.ge.69.5).and.(i.le.70.1))then
        atmnumber=31
      elseif((i.ge.72.5).and.(i.le.73.1))then
        atmnumber=32
      elseif((i.ge.74.5).and.(i.le.75.1))then
        atmnumber=33
      elseif((i.ge.78.5).and.(i.le.79.1))then
        atmnumber=34
      elseif((i.ge.79.5).and.(i.le.80.1))then
        atmnumber=35
      elseif((i.ge.83.5).and.(i.le.84.1))then
        atmnumber=36
      elseif((i.ge.84.9).and.(i.le.86.1))then
        atmnumber=37
      elseif((i.ge.87.5).and.(i.le.88.1))then
        atmnumber=38
      elseif((i.ge.88.5).and.(i.le.89.1))then
        atmnumber=39
      elseif((i.ge.90.9).and.(i.le.91.5))then
        atmnumber=40
      elseif((i.ge.92.5).and.(i.le.93.1))then
        atmnumber=41
      elseif((i.ge.95.5).and.(i.le.96.1))then
        atmnumber=42
      elseif((i.ge.97.9).and.(i.le.98.1))then
        atmnumber=43
      elseif((i.ge.109.9).and.(i.le.101.5))then
        atmnumber=44
      elseif((i.ge.102.5).and.(i.le.103.1))then
        atmnumber=45
      elseif((i.ge.105.9).and.(i.le.106.5))then
        atmnumber=46
      elseif((i.ge.107.5).and.(i.le.108.1))then
        atmnumber=47
      elseif((i.ge.111.9).and.(i.le.112.5))then
        atmnumber=48
      elseif((i.ge.114.5).and.(i.le.115.1))then
        atmnumber=49
      elseif((i.ge.118.5).and.(i.le.119.1))then
        atmnumber=50
      elseif((i.ge.121.5).and.(i.le.122.1))then
        atmnumber=51
      elseif((i.ge.127.5).and.(i.le.128.1))then
        atmnumber=52
      elseif((i.ge.126.5).and.(i.le.127.1))then
        atmnumber=53
      elseif((i.ge.130.9).and.(i.le.131.5))then
        atmnumber=54
      elseif((i.ge.132.5).and.(i.le.133.1))then
        atmnumber=55
      elseif((i.ge.136.9).and.(i.le.137.5))then
        atmnumber=56
      elseif((i.ge.138.5).and.(i.le.139.1))then
        atmnumber=57
      elseif((i.ge.139.9).and.(i.le.140.5))then
        atmnumber=58
      elseif((i.ge.140.6).and.(i.le.141.1))then
        atmnumber=59
      elseif((i.ge.144.0).and.(i.le.144.5))then
        atmnumber=60
      elseif((i.ge.144.9).and.(i.le.145.1))then
        atmnumber=61
      elseif((i.ge.150.0).and.(i.le.150.6))then
        atmnumber=62
      elseif((i.ge.151.5).and.(i.le.152.1))then
        atmnumber=63
      elseif((i.ge.156.9).and.(i.le.157.5))then
        atmnumber=64
      elseif((i.ge.158.5).and.(i.le.159.1))then
        atmnumber=65
      elseif((i.ge.162.0).and.(i.le.163.1))then
        atmnumber=66
      elseif((i.ge.164.5).and.(i.le.165.1))then
        atmnumber=67
      elseif((i.ge.166.5).and.(i.le.167.9))then
        atmnumber=68
      elseif((i.ge.168.0).and.(i.le.169.1))then
        atmnumber=69
      elseif((i.ge.172.9).and.(i.le.173.5))then
        atmnumber=70
      elseif((i.ge.174.0).and.(i.le.175.1))then
        atmnumber=71
      elseif((i.ge.178.0).and.(i.le.179.1))then
        atmnumber=72
      elseif((i.ge.180.0).and.(i.le.181.1))then
        atmnumber=73
      elseif((i.ge.183.0).and.(i.le.184.1))then
        atmnumber=74
      elseif((i.ge.185.9).and.(i.le.186.5))then
        atmnumber=75
      elseif((i.ge.189.9).and.(i.le.190.5))then
        atmnumber=76
      elseif((i.ge.191.9).and.(i.le.192.5))then
        atmnumber=77
      elseif((i.ge.194.9).and.(i.le.195.5))then
        atmnumber=78
      elseif((i.ge.196.5).and.(i.le.197.1))then
        atmnumber=79
      elseif((i.ge.200.0).and.(i.le.201.1))then
        atmnumber=80
      elseif((i.ge.203.9).and.(i.le.204.6))then
        atmnumber=81
      elseif((i.ge.206.9).and.(i.le.207.6))then
        atmnumber=82
      elseif((i.ge.208.5).and.(i.le.209.1))then
        atmnumber=83
c     Po atomic number 84 has the same mass range as Bi (83)
      elseif((i.ge.209.9).and.(i.le.210.1))then
        atmnumber=85
      elseif((i.ge.221.9).and.(i.le.222.1))then
        atmnumber=86
      elseif((i.ge.222.9).and.(i.le.223.1))then
        atmnumber=87
      elseif((i.ge.225.9).and.(i.le.226.1))then
        atmnumber=88
      elseif((i.ge.226.9).and.(i.le.227.1))then
        atmnumber=89
      elseif((i.ge.231.9).and.(i.le.232.1))then
        atmnumber=90
      elseif((i.ge.230.9).and.(i.le.231.1))then
        atmnumber=91
      elseif((i.ge.237.9).and.(i.le.238.1))then
        atmnumber=92
c     Np atomic number 93 has the same mass range as U (92)
      elseif((i.ge.243.9).and.(i.le.244.1))then
        atmnumber=94
      elseif((i.ge.242.9).and.(i.le.243.1))then
        atmnumber=95
      elseif((i.ge.246.9).and.(i.le.247.1))then
        atmnumber=96
c     Bk atomic number 97 has the same mass range as Cm (96)
      elseif((i.ge.250.9).and.(i.le.251.1))then
        atmnumber=98
      elseif((i.ge.251.9).and.(i.le.252.1))then
        atmnumber=99
      elseif((i.ge.256.9).and.(i.le.257.1))then
        atmnumber=100
      elseif((i.ge.257.9).and.(i.le.258.1))then
        atmnumber=101
      elseif((i.ge.258.9).and.(i.le.259.1))then
        atmnumber=102
      elseif((i.ge.261.9).and.(i.le.262.1))then
        atmnumber=103
      elseif((i.ge.266.9).and.(i.le.267.1))then
        atmnumber=104
      elseif((i.ge.267.9).and.(i.le.268.1))then
        atmnumber=105
      elseif((i.ge.270.9).and.(i.le.271.1))then
        atmnumber=106
      elseif((i.ge.269.9).and.(i.le.270.1))then
        atmnumber=107
      elseif((i.ge.268.9).and.(i.le.269.1))then
        atmnumber=108
      elseif((i.ge.277.9).and.(i.le.278.1))then
        atmnumber=109
      elseif((i.ge.280.9).and.(i.le.281.1))then
        atmnumber=110
c     Rg atomic number 112 has the same mass range as Ds (110)
      elseif((i.ge.284.9).and.(i.le.285.1))then
        atmnumber=112
      elseif((i.ge.285.9).and.(i.le.286.1))then
        atmnumber=113
      elseif((i.ge.288.9).and.(i.le.289.1))then
        atmnumber=114
c     Uup atomic number 115 has the same mass range as Fl (114)
      elseif((i.ge.292.9).and.(i.le.293.1))then
        atmnumber=116
      elseif((i.ge.293.9).and.(i.le.294.1))then
        atmnumber=117
c     Uuo atomic number 118 has the same mass range as Uus (117)
      endif
      return
      end function atmnumber

      subroutine timchk(ktim,time)

c***********************************************************************
c     
c     Routine for time elapsed in seconds
c     
c***********************************************************************
      implicit none

      character*12 dat,tim,zon
      integer mynode,ktim,day
      real(8) time,told,tnow
      integer info(8)

      save day

      call date_and_time(dat,tim,zon,info)
       
      if(ktim.eq.0)then

         day=info(3)
         time=3600.d0*dble(info(5))+60.d0*dble(info(6))+
     x         dble(info(7))+0.001d0*dble(info(8))

      else 
         told=time
         tnow=3600.d0*dble(info(5))+60.d0*dble(info(6))+
     x         dble(info(7))+0.001d0*dble(info(8))
         if(day.ne.info(3))then
           told=told-86400.d0
           day=info(3)
         endif
         time=tnow-told
      endif

      return
      end subroutine timchk

      subroutine getroot(ciffile,ind)
c*********************************************************************
c     gather the cif filename before the .pdb
c     this will enable the writing of other output files
c*********************************************************************
      implicit none
      integer i,ind
      character*100 ciffile

      do i=1,100
        if(ciffile(i:i+3).eq.'.cif')then
           ind=i-1
        endif
      enddo
      allocate(rootname(ind))
      do i=1,ind
        rootname(i)=ciffile(i:i)
      enddo
      return
      end subroutine getroot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----Subroutine that assigns the van der Waals radii for the
!-----elements found in the cube file according to the UFF tabulation 
      subroutine VDW_radii_array(atoms_array, vdw_radii)
      integer atoms_array(n_atoms_max), i
      double precision vdw_radii(n_atoms_max)

      double precision vdw_radii_file(n_atoms_type_max) 

!-----Hardcore tabulation taken from the UFF with
!-----the VDW radii for the elements in the periodic table
!-----ordered according to their atomic number
!-----UNITS ARE IN BOHR
      vdw_radii_file(1) = 2.72687
      vdw_radii_file(2) = 2.23177
      vdw_radii_file(3) = 2.31586
      vdw_radii_file(4) = 2.59365
      vdw_radii_file(5) = 3.85788
      vdw_radii_file(6) = 3.63867
      vdw_radii_file(7) = 3.4582
      vdw_radii_file(8) = 3.30702
      vdw_radii_file(9) = 3.17852
      vdw_radii_file(10) = 3.06419
      vdw_radii_file(11) = 2.81853
      vdw_radii_file(12) = 2.85443
      vdw_radii_file(13) = 4.25094
      vdw_radii_file(14) = 4.05819
      vdw_radii_file(15) = 3.91835
      vdw_radii_file(16) = 3.81252
      vdw_radii_file(17) = 3.72937
      vdw_radii_file(18) = 3.65473
      vdw_radii_file(19) = 3.60182
      vdw_radii_file(20) = 3.21159
      vdw_radii_file(21) = 3.11332
      vdw_radii_file(22) = 2.99994
      vdw_radii_file(23) = 2.97065
      vdw_radii_file(24) = 2.85632
      vdw_radii_file(25) = 2.79774
      vdw_radii_file(26) = 2.75144
      vdw_radii_file(27) = 2.71365
      vdw_radii_file(28) = 2.67774
      vdw_radii_file(29) = 3.3023
      vdw_radii_file(30) = 2.61066
      vdw_radii_file(31) = 4.14133
      vdw_radii_file(32) = 4.04401
      vdw_radii_file(33) = 3.99677
      vdw_radii_file(34) = 3.97315
      vdw_radii_file(35) = 3.95803
      vdw_radii_file(36) = 3.91268
      vdw_radii_file(37) = 3.88717
      vdw_radii_file(38) = 3.44025
      vdw_radii_file(39) = 3.16057
      vdw_radii_file(40) = 2.95175
      vdw_radii_file(41) = 2.99049
      vdw_radii_file(42) = 2.88372
      vdw_radii_file(43) = 2.8327
      vdw_radii_file(44) = 2.79963
      vdw_radii_file(45) = 2.7675
      vdw_radii_file(46) = 2.73916
      vdw_radii_file(47) = 2.97443
      vdw_radii_file(48) = 2.69097
      vdw_radii_file(49) = 4.21692
      vdw_radii_file(50) = 4.14984
      vdw_radii_file(51) = 4.17629
      vdw_radii_file(52) = 4.22354
      vdw_radii_file(53) = 4.25188
      vdw_radii_file(54) = 4.16118
      vdw_radii_file(55) = 4.26795
      vdw_radii_file(56) = 3.49883
      vdw_radii_file(57) = 3.32781
      vdw_radii_file(58) = 3.35993
      vdw_radii_file(59) = 3.40718
      vdw_radii_file(60) = 3.37789
      vdw_radii_file(61) = 3.35143
      vdw_radii_file(62) = 3.32592
      vdw_radii_file(63) = 3.30041
      vdw_radii_file(64) = 3.1823
      vdw_radii_file(65) = 3.26072
      vdw_radii_file(66) = 3.23899
      vdw_radii_file(67) = 3.22104
      vdw_radii_file(68) = 3.20403
      vdw_radii_file(69) = 3.18797
      vdw_radii_file(70) = 3.17002
      vdw_radii_file(71) = 3.4393
      vdw_radii_file(72) = 2.96781
      vdw_radii_file(73) = 2.99522
      vdw_radii_file(74) = 2.89978
      vdw_radii_file(75) = 2.79113
      vdw_radii_file(76) = 2.94797
      vdw_radii_file(77) = 2.68341
      vdw_radii_file(78) = 2.60215
      vdw_radii_file(79) = 3.11143
      vdw_radii_file(80) = 2.55585
      vdw_radii_file(81) = 4.10732
      vdw_radii_file(82) = 4.06008
      vdw_radii_file(83) = 4.12905
      vdw_radii_file(84) = 4.44936
      vdw_radii_file(85) = 4.4881
      vdw_radii_file(86) = 4.50227
      vdw_radii_file(87) = 4.62983
      vdw_radii_file(88) = 3.47426
      vdw_radii_file(89) = 3.28623
      vdw_radii_file(90) = 3.20875
      vdw_radii_file(91) = 3.23521
      vdw_radii_file(92) = 3.20781
      vdw_radii_file(93) = 3.23521
      vdw_radii_file(94) = 3.23521
      vdw_radii_file(95) = 3.19458
      vdw_radii_file(96) = 3.14261
      vdw_radii_file(97) = 3.1549
      vdw_radii_file(98) = 3.13033
      vdw_radii_file(99) = 3.1171
      vdw_radii_file(100) = 3.10482
      vdw_radii_file(101) = 3.09348
      vdw_radii_file(102) = 3.06892
      vdw_radii_file(103) = 3.05758

!-----Assigning the corresponding VDW radii to the elements found in
!-----the ESP cube file
      do i=1, n_atoms
       vdw_radii(i) =  vdw_radii_file(atoms_array(i))
      end do

      end subroutine
c######################################################################
c   END OF PROGRAM
c######################################################################

      end program 

