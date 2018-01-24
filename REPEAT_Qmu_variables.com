      
!-----General parameters
      integer n_dim, n_atoms_max, n_atoms_type_max, voxels_max_1d,
     & MAXK, KMAX, KSQMAX, order_type_max, type_list_max
       parameter(n_dim=3,n_atoms_max=900,n_atoms_type_max=107,
     & voxels_max_1d=3000000,MAXK=5000,KMAX=7,KSQMAX=49,
     & order_type_max=50,type_list_max=80)

      integer n_grid(n_dim), n_atoms, atom_index(n_atoms_max),
     & NMAX(n_dim)
      common/int_general/n_grid, n_atoms, atom_index, NMAX

      double precision axis_zero(n_dim), axis_vector(n_dim,n_dim),
     & atom_pos(n_atoms_max,n_dim), real_box_vector(n_dim,n_dim), 
     & Box_volume, recip_box_vector(n_dim,n_dim), KVEC(MAXK)
      common/dbl_general/axis_zero, axis_vector, atom_pos, 
     & real_box_vector, Box_volume, recip_box_vector, KVEC  


!-----Mathematical constants
      double precision pi, TWOPI
      parameter(pi=3.14159265,TWOPI=6.2831853)


!-----Input section
      integer flag_cutoff, fit_RESP, symm_flag, const_flag, 
     & q_tot(4)
      common/int_input/flag_cutoff, fit_RESP, symm_flag, const_flag, 
     & q_tot

      double precision vdw_fact, R_cutoff
      common/dbl_input/vdw_fact, R_cutoff

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

!-----Variables related to VDW params
      double precision vdw_radii(n_atoms_max)
      common/vdws_dbl/vdw_radii

