module nnmod_pes1
     integer(kind=4),parameter :: intype=4,retype=8,atomdim=3
     real(kind=retype) :: pi=dacos(-1.d0)
     character(kind=1,len=80) :: pesdir='./pes1_ag44'
     integer(kind=intype) :: pesid=1
     integer(kind=intype),parameter :: max_iep=6,max_ifail=6
     include "pes_hd.inc1.f90"
end module

module nnmod_pes2
     integer(kind=4),parameter :: intype=4,retype=8,atomdim=3
     real(kind=retype) :: pi=dacos(-1.d0)
     character(kind=1,len=80) :: pesdir='./pes2_4ftp'
     integer(kind=intype) :: pesid=2
     integer(kind=intype),parameter :: max_iep=6,max_ifail=6
     include "pes_hd.inc2.f90"
end module

module nnmod_pes3
     integer(kind=4),parameter :: intype=4,retype=8,atomdim=3
     real(kind=retype) :: pi=dacos(-1.d0)
     character(kind=1,len=80) :: pesdir='./pes3_2ftp'
     integer(kind=intype) :: pesid=1
     integer(kind=intype),parameter :: max_iep=6,max_ifail=6
     include "pes_hd.inc3.f90"
end module

module nnmod_pes12
     integer(kind=4),parameter :: intype=4,retype=8,atomdim=3
     real(kind=retype) :: pi=dacos(-1.d0)
     character(kind=1,len=80) :: pesdir='./pes12_ag44_4ftp'
     integer(kind=intype) :: pesid=1
     integer(kind=intype),parameter :: max_iep=6,max_ifail=6
     include "pes_hd.inc12.f90"
end module

module nnmod_pes13
     integer(kind=4),parameter :: intype=4,retype=8,atomdim=3
     real(kind=retype) :: pi=dacos(-1.d0)
     character(kind=1,len=80) :: pesdir='./pes13_ag44_2ftp'
     integer(kind=intype) :: pesid=1
     integer(kind=intype),parameter :: max_iep=6,max_ifail=6
     include "pes_hd.inc13.f90"
end module

module nnmod_pes22
     integer(kind=4),parameter :: intype=4,retype=8,atomdim=3
     real(kind=retype) :: pi=dacos(-1.d0)
     character(kind=1,len=80) :: pesdir='./pes22_4ftp4ftp'
     integer(kind=intype) :: pesid=2
     integer(kind=intype),parameter :: max_iep=6,max_ifail=6
     include "pes_hd.inc22.f90"
end module

module nnmod_pes23
     integer(kind=4),parameter :: intype=4,retype=8,atomdim=3
     real(kind=retype) :: pi=dacos(-1.d0)
     character(kind=1,len=80) :: pesdir='./pes23_2ftp4ftp'
     integer(kind=intype) :: pesid=2
     integer(kind=intype),parameter :: max_iep=6,max_ifail=6
     include "pes_hd.inc23.f90"
end module

module nnmod_pes33
     integer(kind=4),parameter :: intype=4,retype=8,atomdim=3
     real(kind=retype) :: pi=dacos(-1.d0)
     character(kind=1,len=80) :: pesdir='./pes33_2ftp2ftp'
     integer(kind=intype) :: pesid=1
     integer(kind=intype),parameter :: max_iep=6,max_ifail=6
     include "pes_hd.inc33.f90"
end module

