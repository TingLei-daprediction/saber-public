!===================================================================
module jp_pkind2 ! Integer kinds for helf- and fourth-precision integers
!===================================================================
use mpi
integer,parameter:: hpi=selected_int_kind(3),&
                    fpi=selected_int_kind(2)
end module jp_pkind2
