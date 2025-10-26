module kind_precision_module
    use, intrinsic :: iso_fortran_env, only: i8=>int8, i16=>int16, i32=>int32, i64=>int64
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64, qp=>real128

    implicit none
  
    private
    public :: sp, dp, qp, i8, i16, i32, i64
  
end module kind_precision_module