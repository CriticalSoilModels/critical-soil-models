module mod_test_cmname_parser_suite
   use mod_cmname_parser, only: material_name, integrator_name
   use testdrive, only: new_unittest, unittest_type, error_type, check
   implicit none
   private
   public :: collect_cmname_parser_suite

contains

   subroutine collect_cmname_parser_suite(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("material_name: no integ tag",           test_material_no_tag),        &
         new_unittest("material_name: with integ tag",         test_material_with_tag),      &
         new_unittest("material_name: ortiz_simo suffix",      test_material_ortiz_simo),    &
         new_unittest("integrator_name: no tag gives default", test_integrator_default),     &
         new_unittest("integrator_name: euler suffix",         test_integrator_euler),       &
         new_unittest("integrator_name: ortiz_simo suffix",    test_integrator_ortiz_simo),  &
         new_unittest("material_name: padded cmname",          test_material_padded),        &
         new_unittest("integrator_name: padded cmname",        test_integrator_padded)       &
      ]
   end subroutine collect_cmname_parser_suite

   ! ---------------------------------------------------------------------------

   subroutine test_material_no_tag(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=80) :: cmname
      cmname = "LINEAR_ELASTIC"
      call check(error, trim(material_name(cmname)), "LINEAR_ELASTIC", &
                 more="material_name: no _INTEG_ tag")
   end subroutine test_material_no_tag

   subroutine test_material_with_tag(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=80) :: cmname
      cmname = "MCSS_INTEG_EULER"
      call check(error, trim(material_name(cmname)), "MCSS", &
                 more="material_name: should strip _INTEG_EULER")
   end subroutine test_material_with_tag

   subroutine test_material_ortiz_simo(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=80) :: cmname
      cmname = "MCSS_INTEG_ORTIZ_SIMO"
      call check(error, trim(material_name(cmname)), "MCSS", &
                 more="material_name: should strip _INTEG_ORTIZ_SIMO")
   end subroutine test_material_ortiz_simo

   subroutine test_integrator_default(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=80) :: cmname
      cmname = "LINEAR_ELASTIC"
      call check(error, trim(integrator_name(cmname)), "euler", &
                 more="integrator_name: should default to euler when no _INTEG_ tag")
   end subroutine test_integrator_default

   subroutine test_integrator_euler(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=80) :: cmname
      cmname = "MCSS_INTEG_EULER"
      call check(error, trim(integrator_name(cmname)), "euler", &
                 more="integrator_name: EULER suffix should give lowercase euler")
   end subroutine test_integrator_euler

   subroutine test_integrator_ortiz_simo(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=80) :: cmname
      cmname = "MCSS_INTEG_ORTIZ_SIMO"
      call check(error, trim(integrator_name(cmname)), "ortiz_simo", &
                 more="integrator_name: ORTIZ_SIMO suffix should give lowercase ortiz_simo")
   end subroutine test_integrator_ortiz_simo

   subroutine test_material_padded(error)
      !! Abaqus passes CMNAME as CHARACTER*80 space-padded on the right.
      type(error_type), allocatable, intent(out) :: error
      character(len=80) :: cmname
      cmname = "MCSS                                                                            "
      call check(error, trim(material_name(cmname)), "MCSS", &
                 more="material_name: should handle space-padded cmname")
   end subroutine test_material_padded

   subroutine test_integrator_padded(error)
      !! Abaqus passes CMNAME as CHARACTER*80 space-padded on the right.
      type(error_type), allocatable, intent(out) :: error
      character(len=80) :: cmname
      cmname = "MCSS_INTEG_EULER                                                                "
      call check(error, trim(integrator_name(cmname)), "euler", &
                 more="integrator_name: should handle space-padded cmname")
   end subroutine test_integrator_padded

end module mod_test_cmname_parser_suite
