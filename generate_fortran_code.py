
from model_definition import define_model

# constructing the model
model_definition = define_model(2, 2)

# computing the tensor (might take a long time depending on the resolution
model_definition.compute_tensor(numerical=False, compute_inner_products=True, compute_inner_products_kwargs={'timeout': None})

# generating the tendencies in Fortran format
language_translations = {'β': 'beta', 'sqrt(2)': 'sq2', 'χ_': 'x'}
fortran_tendencies, fortran_jacobian = model_definition.compute_tendencies(language='fortran', lang_translation=language_translations)

# writing to file
fortran_file_part1 = """MODULE var

  REAL(KIND=8) :: n, r, beta
"""

fortran_file_part1 += "\n  REAL(KIND=8) :: "
for i in range(model_definition.ndim):  # must still threat case > 10
    fortran_file_part1 += "x"+str(i)
fortran_file_part1 += "\n"

fortran_file_part1 += """  REAL(KIND=8) :: sq2, pi, tdel
  
  parameter(n = 0.20D0)
  parameter(r = 0.010D0)
  parameter(beta = 1.026671749185127D0)
! chi = (0.0 0.1 0.0 0.0 0.1 0.0 0.0 0.0 0.0 0.0)
! Model coming from Python: chi_{i} <---> x{i-1}
"""

fortran_file_part1 += "\n"
for i in range(model_definition.ndim):  # must still threat case > 10
    fortran_file_part1 += f"  parameter(x{i} = 0.D0)\n"

fortran_file_part1 +="""
  parameter(pi = DACOS(-1.D0))
  parameter(sq2 = SQRT(2.D0))
  parameter(tdel = 0.0001D0)
  
END MODULE var


program WTG_TM_HEUN
USE var
IMPLICIT NONE


!-------------------------------------------------------------------------------------------------------|
!                                                                                                       |
!  HEUN INTEGRATION OF A FORCING (MODEL)                                                                |
!                                                                                                       |
!-------------------------------------------------------------------------------------------------------|

! MODEL STATE

  REAL(KIND=8), DIMENSION(10) :: U, F, FF, UU
  INTEGER :: i, itr, nt, klk
  
! PARAMETERS

  itr=1
  nt=20000000
  klk=itr


  
!
!      OPENING OF THE DATA FILE
!

  open(11,file="param.out")       !parameter file
  open(13,file="evol_field.dat")  !evolution of the reference fields



!
!      WRITING PARAMETER TO SCREEN AND TO FILE

  write(6,*) 'Start of the integration'
       
  write(11,*) 'aspect ratio ',n
  write(11,*) 'beta prime ', beta
  write(11,*) 'bottom friction coeff ', r
  write(11,*) 'chi_1 ', x0
  write(11,*) 'chi_2 ', x1
  write(11,*) 'chi_3 ', x2
  write(11,*) 'chi_4 ', x3
  write(11,*) 'chi_5 ', x4
  write(11,*) 'chi_6 ', x5
  write(11,*) 'chi_7 ', x6
  write(11,*) 'chi_8 ', x7
  write(11,*) 'chi_9 ', x8
  write(11,*) 'chi_10 ', x9
  
  write(11,*) 'pi', pi
  write(11,*) 'sq2', sq2
  write(11,*) 'tdel', tdel
  
  
!      INITIAL CONDITIONS OF THE MODEL

  U(1) = 0.00001D0
  write(6,*) 'Initial condition', U

!
!      TIME LOOP
!


     do i=1,nt
        
!       
!       COMPUTING THE FORCING
!
        
        call forcing(U, F) ! computation of the full-coupled forcing

!
!         EULER STEP OF THE INTEGRATION SCHEME
!

        call step(U, F, UU) ! full-coupled step

!
!         SECOND ESTIMATION OF THE DETERMINISTIC TENDENCIES
!

        call forcing(UU, FF) ! computation of the full-coupled forcing

!
!          FINAL SOLUTION WITH AVERAGED TENDENCIES, HEUN METHOD
!

        call step(U,0.5*(F+FF), UU)
        U=UU

!    
!         WRITING THE EVOLUTION AFTER THE TRANSIENT ITR
!

        if (i.ge.itr) then
           if (i.eq.itr) print*, 'Start of the file writing'
           if (klk.eq.i) then
              klk=klk+1000
              write(13,*) (i-itr)*tdel,U
           endif


        endif
     
     enddo
     
     print*, 'End of the integration'

     close(11)
     close(13)

    END program WTG_TM_HEUN

      



    SUBROUTINE forcing(U, F)
      USE var
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(10), INTENT(IN) :: U
      REAL(KIND=8), DIMENSION(10), INTENT(OUT) :: F
"""

fortran_file_part2 = """

    END SUBROUTINE forcing



    SUBROUTINE step(U, F, UU)
      USE var
      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(10), INTENT(IN) :: U
      REAL(KIND=8), DIMENSION(10), INTENT(IN) :: F
      REAL(KIND=8), DIMENSION(10), INTENT(OUT) :: UU

      UU = U+tdel*F
      
    END SUBROUTINE step

"""

fortran_code = fortran_file_part1
for trend in fortran_tendencies[0]:
    fortran_code += '\n\n\t  ' + trend

fortran_code += '\n' + fortran_file_part2

# # Adapting Chi from Python to Fortran:  x{i} --> chi_{i+1}
# for i in range(model_definition.ndim-1, -1, -1):
#     fortran_code = fortran_code.replace('x'+str(i), 'chi_'+str(i+1))
#
with open('WTG-TM_Heun.f90', 'w') as f:
    f.write(fortran_code)
