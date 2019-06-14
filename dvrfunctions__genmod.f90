        !COMPILER-GENERATED INTERFACE MODULE: Fri Jun 14 15:32:18 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DVRFUNCTIONS__genmod
          INTERFACE 
            SUBROUTINE DVRFUNCTIONS(HOLDER,NODES,NODE,I,SZ,VAL,WEIGHTS)
              INTEGER(KIND=4) :: SZ
              REAL(KIND=8) :: HOLDER(SZ)
              REAL(KIND=8) :: NODES(SZ)
              REAL(KIND=8) :: NODE
              INTEGER(KIND=4) :: I
              REAL(KIND=8) :: VAL
              REAL(KIND=8) :: WEIGHTS(SZ)
            END SUBROUTINE DVRFUNCTIONS
          END INTERFACE 
        END MODULE DVRFUNCTIONS__genmod
