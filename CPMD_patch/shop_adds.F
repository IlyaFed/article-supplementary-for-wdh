C     ==================================================================
C     == ADDITIONAL ROUTINES FOR SURFACE HOPPING                      ==
C     ==================================================================

C     ==================================================================
      SUBROUTINE WRITE_SHMD(IWAHL,IUNIT,INFI,TEMPP)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      INCLUDE 'system.h'
      INCLUDE 'parac.inc'
      INCLUDE 'shop.inc'
      INCLUDE 'shop_rest.inc'

      INTEGER IWAHL,IUNIT
      INTEGER INFI

      REAL*8  TEMPP

      IF (IWAHL.EQ.0) THEN
      IF(IO_PARENT)
     &WRITE(IUNIT,'(5A)') '#___STEP______T[K]',
     .                       '_________E(1)[a.u.]_________E(2)[a.u.]',
     .                       '___ISH_________DE____________PROB',
     .                       '____________RAND_________________NORM',
     .                       '_'

      ELSEIF (IWAHL.EQ.1) THEN
      IF(IO_PARENT)
     &WRITE(IUNIT,'(1X,I7,F10.3,2(1X,F18.10),1X,I5,
     .                 1X,F10.5, 2(1X,D15.7),1X,F20.15)')
     .                     INFI,TEMPP,EC(1),EC(2),IOLDSURF,EC(2)-EC(1),
     .                                     PROB(IOLDSURF),ZZ,D1SQ+D2SQ

      ELSEIF (IWAHL.EQ.-1) THEN
      IF(IO_PARENT)
     &WRITE(IUNIT,'(5A)') '#-----------------',
     .                       '--------------------------------------',
     .                       '---------------------------------',
     .                       '-------------------------------------',
     .                       '-'

      ELSE
      IF(IO_PARENT)
     &WRITE(6,'(A,I5)') ' INVALID OPTION:',IWAHL
        CALL STOPGM('WRITE_SHMD',' ')

      ENDIF

      RETURN
      END

C     ==================================================================
      SUBROUTINE STATE_SELECT(TAG)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      INCLUDE 'system.h'
      INCLUDE 'parac.inc'
      INCLUDE 'shop.inc'
      INCLUDE 'elct.inc'
      INCLUDE 'spin.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'sfac.inc'
      CHARACTER*(*) TAG
      DIMENSION F(*)
      INTEGER I
      REAL*8  FNL(*),DFNL(*)
      REAL*8 XSTATES,XSNOW,XSAIM

c normally commented out !      
c      CALL MY_SYNC(ALLGRP)
C     ==--------------------------------------------------------------==
      IF (TAG.EQ."S0" ) THEN
         if ( parent ) then
      IF(IO_PARENT)
     &WRITE(6,'(/,'' =='',5X,A16,39X,''==''/)') 'STATE SELECT: S0'
         endif
         TLSE=.FALSE.
         IF ( TLSD0 ) THEN
            TLSD=.TRUE.
            NLSD=2
            NLSX=3
            N=NST_S0
         ELSE
            NLSD=1
            NLSX=1
            N=NST_S0
         ENDIF
C-PARALLEL
         XSTATES=DBLE(N)
         XSNOW=0.0D0
         DO I=NPROC,1,-1
            XSAIM = XSNOW + XSTATES/NPROC
            NST12(I-1,1)=NINT(XSNOW)+1
            NST12(I-1,2)=NINT(XSAIM)
            IF ( NINT(XSAIM).GT.N ) NST12(I-1,2)=N
            IF ( I.EQ.1 ) NST12(I-1,2)=N
            XSNOW = XSAIM
         ENDDO        
         NORBPE=NST12(MEPOS,2)-NST12(MEPOS,1)+1 
c      WRITE(6,*)'S0',me,NORBPE,MEPOS,NST12(MEPOS,2),NST12(MEPOS,1)
C-ENDPARALLEL
         CALL DCOPY(NST_S0,FS0,1,F,1)
      ELSEIF(TAG.EQ."S1") THEN
C...    SELECT S1 
         if ( parent ) then
      IF(IO_PARENT)
     &WRITE(6,'(/,'' =='',5X,A16,39X,''==''/)') 'STATE SELECT: S1'
         endif
         TLSD=.FALSE.
         TLSE=.TRUE.
         TROKS=.TRUE.
         NLSD=4
         NLSX=3
         N=NST_S1
C-PARALLEL---------------------------------
         XSTATES=DBLE(N)
         XSNOW=0.0D0
         DO I=NPROC,1,-1
            XSAIM = XSNOW + XSTATES/NPROC
            NST12(I-1,1)=NINT(XSNOW)+1
            NST12(I-1,2)=NINT(XSAIM)
            IF(NINT(XSAIM).GT.N) NST12(I-1,2)=N
            IF(I.EQ.1) NST12(I-1,2)=N
            XSNOW = XSAIM
         ENDDO  
         NORBPE=NST12(MEPOS,2)-NST12(MEPOS,1)+1
c        WRITE(6,*)'S1',me,NORBPE,MEPOS,NST12(MEPOS,2),NST12(MEPOS,1)
C-ENDPARALLEL   
cMcB   copy NST_S0-NST_S1 ZEROS, TOO !
         CALL DCOPY(NST_S0,FS1,1,F,1)
c        CALL DCOPY(NST_S1,FS1,1,F,1)
cMcB
      ELSE
         CALL STOPGM('STATE_SELECT','UNKNOWN TAG')
      ENDIF 
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================

C     ==================================================================
      SUBROUTINE S0_S1_OVERLAP(C0,CM,DET,C12,C21)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      INCLUDE 'system.h'
      INCLUDE 'parac.inc'
#if   defined(__MCB)  
      COMPLEX*16 C0(NGW,*),CM(NGW,*)
      REAL*8     DET,C12,C21
#else
      INCLUDE 'shop.inc'
      INCLUDE 'spin.inc'
      INCLUDE 'geq0.inc'
      INTEGER ISTATE,JSTATE,IISTATE,JJSTATE,JUSTATE,JDSTATE,IG,INFO
      INTEGER IUSTATE,IDSTATE
      COMPLEX*16 C0(NGW,*),CM(NGW,*)
      REAL*8 UM(NELB2,NELB2),DM(NELB2,NELB2),UM2(NELB2,NELB2),
     .       DM2(NELB2,NELB2),DET,DET2,WORK(NELB2,NELB2),DETUM,DETDM,
     .       C12,C21
      REAL*8 UDOTM(NELB2,NELB2),DDOTM(NELB2,NELB2),UDOTM21(NELB2,NELB2),
     .       DDOTM21(NELB2,NELB2)
      REAL*8 U2DOTM(NELB2,NELB2),D2DOTM(NELB2,NELB2),
     .       U2DOTM21(NELB2,NELB2),D2DOTM21(NELB2,NELB2)
      POINTER (IP_UM,UM),(IP_DM,DM),(IP_UDOTM,UDOTM),(IP_DDOTM,DDOTM)
      POINTER (IP_UM2,UM2),(IP_DM2,DM2),(IP_U2DOTM,U2DOTM),
     .        (IP_D2DOTM,D2DOTM)
      POINTER (IP_UDOTM21,UDOTM21),(IP_DDOTM21,DDOTM21),
     .        (IP_U2DOTM21,U2DOTM21),(IP_D2DOTM21,D2DOTM21)
      POINTER (IP_WORK,WORK)
C     ==--------------------------------------------------------------==
      CALL MEMORY(IP_UM,NELB2*NELB2,'UM')
      CALL MEMORY(IP_DM,NELB2*NELB2,'DM')
      CALL MEMORY(IP_UM2,NELB2*NELB2,'UM2')
      CALL MEMORY(IP_DM2,NELB2*NELB2,'DM2')
      CALL MEMORY(IP_UDOTM,NELB2*NELB2,'UDOTM')
      CALL MEMORY(IP_DDOTM,NELB2*NELB2,'DDOTM')
      CALL MEMORY(IP_U2DOTM,NELB2*NELB2,'U2DOTM')
      CALL MEMORY(IP_D2DOTM,NELB2*NELB2,'D2DOTM')
      CALL MEMORY(IP_UDOTM21,NELB2*NELB2,'UDOTM21')
      CALL MEMORY(IP_DDOTM21,NELB2*NELB2,'DDOTM21')
      CALL MEMORY(IP_U2DOTM21,NELB2*NELB2,'U2DOTM21')
      CALL MEMORY(IP_D2DOTM21,NELB2*NELB2,'D2DOTM21')
      CALL MEMORY(IP_WORK,NELB2*NELB2,'WORK')

      IF(TLSD0)THEN

C        if ( parent ) WRITE(6,*)'IN OVERLAP',NELB2

         CALL AZZERO(UM,NELB2*NELB2)
         CALL AZZERO(DM,NELB2*NELB2)
         CALL AZZERO(UM2,NELB2*NELB2)
         CALL AZZERO(DM2,NELB2*NELB2)
         CALL AZZERO(UDOTM,NELB2*NELB2)
         CALL AZZERO(DDOTM,NELB2*NELB2)
         CALL AZZERO(U2DOTM,NELB2*NELB2)
         CALL AZZERO(D2DOTM,NELB2*NELB2)
         CALL AZZERO(UDOTM21,NELB2*NELB2)
         CALL AZZERO(DDOTM21,NELB2*NELB2)
         CALL AZZERO(U2DOTM21,NELB2*NELB2)
         CALL AZZERO(D2DOTM21,NELB2*NELB2)

         DO ISTATE=1,NELB2
            IISTATE=NELB2+ISTATE
            IUSTATE=NST_S0+ISTATE
            IDSTATE=IUSTATE
            IF(ISTATE.EQ.NELB2)IDSTATE=IUSTATE+1

            DO JSTATE=1,NELB2
               JJSTATE=NELB2+JSTATE
               JUSTATE=JSTATE+NST_S0
               JDSTATE=JUSTATE
               IF(JSTATE.EQ.NELB2)JDSTATE=JUSTATE+1
c     
               DO IG=1,NGW
                  IF (IG.EQ.1.AND.GEQ0) THEN
c...<S0|m1>
                     UM(ISTATE,JSTATE)=UM(ISTATE,JSTATE) + 
     .                                            DCONJG(C0(IG,ISTATE))*
     .                                                   C0(IG,JUSTATE)

                     DM(ISTATE,JSTATE)=DM(ISTATE,JSTATE) + 
     .                                           DCONJG(C0(IG,IISTATE))*
     .                                                  C0(IG,JDSTATE)
c...<S0|d/dt|m1>
                     UDOTM(ISTATE,JSTATE)=UDOTM(ISTATE,JSTATE) +
     .                                            DCONJG(C0(IG,ISTATE))*
     .                                                   CM(IG,JUSTATE)

                     DDOTM(ISTATE,JSTATE)=DDOTM(ISTATE,JSTATE) +
     .                                           DCONJG(C0(IG,IISTATE))*
     .                                                  CM(IG,JDSTATE)
c...<m1|d/dt|S0>
                     UDOTM21(ISTATE,JSTATE)=UDOTM21(ISTATE,JSTATE) +
     .                                           DCONJG(C0(IG,IUSTATE))*
     .                                                    CM(IG,JSTATE)

                     DDOTM21(ISTATE,JSTATE)=DDOTM21(ISTATE,JSTATE) +
     .                                           DCONJG(C0(IG,IDSTATE))*
     .                                                   CM(IG,JJSTATE)
c...<S0|m2>
                     UM2(ISTATE,JSTATE)=UM2(ISTATE,JSTATE) + 
     .                                            DCONJG(C0(IG,ISTATE))*
     .                                                   C0(IG,JDSTATE)

                     DM2(ISTATE,JSTATE)=DM2(ISTATE,JSTATE) + 
     .                                           DCONJG(C0(IG,IISTATE))*
     .                                                  C0(IG,JUSTATE)
c...<S0|d/dt|m2>
                     U2DOTM(ISTATE,JSTATE)=U2DOTM(ISTATE,JSTATE) +
     .                                           DCONJG(C0(IG,IISTATE))*
     .                                                  CM(IG,JDSTATE)

                     D2DOTM(ISTATE,JSTATE)=D2DOTM(ISTATE,JSTATE) +
     .                                           DCONJG(C0(IG,IISTATE))*
     .                                                  CM(IG,JUSTATE)
c...<m2|d/dt|S0>
                     U2DOTM21(ISTATE,JSTATE)=U2DOTM21(ISTATE,JSTATE) +
     .                                           DCONJG(C0(IG,IDSTATE))*
     .                                                  CM(IG,JJSTATE)

                     D2DOTM21(ISTATE,JSTATE)=D2DOTM21(ISTATE,JSTATE) +
     .                                           DCONJG(C0(IG,IUSTATE))*
     .                                                  CM(IG,JJSTATE)
                  ELSE
c...<S0|m1>
                     UM(ISTATE,JSTATE)=UM(ISTATE,JSTATE) +
     .                                            DCONJG(C0(IG,ISTATE))*
     .                                                   C0(IG,JUSTATE)+
     .                                                   C0(IG,ISTATE)*
     .                                            DCONJG(C0(IG,JUSTATE))

                     DM(ISTATE,JSTATE)=DM(ISTATE,JSTATE) +
     .                                           DCONJG(C0(IG,IISTATE))*
     .                                                  C0(IG,JDSTATE)+
     .                                                   C0(IG,IISTATE)*
     .                                            DCONJG(C0(IG,JDSTATE))
c...<S0|d/dt|m1>
                     UDOTM(ISTATE,JSTATE)=UDOTM(ISTATE,JSTATE) +
     .                                            DCONJG(C0(IG,ISTATE))*
     .                                                   CM(IG,JUSTATE)+
     .                                                   C0(IG,ISTATE)*
     .                                            DCONJG(CM(IG,JUSTATE))

                     DDOTM(ISTATE,JSTATE)=DDOTM(ISTATE,JSTATE) +
     .                                           DCONJG(C0(IG,IISTATE))*
     .                                                  CM(IG,JDSTATE)+
     .                                                   C0(IG,IISTATE)*
     .                                            DCONJG(CM(IG,JDSTATE))
c...<m1|d/dt|S0>
                     UDOTM21(ISTATE,JSTATE)=UDOTM21(ISTATE,JSTATE) +
     .                                           DCONJG(C0(IG,IUSTATE))*
     .                                                  CM(IG,JSTATE)+
     .                                                   C0(IG,IUSTATE)*
     .                                            DCONJG(CM(IG,JSTATE))

                     DDOTM21(ISTATE,JSTATE)=DDOTM21(ISTATE,JSTATE) +
     .                                           DCONJG(C0(IG,IDSTATE))*
     .                                                  CM(IG,JJSTATE)+
     .                                                   C0(IG,IDSTATE)*
     .                                            DCONJG(CM(IG,JJSTATE))
c...<S0|m2>
                     UM2(ISTATE,JSTATE)=UM2(ISTATE,JSTATE) + 
     .                                            DCONJG(C0(IG,ISTATE))*
     .                                                   C0(IG,JDSTATE)+
     .                                                   C0(IG,ISTATE)*
     .                                            DCONJG(C0(IG,JDSTATE))


                     DM2(ISTATE,JSTATE)=DM2(ISTATE,JSTATE) + 
     .                                           DCONJG(C0(IG,IISTATE))*
     .                                                  C0(IG,JUSTATE)+
     .                                                   C0(IG,IISTATE)*
     .                                            DCONJG(C0(IG,JUSTATE))
c...<S0|d/dt|m2>
                     U2DOTM(ISTATE,JSTATE)=U2DOTM(ISTATE,JSTATE) +
     .                                           DCONJG(C0(IG,IISTATE))*
     .                                                  CM(IG,JDSTATE)+
     .                                                  C0(IG,IISTATE)*
     .                                           DCONJG(CM(IG,JDSTATE))

                     D2DOTM(ISTATE,JSTATE)=D2DOTM(ISTATE,JSTATE) +
     .                                           DCONJG(C0(IG,IISTATE))*
     .                                                  CM(IG,JUSTATE)+
     .                                                   C0(IG,IISTATE)*
     .                                            DCONJG(CM(IG,JUSTATE))
c...<m2|d/dt|S0>
                     U2DOTM21(ISTATE,JSTATE)=U2DOTM21(ISTATE,JSTATE) +
     .                                           DCONJG(C0(IG,IDSTATE))*
     .                                                  CM(IG,JJSTATE)+
     .                                                   C0(IG,IDSTATE)*
     .                                            DCONJG(CM(IG,JJSTATE))

                     D2DOTM21(ISTATE,JSTATE)=D2DOTM21(ISTATE,JSTATE) +
     .                                           DCONJG(C0(IG,IUSTATE))*
     .                                                  CM(IG,JJSTATE)+
     .                                                   C0(IG,IUSTATE)*
     .                                            DCONJG(CM(IG,JJSTATE))

                  END IF
               END DO
C              if ( parent ) then 
C                  WRITE(6,*)'UM',istate,jstate,udotm(istate,jstate)
C     $                  ,u2dotm(istate,jstate)
C                  WRITE(6,*)'DM',istate,jstate,ddotm(istate,jstate)
C     $                  ,d2dotm(istate,jstate)
C              endif
            END DO
         END DO
         
         CALL MY_SUM_D(UM,NELB2*NELB2,ALLGRP)
         CALL MY_SUM_D(DM,NELB2*NELB2,ALLGRP)
         CALL MY_SUM_D(UM2,NELB2*NELB2,ALLGRP)
         CALL MY_SUM_D(DM2,NELB2*NELB2,ALLGRP)
         CALL MY_SUM_D(UDOTM,NELB2*NELB2,ALLGRP)
         CALL MY_SUM_D(DDOTM,NELB2*NELB2,ALLGRP)
         CALL MY_SUM_D(U2DOTM,NELB2*NELB2,ALLGRP)
         CALL MY_SUM_D(D2DOTM,NELB2*NELB2,ALLGRP)
         CALL MY_SUM_D(UDOTM21,NELB2*NELB2,ALLGRP)
         CALL MY_SUM_D(DDOTM21,NELB2*NELB2,ALLGRP)
         CALL MY_SUM_D(U2DOTM21,NELB2*NELB2,ALLGRP)
         CALL MY_SUM_D(D2DOTM21,NELB2*NELB2,ALLGRP)

CNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
c...<S0|m1>
         CALL DET_AND_INV(NELB2,UM,WORK,INFO,DETUM)
         CALL DET_AND_INV(NELB2,DM,WORK,INFO,DETDM)
         DET=DETUM*DETDM/SQRT(2.0D0)
C         WRITE(6,*)'DET 1',det
C.................................................................
         C12=0.0d0
         C21=0.0d0

         DO ISTATE=1,NELB2
            DO JSTATE=1,NELB2
               C12=C12+DET*( UDOTM(ISTATE,JSTATE)*UM(JSTATE,ISTATE)+ 
     .                       DDOTM(ISTATE,JSTATE)*DM(JSTATE,ISTATE) )
               C21=C21+DET*( UDOTM21(ISTATE,JSTATE)*UM(ISTATE,JSTATE)+
     .                       DDOTM21(ISTATE,JSTATE)*DM(ISTATE,JSTATE) )
            ENDDO
         ENDDO
      IF (IO_PARENT)
     &WRITE(6,*)'D_21, D_12 [1]',C21,C12 
c...<S0|m2>
         CALL DET_AND_INV(NELB2,UM2,WORK,INFO,DETUM)
         CALL DET_AND_INV(NELB2,DM2,WORK,INFO,DETDM)
         DET2=DETUM*DETDM/SQRT(2.0D0)
C         WRITE(6,*)'DET 2',det2
C.................................................................

         DO ISTATE=1,NELB2
            DO JSTATE=1,NELB2
               C12=C12+DET2*( U2DOTM(ISTATE,JSTATE)*UM2(JSTATE,ISTATE)+ 
     .                        D2DOTM(ISTATE,JSTATE)*DM2(JSTATE,ISTATE) )
               C21=C21+DET2*(
     .                      U2DOTM21(ISTATE,JSTATE)*UM2(ISTATE,JSTATE)+
     .                      D2DOTM21(ISTATE,JSTATE)*DM2(ISTATE,JSTATE) )
            ENDDO
         ENDDO
      IF (IO_PARENT)
     &WRITE(6,*)'D_21, D_12 [2]',C21,C12     
         DET=DET+DET2
CNNNNNN UM AND DM ARE NOW THE INVERSE MATRICES !!!!!!!!!!! 
C         STOP
      ELSE
      IF (IO_PARENT)
     &WRITE(6,*)'ACHTUNG, KEIN LSD!!'

         CALL AZZERO(UM,NST_S0*NST_S0)
         CALL AZZERO(DM,NST_S0*NST_S0)
         CALL AZZERO(UDOTM,NST_S0*NST_S0)
         CALL AZZERO(DDOTM,NST_S0*NST_S0)
         CALL AZZERO(UDOTM21,NST_S0*NST_S0)
         CALL AZZERO(DDOTM21,NST_S0*NST_S0)

         DO ISTATE=1,NST_S0
            IUSTATE=NST_S0+ISTATE
            IDSTATE=IUSTATE
            IF ( ISTATE.EQ.NELB2 ) IDSTATE=IUSTATE+1

            DO JSTATE=1,NST_S0
               JUSTATE=JSTATE+NST_S0
               JDSTATE=JUSTATE
               IF ( JSTATE.EQ.NELB2 ) JDSTATE=JUSTATE+1
      
               DO IG=1,NGW
                  IF ( (IG.EQ.1).AND.GEQ0 ) THEN

                     UM(ISTATE,JSTATE)=UM(ISTATE,JSTATE) +
     .                                            DCONJG(C0(IG,ISTATE))*
     .                                                   C0(IG,JUSTATE)

                     DM(ISTATE,JSTATE)=DM(ISTATE,JSTATE) +
     .                                            DCONJG(C0(IG,ISTATE))*
     .                                                   C0(IG,JDSTATE)

                     UDOTM(ISTATE,JSTATE)=UDOTM(ISTATE,JSTATE) +
     .                                            DCONJG(C0(IG,ISTATE))*
     .                                                   CM(IG,JUSTATE)

                     DDOTM(ISTATE,JSTATE)=DDOTM(ISTATE,JSTATE) +
     .                                            DCONJG(C0(IG,ISTATE))*
     .                                                   CM(IG,JDSTATE)

                     UDOTM21(ISTATE,JSTATE)=UDOTM21(ISTATE,JSTATE) +
     .                                           DCONJG(C0(IG,IUSTATE))*
     .                                                    CM(IG,JSTATE)

                     DDOTM21(ISTATE,JSTATE)=DDOTM21(ISTATE,JSTATE) +
     .                                           DCONJG(C0(IG,IDSTATE))*
     .                                                    CM(IG,JSTATE)
                  ELSE

                     UM(ISTATE,JSTATE)=UM(ISTATE,JSTATE) +
     .                                            DCONJG(C0(IG,ISTATE))*
     .                                                   C0(IG,JUSTATE)+
     .                                                   C0(IG,ISTATE)*
     .                                            DCONJG(C0(IG,JUSTATE))

                     DM(ISTATE,JSTATE)=DM(ISTATE,JSTATE) +
     .                                            DCONJG(C0(IG,ISTATE))*
     .                                                   C0(IG,JDSTATE)+
     .                                                   C0(IG,ISTATE)*
     .                                            DCONJG(C0(IG,JDSTATE))
                     
                     UDOTM(ISTATE,JSTATE)=UDOTM(ISTATE,JSTATE) +
     .                                            DCONJG(C0(IG,ISTATE))*
     .                                                   CM(IG,JUSTATE)+
     .                                                   C0(IG,ISTATE)*
     .                                            DCONJG(CM(IG,JUSTATE))
                     
                     DDOTM(ISTATE,JSTATE)=DDOTM(ISTATE,JSTATE) +
     .                                            DCONJG(C0(IG,ISTATE))*
     .                                                   CM(IG,JDSTATE)+
     .                                                   C0(IG,ISTATE)*
     .                                            DCONJG(CM(IG,JDSTATE))
                     
                     UDOTM21(ISTATE,JSTATE)=UDOTM21(ISTATE,JSTATE) +
     .                                           DCONJG(C0(IG,IUSTATE))*
     .                                                    CM(IG,JSTATE)+
     .                                                   C0(IG,IUSTATE)*
     .                                             DCONJG(CM(IG,JSTATE))

                     DDOTM21(ISTATE,JSTATE)=DDOTM21(ISTATE,JSTATE) +
     .                                           DCONJG(C0(IG,IDSTATE))*
     .                                                  CM(IG,JSTATE)  +
     .                                                  C0(IG,IDSTATE)*
     .                                           DCONJG(CM(IG,JSTATE))
                     
                  ENDIF
               ENDDO

            ENDDO
         ENDDO

         CALL MY_SUM_D(UM,NST_S0*NST_S0,ALLGRP)
         CALL MY_SUM_D(DM,NST_S0*NST_S0,ALLGRP)
         CALL MY_SUM_D(UDOTM,NST_S0*NST_S0,ALLGRP)
         CALL MY_SUM_D(DDOTM,NST_S0*NST_S0,ALLGRP)
         CALL MY_SUM_D(UDOTM21,NST_S0*NST_S0,ALLGRP)
         CALL MY_SUM_D(DDOTM21,NST_S0*NST_S0,ALLGRP)
CNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
         CALL DET_AND_INV(NST_S0,UM,WORK,INFO,DETUM)
         CALL DET_AND_INV(NST_S0,DM,WORK,INFO,DETDM)
CNNNNNN UM AND DM ARE NOW THE INVERSE MATRICES !!!!!!!!!!! 
         DET=SQRT(2.0D0)*DETUM*DETDM
C.................................................................
         C12=0.0D0
         C21=0.0D0
         DO ISTATE=1,NELB2
            DO JSTATE=1,NELB2
               C12=C12+DET*( UDOTM(ISTATE,JSTATE)*UM(JSTATE,ISTATE)+ 
     .                       DDOTM(ISTATE,JSTATE)*DM(JSTATE,ISTATE) )
               C21=C21+DET*( UDOTM21(ISTATE,JSTATE)*UM(ISTATE,JSTATE)+
     .                       DDOTM21(ISTATE,JSTATE)*DM(ISTATE,JSTATE) )
            ENDDO
         ENDDO
      IF (IO_PARENT)
     &WRITE(6,*)'D_21, D_12',C21,C12         
      ENDIF

C...FREE MEMORY
      CALL FREEM(IP_UM)
      CALL FREEM(IP_DM)
      CALL FREEM(IP_UM2)
      CALL FREEM(IP_DM2)
      CALL FREEM(IP_UDOTM)
      CALL FREEM(IP_DDOTM)
      CALL FREEM(IP_U2DOTM)
      CALL FREEM(IP_D2DOTM)
      CALL FREEM(IP_UDOTM21)
      CALL FREEM(IP_DDOTM21)
      CALL FREEM(IP_U2DOTM21)
      CALL FREEM(IP_D2DOTM21)
      CALL FREEM(IP_WORK)
#endif
C     ==--------------------------------------------------------------==
       RETURN
       END
C     ==================================================================

C     ==================================================================
      SUBROUTINE DET_AND_INV(N,A,B,INFO,DET)
C     ==--------------------------------------------------------------==
C     == Inverse Matrix A(N,N)                                        ==
C     ==--------------------------------------------------------------==
C     == INPUT:  N Dimension                                          ==
C     ==         A(N,N) Matrix                                        ==
C     == OUTPUT: A(N,N) inverse matrix of A                           ==
C     ==         B is a work array                                    ==
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      INTEGER N,INFO,I,IPIVOT(N)
      REAL*8  A(N,N),B(N,N),DET,WORK(N),SIGN
C     Variables
      INTEGER LB
C     ==--------------------------------------------------------------==
      IF(N.EQ.1) THEN
        A(1,1)=1.D0/A(1,1)
      ELSE
        LB=(N-1)*N
C       Compute an LU factorization
c           WRITE(6,*)'A',A(1,1),A(1,2)
c           WRITE(6,*)'A',A(2,1),A(2,2)
        CALL DGETRF(N,N,A,N,IPIVOT,INFO)
c        WRITE(6,*)'INFO',INFO
        DET=1.0D0
        DO I=1,N
        SIGN=1.0D0
           IF(IPIVOT(I).NE.I)SIGN=-1.0D0
c           WRITE(6,*)'A(I,I)',A(I,I),IPIVOT(I)
           DET=DET*SIGN*A(I,I)
C           WRITE(6,*)'det',DET,SIGN
        END DO
        IF(INFO.EQ.0) THEN
C         Compute the inverse of a matrix using the LU factorization
          CALL DGETRI(N,A,N,IPIVOT,WORK,N,INFO)
        ENDIF
      ENDIF
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================

C     ==================================================================
      SUBROUTINE DECIDE(E,FION0,FION1,JUMP)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      INCLUDE 'system.h'
      INCLUDE 'parac.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'coor.inc'
      INCLUDE 'ropt.inc'
      INCLUDE 'shop.inc'
      INCLUDE 'shop_rest.inc'
      INCLUDE 'shop_ekinqm.inc'
      INCLUDE 'mm_dim.inc'
cMcB
      INCLUDE 'ions.inc'
cMcB
      INCLUDE 'fileopen.inc'

c      INCLUDE 'sfac.inc'
c     REAL*8     FNL(*),DFNL(*)
      REAL*8 FION0(3,NAX,*), FION1(3,NAX,*) 
      REAL*8 E(2),PPP(2),DE,SKAL,EKINP,ARG,EKINSH
      INTEGER I,J,K,JSURF,TMPSURF,IA,IS
      LOGICAL JUMP,FERROR
      CHARACTER*2 TAG(2)
c     COMMON /PROB1/D1SQ,D2SQ,D11DOT,D22DOT 
c     REAL*8 PROB(2)
c     COMMON /PROB2/ PROB,ZZ
      COMMON /EKIN/EKINP
CMcB 
      REAL*8     DELT
CMcB 
CMcB
      DELT=DELT_IONS
c     DELT=DELT_ELEC
CMcB
      TAG(1)='S0'
      TAG(2)='S1'
      IF(ISURF.EQ.1) JSURF=2
      IF(ISURF.EQ.2) JSURF=1
      IOLDSURF=ISURF
      IF(.NOT.JUMP) THEN
         CALL STATE_SELECT(TAG(ISURF))
         ETOT=E(ISURF)
         IF(ISURF.EQ.2) THEN
            DO I=1,3
               DO J=1,NAX
                  DO K=1,NSX
                     FION(I,J,K)= FION1(I,J,K)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            DO I=1,3
               DO J=1,NAX
                  DO K=1,NSX
                     FION(I,J,K)= FION0(I,J,K)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ELSE
CNNNN IF SURFACE JUMP ALLOWED:
      PPP(2)=-DELT*D22DOT/D2SQ
      PPP(1)=-DELT*D11DOT/D1SQ
      PROB(1)=PPP(1)
      PROB(2)=PPP(2)
C---- WRITE COEFFICIENTS TO FILE --------------------------
c     if (parent) then
c       call fileopen(789,'prob.dat',FO_APP,FERROR)
c       write(789,12)infi,d1sq,d2sq
c12     format(i6,4(1x,f12.8))     
c       call fileclose(789)
c     end if
C--------------------------------------------------------      
      CALL SHOP_GRND(ZZ)
cMcB
c     ZZ=ZZ*1.0D-5
cMcB
      IF (IO_PARENT)
     &WRITE(6,*)'tot',d1sq+d2sq
c       IF(PARENT) write(6,*)'d1,d2',d1sq,d2sq
c       IF(PARENT) write(6,*)'d11dot,d22dot',d11dot,d22dot
c     IF(PARENT) WRITE(6,*)'ZZ',ZZ,'PROB',P(ISURF),ISURF
      IF (IO_PARENT)
     &WRITE(6,*)'ZZ',ZZ,'PROB',PPP(ISURF),ISURF
CNNNN IF RANDOM NUMBER SMALLER THAN TRANSITION PROBABILITY:
c     IF (ZZ.LT.P(ISURF)) THEN
      IF (ZZ.LT.PPP(ISURF)) THEN
        IF (TQMMM) THEN
           EKINSH=EKINQMSH
        ELSE
           EKINSH=EKINP
        ENDIF
        ARG=(E(ISURF)-E(JSURF))/EKINSH+1.0D0
        IF(ARG.GT.0.0D0) THEN  !classically allowed transition
           CALL STATE_SELECT(TAG(JSURF))
           ETOT=E(JSURF)
           SKAL=SQRT(ARG)
      IF (IO_PARENT)
     &WRITE(6,*)'SKAL',SKAL,EKINSH
           DO I=1,3
c           scale QM centres only!
c             DO J=1,NAX
c                DO K=1,NSX
c my_code begin
c             DO IS=1,NAX
c                DO IA=1,NSX
c this part don't work
c my_code end
c my_code begin
              DO IS=1,NSP
                 DO IA=1,NA(IS)
c this part don't test
c my_code end
c             DO IS=1,NSPq
c                DO IA=1,NAq(IS)
                    VELP(I,IA,IS)= SKAL*VELP(I,IA,IS) !rescale velos
                 END DO
              END DO
           END DO  
           IF(JSURF.EQ.2)THEN
           DO I=1,3
              DO J=1,NAX
                 DO K=1,NSX
                    FION(I,J,K)= FION1(I,J,K)
                 END DO
              END DO
           END DO
           ELSE
           DO I=1,3
              DO J=1,NAX
                 DO K=1,NSX
                     FION(I,J,K)= FION0(I,J,K)
                 END DO
              END DO
           END DO
           ENDIF
           ISURF=JSURF
        ELSE ! classically forbidden transition (arg<0)
CNNNN IF RANDOM NUMBER NOT SMALLER THAN TRANSITION PROBABILITY:
        CALL STATE_SELECT(TAG(ISURF))
        ETOT=E(ISURF)
            IF(ISURF.EQ.2)THEN !isurf=2
            DO I=1,3
               DO J=1,NAX
                  DO K=1,NSX
                  FION(I,J,K)= FION1(I,J,K)
                  END DO
               END DO
            END DO
            ELSE
            DO I=1,3
               DO J=1,NAX
                  DO K=1,NSX
                  FION(I,J,K)= FION0(I,J,K)
                  END DO
               END DO
            END DO
            ENDIF              !isurf=2     
        ENDIF               !arg<0
      ENDIF                 !zz>p
      ENDIF                 !jump=true
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================

C     common block init for random number generator
      BLOCK DATA shop_rng_init
      INTEGER N,N1
      PARAMETER (N=624,N1=N+1)
      INTEGER mti,mt(0:N-1)
      COMMON /rng_block/mti,mt

C     mti==N+1 means mt[N] is not initialized
      DATA mti/N1/
      END

C     ==================================================================
      SUBROUTINE shop_grnd(RNV)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      INCLUDE 'system.h'
      INCLUDE 'parac.inc'
      INCLUDE 'fileopen.inc'
C     arguments
      REAL*8 RNV

C     Period parameters
      INTEGER N,N1,M,MATA,UMASK,LMASK
      PARAMETER (N=624, N1=N+1, M=397, MATA=-1727483681)
      PARAMETER (UMASK=-2147483647, LMASK=2147483647)

C     Tempering parameters
      INTEGER TMASKB,TMASKC
      PARAMETER(TMASKB= -1658038656, TMASKC=-272236544 )

C     array for the state vector
      INTEGER mti,mt(0:N-1)
      common /rng_block/mti,mt
      save   /rng_block/

      INTEGER  mag01(0:1)
      data mag01/0, MATA/
      save mag01

c     Locals
      INTEGER y,kk,iunit
      LOGICAL FERROR
      CHARACTER*15 rng_state_file

C-----------------------------------------------------------------------
C     generate N words at one time
      IF(mti.ge.N) THEN

C     initialize or restore random number generator on first call.
        IF(mti.eq.N+1) THEN
          CALL shop_sgrnd
        ENDIF
C
        DO kk=0,N-M-1
            y=IOR(IAND(mt(kk),UMASK),IAND(mt(kk+1),LMASK))
            mt(kk)=IEOR(IEOR(mt(kk+M),ISHFT(y,-1)),mag01(IAND(y,1)))
        enddo

        do kk=N-M,N-2
            y=IOR(IAND(mt(kk),UMASK),IAND(mt(kk+1),LMASK))
            mt(kk)=IEOR(IEOR(mt(kk+(M-N)),ISHFT(y,-1)),mag01(IAND(y,1)))
        ENDDO

        y=IOR(IAND(mt(N-1),UMASK),IAND(mt(0),LMASK))
        mt(N-1)=IEOR(IEOR(mt(M-1),ISHFT(y,-1)),mag01(IAND(y,1)))
        mti = 0
      ENDIF
C
      y=mt(mti)
      mti=mti+1
C
      y=IEOR(y,ISHFT(y,-11))
      y=IEOR(y,IAND(ISHFT(y,7),TMASKB))
      y=IEOR(y,IAND(ISHFT(y,15),TMASKC))
      y=IEOR(y,ISHFT(y,-18))
C
      IF (y.LT.0) THEN
        RNV=(DBLE(y)+2.0d0**32)/(2.0d0**32-1.0d0)
      ELSE
        RNV=DBLE(y)/(2.0d0**32-1.0d0)
      ENDIF

      IF(PARENT) THEN

C     record the state of the random number generator in the file
C     SHOP_RNG_STATE. this should not pose a performance problem,
C     since this is only done once per MD step.
        FERROR=.FALSE.
        iunit=69
        rng_state_file='SHOP_RNG_STATE'
        IF(IO_PARENT)
     &CALL FILEOPEN(iunit,rng_state_file,FO_DEF,FERROR)
        IF(.NOT.FERROR) THEN
      IF(IO_PARENT)
     &WRITE(iunit,'(I14)')   mti
      IF(IO_PARENT)
     &WRITE(iunit,'(5I14)') (mt(kk),kk=0,N-1)
          IF(IO_PARENT)
     &CALL FILECLOSE(iunit)
        ENDIF
      ENDIF
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================

C     ==================================================================
C     seed or restore the random number generator.
C     for surface hopping the seed is supposed to be new on every
C     new run, so we try to read the recorded state of the last run first.
      SUBROUTINE shop_sgrnd
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      INCLUDE 'system.h'
      INCLUDE 'parac.inc'
      INCLUDE 'irat.inc'
      INCLUDE 'fileopen.inc'
C     Period parameter
      INTEGER N
      PARAMETER(N=624)

      LOGICAL FERROR
      REAL*8 TIMEC,TIMEF,REPPRNGU
      INTEGER SEED,MSGLEN,I,iunit
      CHARACTER*15 rng_state_file
C
C     the array for the state vector
      INTEGER mti,mt(0:N-1)
      COMMON /rng_block/mti,mt
      SAVE   /rng_block/
cMcB
      integer mttmp
cMcB

      IF(PARENT) THEN
      IF(IO_PARENT)
     &WRITE(6,200)
        FERROR=.FALSE.
        iunit=69
        rng_state_file='SHOP_RNG_STATE'
        IF(IO_PARENT)
     &CALL FILEOPEN(iunit,rng_state_file,FO_OLD,FERROR)
        IF(.NOT.FERROR) THEN
C     read state vector from restart in a fault tolerant way.
      IF(IO_PARENT)
     &WRITE(6,'(1X,2A,$)')'READING OLD STATE FROM: ',
     &           rng_state_file
          IF(IO_PARENT)
     &READ(iunit,'(I14)',ERR=100,END=100)   mti
          IF(IO_PARENT)
     &READ(iunit,'(5I14)',ERR=100,END=100) (mt(i),i=0,N-1)
          FERROR=.FALSE.
      IF(IO_PARENT)
     &WRITE(6,'('' OK.'',/)')
          GOTO 110
 100      CONTINUE
      IF(IO_PARENT)
     &WRITE(6,'('' ERROR.'',/)')
          FERROR=.TRUE.
 110      CONTINUE
          IF(IO_PARENT)
     &CALL FILECLOSE(iunit)
        ENDIF

C     no previous random number state was available or readable. 
C     setting initial seeds for the mt() array using 
C     a simple random number generator.
        IF(FERROR) THEN
c         SEED=ANINT(DABS(TIMEC()+TIMEF()+REPPRNGU())*1.0D8)+42
cMcB
          SEED=ANINT(DABS(TIMEC()+TIMEF()+REPPRNGU()*1.0D8))+42
c         SEED=ANINT(DABS(TIMEC()+TIMEF()+REPPRNGU())*1.0D2)+42
cMcB

cMcB
c         if (parent) then
c            write (6,*) 'sgrnd> ',' TIMEC() = ',TIMEC()
c            write (6,*) 'sgrnd> ',' TIMEF() = ',TIMEF()
c            write (6,*) 'sgrnd> ',' REPPRNGU() = ',REPPRNGU()
c         endif
cMcB
          IF (PARENT) THEN
      IF(IO_PARENT)
     &WRITE(6,'(1X,A,I12,/)') 'SEEDING RNG WITH SEED:',SEED
          ENDIF
          mt(0)= IAND(SEED,-1)
          DO mti=1,N-1
cMcB
            mttmp=69069 * mt(mti-1)
c           mt(mti) = IAND(69069 * mt(mti-1),-1)
            mt(mti) = IAND(mttmp,-1)

c           if (parent) write (6,*) 'sgrnd> ',mti,mttmp,mt(mti)
cMcB
          ENDDO
cMcB
c         stop 
cMcB
        ELSE
        ENDIF  
      ENDIF

C     distribute the RNG state vector across all nodes
C     so we get the same sequence of numbers on all of them.
      MSGLEN=(N+1) * 8/IRAT     ! mti + mt() array
      CALL MY_BCAST(mti,MSGLEN,SOURCE,ALLGRP)

 200  FORMAT(/' SEEDING RANDOM NUMBER GENERATOR FOR SURFACE HOPPING.',
     &     /,' REF: M. Matsumoto and T. Nishimura,',
     &     /,' ACM Transactions on Modeling and Computer Simulation,',
     &     /,' Vol. 8, No. 1, January 1998, pp 3--30.',/)

C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
