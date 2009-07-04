








C
C YPRIMEG.FOR - Gateway function for YPRIME.FOR
C
C This is an example of the FORTRAN code required for interfacing
C a .MEX file to MATLAB.
C
C This subroutine is the main gateway to MATLAB.  When a MEX function
C  is executed MATLAB calls the MEXFUNCTION subroutine in the corresponding
C  MEX file.  
C
C Copyright 1984-2004 The MathWorks, Inc.
C $Revision: 1.9.2.1 $
C
      SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on 64-bit platforms
C
      INTEGER PLHS(*), PRHS(*)

C-----------------------------------------------------------------------
C

      INTEGER NLHS, NRHS
C
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on 64-bit platforms
C
      INTEGER MXCREATEDOUBLEMATRIX, MXGETPR

C-----------------------------------------------------------------------
C

      INTEGER MXGETM, MXGETN
C
C KEEP THE ABOVE SUBROUTINE, ARGUMENT, AND FUNCTION DECLARATIONS FOR USE
C IN ALL YOUR FORTRAN MEX FILES.
C---------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on 64-bit platforms
C
c      INTEGER YPP, TP, YP
      INTEGER xnutp, xprotp, xalfap, xhp, ap, xp
      INTEGER abarnump,abarp,zbarp,yelocalp

C-----------------------------------------------------------------------
C

      INTEGER M, N
c      REAL*8 RYPP(4), RTP, RYP(4)
      real*8 rxnut, rxprot, rxalfa, rxh, ra, rx
      real*8 rabarnum,rabar,rzbar,ryelocal

C
C CHECK FOR PROPER NUMBER OF ARGUMENTS
C
      IF (NRHS .NE. 6) THEN
        CALL MEXERRMSGTXT('requires 6 input arguments')
      ELSEIF (NLHS .NE. 4) THEN
        CALL MEXERRMSGTXT('requires 4 output argument')
      ENDIF
C
C CHECK THE DIMENSIONS OF Y.  IT CAN BE 4 X 1 OR 1 X 4.
C
      do iii=1,NRHS

         M = MXGETM(PRHS(iii))
         N = MXGETN(PRHS(iii))
C
         IF ((M .NE. 1) .OR. (N .NE. 1)) THEN
            CALL MEXERRMSGTXT('All required to be a 1 x 1 vector')
         ENDIF

      end do
C
C CREATE A MATRIX FOR RETURN ARGUMENT
C
c      PLHS(1) = MXCREATEDOUBLEMATRIX(1,1,0)
      do iii=1,NLHS
         PLHS(iii) = MXCREATEDOUBLESCALAR(0.0)
      end do
C
C ASSIGN POINTERS TO THE VARIOUS PARAMETERS
C
c      YPP = MXGETPR(PLHS(1))
      abarnump = MXGETPR(PLHS(1))
      abarp = MXGETPR(PLHS(2))
      zbarp = MXGETPR(PLHS(3))
      yelocalp = MXGETPR(PLHS(4))
C
c mxGetScalar
c      TP = MXGETPR(PRHS(1))
c      YP = MXGETPR(PRHS(2))
      xnutp = MXGETPR(PRHS(1))
      xprotp = MXGETPR(PRHS(2))
      xalfap = MXGETPR(PRHS(3))
      xhp = MXGETPR(PRHS(4))
      ap = MXGETPR(PRHS(5))
      xp = MXGETPR(PRHS(6))
C
C COPY RIGHT HAND ARGUMENTS TO LOCAL ARRAYS OR VARIABLES
      CALL MXCOPYPTRTOREAL8(xnutp, rxnut, 1)
      CALL MXCOPYPTRTOREAL8(xprotp, rxprot, 1)
      CALL MXCOPYPTRTOREAL8(xalfap, rxalfa, 1)
      CALL MXCOPYPTRTOREAL8(xhp, rxh, 1)
      CALL MXCOPYPTRTOREAL8(ap, ra, 1)
      CALL MXCOPYPTRTOREAL8(xp, rx, 1)
C
C DO THE ACTUAL COMPUTATIONS IN A SUBROUTINE
C       CREATED ARRAYS.  
C
      call compute_nuclear_azbar(rxnut, rxprot, rxalfa, rxh, ra, rx, rabarnum, rabar, rzbar, ryelocal)
c      CALL YPRIME(RYPP,RTP,RYP)
C
C COPY OUTPUT WHICH IS STORED IN LOCAL ARRAY TO MATRIX OUTPUT
c      CALL MXCOPYREAL8TOPTR(RYPP, YPP, 4)
      CALL MXCOPYREAL8TOPTR(rabarnum, abarnump, 1)
      CALL MXCOPYREAL8TOPTR(rabar, abarp, 1)
      CALL MXCOPYREAL8TOPTR(rzbar, zbarp, 1)
      CALL MXCOPYREAL8TOPTR(ryelocal, yelocalp, 1)
C
      RETURN
      END
