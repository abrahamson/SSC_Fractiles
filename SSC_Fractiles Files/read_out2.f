
      subroutine read_out2Head ( nFlt, nWidth)

      implicit none
      include 'fract.h'

      integer nFlt, nWidth(MAX_FLT)
      integer nRD, iFlt, k
      real version
      character*80 file1
      
      nRD = 12
                     
c     Open output file
      read (31,'( a80)',ERR=2000) file1

      write (*,*) 'Opening output file from the hazard runs.'
      write (*,*) file1
      open (nRD,file=file1,status='old')
      
C     Check for version compatibility with hazard code
        read (nRD,*,ERR=2001) version
         if (version .ne. 45.2) then
           write (*,*) 'out2 from incompatible version of Haz45, use Haz45.2'
           stop 99
         endif      
      
      read (nRD,*,ERR=2002) nFlt
      read (nRD,*,ERR=2003) (nWidth(k),k=1,nFlt)
      
c     Check dimensions
      do iFlt=1,nFlt
        if ( nWidth(iFlt) .gt. MAX_WIDTH) then
          write (*,'( 2x,''Increase MAX_WIDTH to'',i5)') nWidth(iFlt)
          stop 99
        endif
      enddo
      return
      
 2000 write (*,'( 2x,''error reading Hazard output file name'' )') 
      write (*,'( 2x,''file: '',a80)') file1
      stop
 2001 write (*,'( 2x,''error reading version'' )') 
      stop
 2002 write (*,'( 2x,''error reading nFlt'' )') 
      stop
 2003 write (*,'( 2x,''error reading nWidth'' )') 
      stop
      end


c -----------------------------------------------------------

      subroutine read_out2 ( iFlt, nWidth, nParamVar, nMag, rate1,
     1         cumwt_param, cumwt_width, mag, fname )

      implicit none
      include 'fract.h'

      integer nWidth(MAX_FLT), nMag, nParamVar(MAX_WIDTH)
      integer nRD, iMag, iFlt, iWidth, j, jFlt, jWidth, k
      character*80 fName
      real rate1(MAX_NMAG, MAX_WIDTH, MAXPARAM)
      
      real wt(MAX_WIDTH,MAXPARAM)
      real cumwt_param(MAX_WIDTH,MAXPARAM)
      real wt_width(MAX_WIDTH)
      real cumwt_width(MAX_WIDTH)
      real mag(MAX_NMAG), temp(MAXPARAM)
      real sum
      
      nRD = 12
      
      read (nRd,*) (wt_width(k),k=1,nWidth(iFlt))
      
      do iWidth=1,nWidth(iFlt)
        read (nRd,'( i5,2x,a80)') jFlt, fname        
        read (nRD,*) nMag, nParamVar(iWidth)

        if ( nParamVar(iWidth) .gt. MAXPARAM ) then
          write (*,'( 2x,''increase MAXPARAM to'',i5)') nParamVar(iWidth)
          stop 99
        endif
        if ( nMag .gt. MAX_NMAG ) then
          write (*,'( 2x,''increase MAX_NMAG to'',i5)') nMag 
          stop 99
        endif
        read (nRD,*) (wt(iWidth,j),j=1,nParamVar(iWidth)) 

        do iMag=1,nMag
          read (nRD,*) jFlt, jWidth, mag(imag),(temp(j),j=1,nParamVar(iWidth))
          do j=1,nParamVar(iWidth)
            rate1(iMag,iWidth,j) = temp(j)
            if ( temp(j) .lt. 0. ) then
              write (*,'( 4i5,e12.4)') iFlt, iWidth, iMag,j, temp(j)
c              pause ‘bad rate (negative)’
            endif
          enddo
        enddo
      enddo
      
c     Compute cumulative wts for the parameter variations
      do iWidth=1,nWidth(iFlt)

c       First, check that weights sum to unity
        sum = 0.
        do j=1,nParamVar(iWidth)
          sum = sum + wt(iWidth,j)
        enddo
        if ( sum .lt. 0.99 .or. sum .gt. 1.01 ) then
          write (*,'( 2x,''Parameter weights do not sum to 1'')') 
          write (*,'( 2x,''iflt='',i4,2x,''iWidth='',i4,2x,''sum wt='',f10.4)') 
     1         iFlt, iWidth, sum 
          stop 99
        endif

c       Set cumulative weights , forcing sum to unity       
        cumwt_param(iWidth,1) = wt(iWidth,1)
        do j=2,nParamVar(iWidth)
          cumwt_param(iWidth,j) = wt(iWidth,j)+ cumwt_param(iWidth,j-1)
        enddo
        cumwt_param(iWidth,nParamVar(iWidth)) = 1.0
      enddo

c     Check if width weights sum to unity
      sum = 0.
      do iWidth=1,nWidth(iFlt)
        sum = sum + wt_width(iWidth)
      enddo     
      if ( sum .lt. 0.99 .or. sum .gt. 1.01 ) then
        write (*,'( 2x,''Width weights do not sum to 1'')') 
        write (*,'( 2x,''iflt='',i4,2x,''sum wt='',f10.4)') iFlt, sum 
        stop 99     
      endif

c     Compute cumulative wts for the width variations
c     forcing sum to unity       
      cumwt_width(1) = wt_width(1)
      do iWidth=2,nWidth(iFlt)
        cumwt_width(iWidth) = wt_width(iWidth)+ cumwt_width(iWidth-1)
      enddo
      cumwt_width(nWidth(iFlt)) = 1.0

      return
  
      end
