      program Fract_Haz45
 
c     compatible with Haz45j
c     Last modified: 10/11/15

      implicit none

      include 'fract.h'
      
      integer nFlt, iflt
      integer nParamVar(MAX_WIDTH)
      integer  nWidth(MAX_FLT)
      character*80 filein, file1, fname, name1
      real cumWt_param(MAX_WIDTH,MAXPARAM),
     1     cumWt_Width(MAX_WIDTH)
      real mag(MAX_NMAG)
      real ran1, N(MAX_SAMPLE,MAX_NMAG)
      real sortarray(MAX_SAMPLE)
      integer i, nMag
      real tempx
      integer iSample, mcWidth, mcParam, iMag, iPerc, i1
      integer nperc
      real sum
      real perc(100,MAX_NMAG),mean(MAX_NMAG) 
      integer iseed, nSample, j, k
      real rate1(MAX_NMAG, MAX_WIDTH, MAXPARAM)
      real perc1(9)
      
      data perc1 / 0.01, 0.05, 0.1, 0.15, 0.50, 0.85, 0.90, 0.95, 0.99 /
    
      write (*,*) '*************************'
      write (*,*) '* SSC Mag Recur Fractile Code for use *'
      write (*,*) '*  with Hazard_v45j code  *'
      write (*,*) '*    Oct, 2015, NAA     *'
      write (*,*) '*************************'

c     Open run file
      write (*,*) 'Enter the input filename.'
      read (*,'(a80)') filein
      open (31,file=filein,status='old')

c     read run file parameters
      read (31,*) iseed
      read (31,*) nSample

c     Check array dimensions
      name1 = 'MAX_SAMPLE '
      call CheckDim ( nSample, MAX_SAMPLE, name1 )

c     Open output file
      read (31,'( a80)') file1
      open (30,file=file1,status='unknown')
      write (30,'( 2x,''iFlt  mag '',5x, ''mean'',8x,'' 1%'',9x,'' 5%'',9x,''10%'',
     1               9x,''15%'',9x,''50%'',9x,''85%'',9x,''90%'',9x,
     2               ''95%'',9x,''99%'',2x,''   name'')')   


c     Initialize random number generator
      do i=1,500
        tempx = Ran1 (iseed)
      enddo

c     Read header of the out2 file      
      write (*,'( 2x,''reading SSC logic tree file (out2)'')')
      call read_out2Head ( nFlt, nWidth, nParamVar,nMag, rate1,
     1         cumwt_param, cumwt_width, mag )
      
c     Loop over each source (this treats each segment separately) 
      do iFlt=1,nFlt

c       Read out2 file for this source      
        call read_out2 ( iFlt, nWidth, nParamVar,nMag, rate1,
     1         cumwt_param, cumwt_width, mag, fname )
        
c       Initialize mag Recur for this source
        do i=1,nSample
          do j=1,MAX_NMAG
            N(i,j) = 0.
          enddo
        enddo

c       Monte Carlo Sampling of Rates
        do iSample=1,nSample

c         Select the fault width for this source 
          call GetRandom0 ( iseed, nWidth(iFlt), cumWt_width, mcWidth )

c         Select the parameter variation for this source and fault width
          call GetRandom1b (iseed, nParamVar(mcWidth), 
     1          cumWt_param, mcWidth, mcParam, MAX_WIDTH, MAXPARAM)

c         Cumpute the cumulative rate (N)
c         (rate1 array is the incremental rate)
          do iMag=nMag,1,-1
            if ( iMag .eq. nMag) then
              N(iSample,iMag) = rate1(imag,mcWidth,mcParam)
            else 
              N(iSample,iMag) = N(iSample,iMag+1) + rate1(imag,mcWidth,mcParam)
            endif
          enddo

c       end loop over number of monte carlo samples
        enddo

c       Compute mean N
        do iMag = 1,nMag
          sum = 0.
          do iSample=1,nSample
            sum = sum + N(iSample,iMag)
          enddo
          mean(iMag) = sum/nSample
        enddo

c       Sort N data so that cummulative distribution function
c       can be computed  (sorts in place)
c       sortarray is a working array
        do iMag = 1,nMag
          call sort(N(1,iMag),sortarray,nSample)
        enddo

        nperc = 9
c       Compute fractile values
c        step = 1./(nperc+1)
        do iMag=1,nMag
          do iperc = 1,nPerc
            i1 = int( perc1(iPerc) * nSample )
  	    perc(iperc,iMag) = N(i1,iMag)
          enddo
	enddo

c       Write output
c        write (30,'( i5, 2x,a80)') iFlt, fname
        do iMag=1,nMag
          write (30,'(i5,2x,f4.2,10e12.4,2x,a80)') iFLt, mag(iMag), mean(iMag),
     1       (perc(i,iMag),i=1,nperc), fname
        enddo

      enddo
c     end loop on flts
      close (30)

      write (*,*) 
      write (*,*) '*** SSC Fractile Code (45) Completed with Normal Termination ***'

      stop
      end

c -----------------------------------------------------------

      subroutine read_out2Head ( nFlt, nWidth)

      implicit none
      include 'fract.h'

      integer nFlt, nWidth(MAX_FLT)
      integer nRD, iFlt, k
      character*80 fName, file1, dummy
      
      nRD = 12
                     
c     Open output file
      read (31,'( a80)') file1

      write (*,*) 'Opening output file from the hazard runs.'
      write (*,*) file1
      open (nRD,file=file1,status='old')
      read (nRD,*) nFlt
      read (nRD,*) (nWidth(k),k=1,nFlt)
      
c     Check dimensions
      do iFlt=1,nFlt
        if ( nWidth(iFlt) .gt. MAX_WIDTH) then
          write (*,'( 2x,''Increase MAX_WIDTH to'',i5)') nWidth(iFlt)
          stop 99
        endif
      enddo
      return


      end


c -----------------------------------------------------------

      subroutine read_out2 ( iFlt, nWidth, nParamVar, nMag, rate1,
     1         cumwt_param, cumwt_width, mag, fname )

      implicit none
      include 'fract.h'

      integer nFlt, nWidth(MAX_FLT), nMag, nParamVar(MAX_WIDTH)
      integer nRD, iMag, iFlt, iWidth, j, jFlt, jWidth, k
      character*80 fName, file1, dummy
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
