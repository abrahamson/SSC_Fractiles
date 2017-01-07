      program Fract_Haz45
 
c     compatible with Haz45.2

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
      integer iseed, nSample, j
      real rate1(MAX_NMAG, MAX_WIDTH, MAXPARAM)
      real perc1(9)
      
      data perc1 / 0.01, 0.05, 0.1, 0.15, 0.50, 0.85, 0.90, 0.95, 0.99 /
    
      write (*,*) '*******************************'
      write (*,*) '* SSC Mag Recur Fractile Code *'
      write (*,*) '*  compatible with Haz45.2    *'
      write (*,*) '*      Sept, 2016, NAA         *'
      write (*,*) '*******************************'

c     Open run file
      write (*,*) 'Enter the input filename.'
      read (*,'(a80)') filein
      open (31,file=filein,status='old',ERR=2000)

c     read run file parameters
      read (31,*,ERR=2001) iseed
      read (31,*,ERR=2002) nSample

c     Check array dimensions
      name1 = 'MAX_SAMPLE '
      call CheckDim ( nSample, MAX_SAMPLE, name1 )

c     Open output file
      read (31,'( a80)',ERR=2003) file1
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
            if ( N(iSample,iMag) .lt. 0. ) write (*,'( 2i5,e12.4)') iMag, iSample, n(iSample,iMag)
          enddo
          mean(iMag) = sum/nSample
          if ( iFlt .eq. 119 ) then
            write (*,'( i5,2e12.5)') iMag, sum, mean(iMag)
          endif
        enddo
        if ( iFlt .eq. 119 ) pause

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
      write (*,*) '*** SSC Fractile Code Completed with Normal Termination ***'

      stop
 2000 write (*,'( 2x,''error reading input file'' )') 
      write (*,'( 2x,''file: '',a80)') filein
      stop
 2001 write (*,'( 2x,''error reading iSeed'' )') 
      stop
 2002 write (*,'( 2x,''error reading nSample'' )') 
      stop
 2003 write (*,'( 2x,''error reading output file name'' )') 
      stop

      end      
      