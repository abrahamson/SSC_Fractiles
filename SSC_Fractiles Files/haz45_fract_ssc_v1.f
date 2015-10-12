      program Fract_Haz45
 
c     compatible with Haz45j
c     Last modified: 10/11/15

      implicit none

      include 'fract.h'
      
      integer nFlt, iflt
      integer nParamVar(MAX_FLT,MAX_WIDTH)
      integer  nWidth(MAX_FLT)
      character*80 filein, file1, fname(MAX_FLT), name1
      real cumWt_param(MAX_FLT,MAX_WIDTH,MAXPARAM),
     1     cumWt_Width(MAX_FLT,MAX_WIDTH)
      real mag(MAX_FLT,MAXPARAM)
      real ran1, N(MAX_SAMPLE,MAX_NMAG)
      real sortarray(MAX_SAMPLE)
      integer i, nMag(MAX_FLT)
      real tempx
      integer iSample, mcWidth, mcParam, iMag, iPerc, i1
      integer nperc
      real sum

      real perc(100,MAX_NMAG),mean(MAX_NMAG) 
      integer iseed, nSample, j
      real rate1(MAX_FLT,MAX_NMAG,MAXPARAM)
      real perc1(9)
      data perc1 / 0.01, 0.05, 0.1, 0.15, 0.50, 0.85, 0.90, 0.95, 0.99 /
    
      write (*,*) '*************************'
      write (*,*) '* SSC Mag Recur Fractile Code for use *'
      write (*,*) '*  with Hazard_v45j code  *'
      write (*,*) '*    Oct, 2015, NAA     *'
      write (*,*) '*************************'

      write (*,*) 'Enter the input filename.'
      
      read (*,'(a80)') filein
      open (31,file=filein,status='old')

      read (31,*) iseed
      read (31,*) nSample

      name1 = 'MAX_SAMPLE '
      call CheckDim ( nSample, MAX_SAMPLE, name1 )

c     Open output file
      read (31,'( a80)') file1
      open (30,file=file1,status='new')
      write (30,'( i5, 2x,''number of sources'')') nFlt

c     Initialize random number generator
      do i=1,500
        tempx = Ran1 (iseed)
      enddo
       
      write (*,'( 2x,''reading SSC logic tree file (out2)'')')
      call read_out2 ( nFlt, nWidth, nParamVar,nMag, rate1,
     1         cumwt_param, cumwt_width, mag )
      write (*,'( 2x,''out of read out2'')')
      pause
      
c     Loop over each source (this treats each segment separately) 
      do iFlt=1,nFlt

c       Initialize mag Recur for this fault
        do i=1,nSample
          do j=1,MAX_NMAG
            N(i,j) = 0.
          enddo
        enddo

c       Monte Carlo Sampling of Rates
        do iSample=1,nSample

c         Select the fault width for this subsource 
          call GetRandom1 ( iseed, nWidth(iFlt), cumWt_width, iflt, mcWidth, MAX_FLT )

c         Select the parameter variation for this subsource and fault width
          call GetRandom2 (iseed, nParamVar(iFlt,mcWidth), 
     1          cumWt_param, iFlt, mcWidth, mcParam, MAX_FLT, MAX_WIDTH)

c         Cumpute the cumulative rate (N)
c         (the rate1 is the incremental rate)
          do iMag=nMag(iFlt),1,-1
            if ( iMag .eq. nMag(iFlt)) then
              N(iSample,iMag) = rate1(iFlt,mcWidth,mcParam)
            else 
              N(iSample,iMag) = N(iSample,iMag+1) + rate1(iFlt,mcWidth,mcParam)
            endif
          enddo

c       end loop over number of monte carlo samples
        enddo

c       Compute mean N
        do iMag = 1,nMag(iFlt)
          sum = 0.
          do iSample=1,nSample
            sum = sum + N(iSample,iMag)
          enddo
          mean(iMag) = sum/nSample
        enddo

c       Sort N data so that cummulative distribution function
c       can be computed  (sorts in place)
c       sortarray is a working array
        do iMag = 1,nMag(iFlt)
          call sort(N(1,iMag),sortarray,nSample)
        enddo

        nperc = 9
c       Compute fractile values
c        step = 1./(nperc+1)
        do iMag=1,nMag(iFlt)
          do iperc = 1,nPerc
            i1 = int( perc1(iPerc) * nSample )
  	    perc(iperc,iMag) = N(i1,iMag)
          enddo
	enddo

c       Write output
        write (30,'( i5, 2x,a80)') iFlt, fname(iFlt)
        write (30,'( 2x,''mag '', ''mean'',8x,'' 1%'',9x,'' 5%'',9x,''10%'',
     1               9x,''15%'',9x,''50%'',9x,''85%'',9x,''90%'',9x,
     2               ''95%'',9x,''99%'')')   
        do iMag=1,nMag(iFlt)
          write (30,'(2x,f4.2,10e12.4)') mag(iFlt,iMag), mean(iMag),
     1       (perc(i,iMag),i=1,nperc)
        enddo

      enddo
c     end loop on flts
      close (30)

      write (*,*) 
      write (*,*) '*** SSC Fractile Code (42) Completed with Normal Termination ***'

      stop
      end

c -----------------------------------------------------------

      subroutine read_out2 ( nFlt, nWidth, nParamVar,nMag, rate1,
     1         cumwt_param, cumwt_width, mag )

      implicit none
      include 'fract.h'

      integer nFlt, nWidth(MAX_FLT), nMag(MAX_FLT), 
     1        nParamVar(MAX_FLT,MAX_WIDTH)
      integer nRD, iMag, iFlt, iFltWidth, j, jFlt, jFltWidth, k
      character*80 fName, file1, dummy
      real rate1(MAX_FLT, MAX_NMAG, MAX_WIDTH, MAXPARAM)
      real wt(MAX_FLT,MAX_WIDTH,MAXPARAM)
      real cumwt_param(MAX_FLT,MAX_WIDTH,MAXPARAM)
      real wt_width(MAX_WIDTH)
      real cumwt_width(MAX_FLT,MAX_WIDTH)
      real mag(MAX_FLT,MAXPARAM), temp(MAXPARAM)
      
      nRD = 12
                     
c     Open output file
      read (31,'( a80)') file1

      write (*,*) 'Opening output file from the hazard runs.'
      write (*,*) file1
      open (nRD,file=file1,status='old')
      read (nRD,*) nFlt
      read (nRD,*) (nWidth(k),k=1,nFlt)
 
      do iFlt=1,nFlt
        read (nRd,*) (wt_width(k),k=1,nWidth(iFlt))
        do iFltWidth=1,nWidth(iFlt) 
          read (nRD,*) nMag(iFlt), nParamVar(iFlt,iFltWidth)
          read (nRD,*) (wt(iFlt,iFltWidth,j),j=1,nParamVar(iFlt,iFltWidth)) 
          do iMag=1,nMag(iFlt)
            read (nRD,*) jFlt, jFltWidth, mag(iFlt,imag),
     1       (temp(j),j=1,nParamVar(iFlt,iFltWidth))
            do j=1,nParamVar(iFlt,iFltWidth)
              rate1(iFlt,iMag,iFltWidth,j) = temp(j)
            enddo
          enddo
        enddo
      enddo
      
c     Compute cumulative wts
      do iFlt=1,nFlt
        do iFltWidth=1,nWidth(iFlt)
          cumwt_param(iFlt,iFltWidth,1) = wt(iFlt,iFltWidth,1)
          do j=2,nParamVar(iFlt,iFltWidth)
            cumwt_param(iFlt,iFltWidth,j) = wt(iFlt,iFltWidth,j)+
     1             cumwt_param(iFlt,iFltWidth,j-1)
          enddo
        enddo
      enddo      

      close (nRD)
c      enddo

  100 continue

      return
  
      write (*,'( 2x,''bad site number'')')
      stop 98

      end

c -------------------------------------------------------------------------
