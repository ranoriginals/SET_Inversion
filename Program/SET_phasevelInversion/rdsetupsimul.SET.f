*  modified rdsetupsimulf from DWF.
 
c  rdsetupsimulf.1.f
c  This program is the part of simannerr that reads in the data and does the
c  Fourier analysis.  It is separated so that simannerr# can run on katmai or
c  other machines that do not have sac.  Writes out the data ready for use in
c  simannerr#  
c  Pipe in data, e.g.,   rdsetupsimul < eqlistper50
c
c  Current version does not make station corrections.
c  "ranor station number > 300 , 300 -> 500"
      parameter (maxnfreq=20, maxnsta=500, 
     1     maxpts = 60000, maxnodes = 1000, maxevnts = 500)

      real*4 staph(maxevnts,maxnsta,maxnfreq)
      real*4 staamp(maxevnts,maxnsta,maxnfreq)
      real*4 freq(maxnfreq)
      real*4 stadist(maxevnts,maxnsta), staazi(maxevnts,maxnsta)
      real*4 stacor(maxnsta), tempcor(maxnsta),geomsprd(maxnsta)
      real*4 bazi(maxevnts,maxnsta), stadelt(maxevnts,maxnsta)
      real*4 tdata(maxpts),beg(maxevnts,maxnsta),delt(maxnsta)
      real*4 stalat(maxevnts,maxnsta),stalon(maxevnts,maxnsta)
      real*4 attnfac(maxnfreq),fattnfac(10,2)
      integer*4 nsta(maxevnts), dot, idot
      integer*4 nfreq,nstapts(maxnsta),idnum(maxnsta)
      integer*4 istanum(maxnsta),istacor(maxnsta), nevents
      character*70 foutput, fsummary, fstalist

c   modify here       
      character*125 fn(maxevnts,maxnsta)
      character*70 finvrsnodes, fftoutput, fvariance,fmaxavamp
      character*70 ftemp,fvelarea
      character*16 staampcor
c      非常重要！！！
c      character*69 dummy1   ! cutting length depends on the path
      character*2 nettemp
c     ranor : character*4 staname(maxnsta)表示定义一个一维字符数组，有maxnsta个字符串，每个可容纳4个字符串
c     由于SET的台站名都是5个，应该改成5
c      character*4 statemp, staname(maxnsta)
      character*5 statemp, staname(maxnsta)
      character*8 statemp2
      

      logical debug
      debug = .true.

      pi = 3.1415928
      convdeg = 3.1415928/180.
      circ = 6371.*3.1415928/180.
      twopi = 3.1415928*2.

c  **************************************************
c  WARNING:  The following two statements need to be switched depending on which
c  compiler is used
c  ************************************************

c200         format(a75, i2)
c201         format(a70, a2)
c202         format(a69, a4)
c210	       format(a49, i2)
c200         format(1H ,a44, i2)
c200         format(a44, i2)

c  read in frequency and its corresponding attenuation coefficient ( data from 
c  Brian Mitchell, 1995)
c      data(fattnfac(i,1), fattnfac(i,2), i=1,6)/0.05, 0.15e-3, 0.025, 
c     * 0.15e-3, 0.01667, 0.15e-3, 0.01429, 0.1e-3, 0.0125, 0.05e-3, 
c     * 0.0111, 0.05e-3/  

      cattnfac=0.25e-3

c  read list of files and frequencies to be analyzed and file to output results
c  Usually will pipe in data from some file like fasearrayinp
      read(*,*) nevents
c     读入对应频带的总事件数
c     ranor:debug
c      write(*,*) nevents
c      stop
c      nobs = 0

      do iev = 1, nevents
        read(*,*) nsta(iev), idnum(iev)
c        读入每个事件有多少个台以及每个事件的idnum
c        ranor:debug
c        write(*,*) nsta(iev),idnum(iev)
c        nobs = nobs + 2*nsta(iev)
        read(*,'(a)') (fn(iev,i), i=1,nsta(iev))
c        读入每个事件每个台站记录的路径
c        ranor:debug
c        write(*,'(a)') (fn(iev,1),i=1,nsta(iev))
c       if(debug) write(*,*) iev, nsta(iev), idnum(iev)
      enddo
c      stop
c  a lot of the following input is unneeded, just here to parallel simannerr#
c  input, so can use same fasearrayinp# or eqlistper# file for both
      read(*,*) nfreq
c      write(*,*) nfreq
c     对应第一行，表示这次计算有几个主频的数据
      read(*,*) (freq(j), j= 1, nfreq)
c      write(*,*) (freq(j), j= 1, nfreq)
c      对应第二行，表示这次计算的主频
      read(*,'(a)') foutput
c      write(*,'(a)') foutput
c      对应第三行：detail.040.591.080.JDF1
      read(*,'(a)') fsummary
c      write(*,'(a)') fsummary
c      对应第四行：summar.040.591.080.JDF1
      read(*,'(a)') finvrsnodes
c      write(*,'(a)') finvrsnodes
c      对应第五行：JDFgridinp.dat 研究区域的网格化文件
      read(*,'(a)') fftoutput
c      write(*,'(a)') fftoutput
c      对应第六行：phampcor.040
      read(*,'(a)') fvariance
c      write(*,'(a)') fvariance
c      对应第七行：covar.040.591.080.JDF1
      read(*,'(a)') fmaxavamp
c      write(*,'(a)') fmaxavamp
c      对应第八行：mavamp.040.591.080.JDF1
      read(*,'(a)') ftemp
c      write(*,'(a)') ftemp
c      对应第九行：tempd
      read(*,'(a)') fvelarea
c      write(*,'(a)') fvelarea
c      对应第十行：JDFrgnlBndryPts
      read(*,'(a)') fstalist
c      write(*,'(a)') fstalist
c      对应第十一行：stanumlist 台站列表文件
      read(*,*) unifvel
c      write(*,*) unifvel
c      对应第十二行：该主频对应的相速度平均值
      read(*,*) iterlimit, dampvel, dampaniso
c      write(*,*) iterlimit, dampvel, dampaniso
c      对应第十三行：迭代次数，阻尼值
c      stop
c      ranor: for debug      
c      write(*,*) foutput,fsummary,finvrsnodes,fstalist
c      stop

c      idot=dot(foutput)
c      write(*,*) idot
c      staampcor = 'staampcor'//foutput(idot+1:idot+4)
c      write(14,*) staampcor
c      staampcor = 'staampcorBR3.dat'
c     file11 - phampcor.040
c     file13 - tempd
c     file14 - followit13
c     file18-stanumlist
      open(11, file = fftoutput)
c      open(12, file = staampcor, status='old')

      open(13, file = ftemp)
      open(14, file = 'followit13')
      open(18, file = fstalist)

c  fetch station corrections for amplitudes and assign to correction array
c      read(12,*) nstacor
c      do i = 1,nstacor
c        read(12,*) istacor(i),tempcor(i)
c      enddo
c      write(*,*) nevents
c      write(13,*) nevents
      
c      do i = 1, nstacor
c        stacor(istacor(i)) = tempcor(i)
c      write(14,*)  istacor(i), stacor(istacor(i))
c      enddo
c

c  find the attenuation factor for a given frequency
c      do i=1, nfreq
c       do j=1, 6
c        if ((freq(i).le.fattnfac(j,1)).and.(freq(i).ge.
c     1  fattnfac(j+1,1))) then
c         attnfac(i) = fattnfac(j+1,2)+(freq(i)-fattnfac(j+1,1))* 
c     1   (fattnfac(j,2)-fattnfac(j+1,2))/(fattnfac(j,1)-
c     1   fattnfac(j+1,1))
c	end if
c       enddo		!j

c       if (freq(i).gt.fattnfac(1,1)) then
c	  attnfac(i) = fattnfac(1,2)
c       end if

c       if (freq(i).lt.fattnfac(6,1)) then
c	  attnfac(i) = fattnfac(6,2)
c       end if

c	write(*,*) "attnfac:  ", attnfac(i)
c      enddo		!i 	

c  generate attenuation factor assuming Q = 100 and vel = 3.80km/s
c       do i = 1, nfreq
c         attnfac(i) = pi*freq(i)/380.
c       enddo

c  new attenuation factor based on 25 s inversion with multipathing
c  taken into account is much smaller with Q about 385.  Average 
c  velocity more like 3.80

c	do i = 1, nfreq
c           attnfac(i) = pi*freq(i)/(3.80*385.)
c       enddo

c  first, read in master list of stations  
c  这一段读取stanumlist中的台站名  
      do ista2 = 1, maxnsta
c       ranor: debug
c        read(*,*) staname(ista2)
        read(18,'(5a)') staname(ista2)
c       ranor:debug
c        write(*,'(5a)') staname(ista2)
c        ranor: debug
c        write(*,*) staname(ista2)
c        if(debug) write(*,*) ista2, staname(ista2)
        if (staname(ista2).eq.'nope') then
c         ranor: error maxnsta -> ista2
          mxnsta = ista2 - 1
c          mxnsta = maxnsta - 1
          continue 
          go to 1111
        endif
      enddo
1111  continue
c      stop

c      ranor: debug      
c      write(*,*) mxnsta,ista2
c      stop

c  start input loop over events
      do iev = 1, nevents
c        每个事件写入followit13，分别是idnum和一个事件所含的台站记录数
        write(14,*) idnum(iev), nsta(iev)
c       ranor:debug
c       write(*,*) idnum(iev), nsta(iev)

        do ista = 1, nsta(iev)
c          关于rsac1 & getfhv & getkhv：http://blog.sina.com.cn/s/blog_872ecf9f0100vckx.html
          call rsac1(fn(iev,ista),tdata,nstapts(ista),beg(iev,ista),
     1       delt(ista),maxpts,nerr)
c         ranor: debug
c          write(*,*) fn(iev,ista),nstapts(ista),beg(iev,ista),delt(ista),maxpts
c          stop
c         rsac1---fn(iev,ista): 每个台站记录所在的路径
c                 tdata: 数组用于储存读入的文件中每个数据点的振幅值
c                 nstapts: 数组的长度  10001(2000*5+1)
c                 beg: 数据开始时间 每个事件的起始点
c                 delt: 数据采样间隔 0.2
c                 maxpts: 所读入的最大的数据点数 60000
c                 nerr: 错误返回标记，0代表成功，非零代表失败
c         rsac1: 读取等间隔文件
c         getfhv: 获取浮点型头段变量的值，nerr: error return flag
c                 三个参数: 1.字符串，代表要读取的头段变量名 2. 表示将接收的头段变量存在什么变量内 3.错误返回标记
c         getkhv: 获取字符串头段变量值
c         call后所有的write为了检查输出添加--输出正确
          call getfhv('DIST',stadist(iev,ista),nerr)
c          write(*,*) stadist(iev,ista)
          call getfhv('AZ', staazi(iev,ista),nerr)
c          write(*,*) staazi(iev,ista)
          call getfhv('BAZ', bazi(iev,ista),nerr)
c          write(*,*) bazi(iev,ista)
          call getfhv('GCARC', stadelt(iev,ista), nerr)
c          write(*,*) stadelt(iev,ista)
          call getfhv('STLA', stalat(iev,ista),nerr)
c          write(*,*) stalat(iev,ista)
          call getfhv('STLO', stalon(iev,ista),nerr)
c          write(*,*) stalon(iev,ista)
          call getkhv('KSTNM', statemp2,nerr)
c          write(*,*) statemp2
c          stop

c         第一个事件第一个台站的baz
          rembazi = bazi(iev,1)
c          write(*,*) rembazi
c          第一个事件第一个台站的dist
          remdist = stadist(iev,1)
c          write(*,*) stadist(iev,1)
c         stop

c         .gt.在Fortran中表示大于
c         ranor:debug
c          write(*,*) rembazi-bazi(iev,ista)

          if ((abs(rembazi-   bazi(iev,ista)).gt.4.0).or.
     1        (abs(remdist-stadist(iev,ista)).gt.600.)) then 
            write(14,*) bazi(iev,ista),stadist(iev,ista),rembazi,remdist
            write(14,*) fn(iev,ista)
c           ranor:debug 
c            write(*,*) bazi(iev,ista),stadist(iev,ista),rembazi,remdist
c            write(*,*) fn(iev,ista)
          endif
c          stop
          geomsprd(ista) = sqrt(sin(stadelt(iev,ista)*convdeg))

          statemp = trim(statemp2)
          write(14,*) statemp2, statemp

          istanum(ista) = 0
          
          do ista2 = 1, mxnsta
             if (statemp.eq.staname(ista2)) then
               istanum(ista) = ista2
               write(*,2070)stalat(iev,ista),stalon(iev,ista),statemp,ista2
             endif
          enddo
      
c      stop
c 2070   format(f10.4, f10.4,2x, a4, 2x, i3)
c ranor: 台站名4-5
2070   format(f10.4, f10.4,2x, a5, 2x, i3)

          if (istanum(ista).eq.0) then
            write(*,*) statemp, ista
            write(*,*) 'WARNING ', statemp,' not in station list'
          endif
  
c ********************************modify here***********************************

c          read (13,201) dummy1, nettemp
c	  if (index(nettemp,'NR')) then
c	     rewind(13)
c	     read(13,200) dummy1, istanum(ista)
c	     istanum(ista)=istanum(ista)-69
c	  else 
c	     rewind(13)
c	     read(13,202) dummy1,statemp
c	    
c	     if (index(statemp,'DGR')) istanum(ista)=15
c	     if (index(statemp,'PLM')) istanum(ista)=16
c	     if (index(statemp,'JCS')) istanum(ista)=17
c	     if (index(statemp,'BAR')) istanum(ista)=18
c	     if (index(statemp,'BC3')) istanum(ista)=19
c	     if (index(statemp,'SWS')) istanum(ista)=20
c	     if (index(statemp,'GLA')) istanum(ista)=21
c	     if (index(statemp,'SNCC')) istanum(ista)=22
c	     if (index(statemp,'CHXB')) istanum(ista)=23
c	     if (index(statemp,'PPXB')) istanum(ista)=24
c	     if (index(statemp,'BAHB')) istanum(ista)=25
c	     if (index(statemp,'TOPB')) istanum(ista)=26
c	     if (index(statemp,'GUAY')) istanum(ista)=27
c	  end if

c	  write(14,*) dummy1, istanum(ista)
c          write(14,*) fn(iev,ista)
c          rewind (13)
c          write(14,*) nstapts(ista),delt(ista)
c
c  get phases and amplitudes at desired frequencies
c
      
         
          do ifreq = 1, nfreq
            call frt(tdata, freq(ifreq), staamp(iev,ista,ifreq),
     1               staph(iev,ista,ifreq), nstapts(ista),delt(ista))

c  correct amplitudes for geometrical spreading, attenuation and station factor
            attneffect= exp(cattnfac*(stadist(iev,ista)-
     1                                         stadist(iev,1)))
c            attneffect= exp(attnfac(ifreq)*(stadist(iev,ista)-
c     1                                         stadist(iev,1)))

c            staamp(iev,ista,ifreq)= staamp(iev,ista,ifreq)*attneffect
c     1         *(geomsprd(ista)/stacor(istanum(ista)))
            staamp(iev,ista,ifreq)= staamp(iev,ista,ifreq)*attneffect
     1         *(geomsprd(ista))

          enddo  
        enddo
      enddo
c      stop
c  end of input loop over events
c  output loop over events
      do iev = 1, nevents
c       11-ffoutput
        write (11,*) iev, idnum(iev)
c       ranor :debug 
c        write (*,*) iev, idnum(iev)
        do ista = 1, nsta(iev)
          write(11,*) beg(iev,ista)
c         ranor:debug
c          write(*,*) beg(iev,ista)
c         stadist -- dist
c         staazi  -- az
c         bazi    -- baz
c         stadelt -- gcarc
c         stalat  -- stla
c         stalon  -- stlo
c         查看phampcor.040发现第一个write写进去文件出现了换行
          write (11,111) stadist(iev,ista), staazi(iev,ista), bazi(iev,ista),
     1                 stadelt(iev,ista), stalat(iev,ista), stalon(iev,ista)
c          write (*,111) stadist(iev,ista)
          write(11,*) (staamp(iev,ista,ifreq),staph(iev,ista,ifreq),
     1                 ifreq=1,nfreq)
111       format(f10.4,f10.4,f10.4,f10.4,f10.4,f10.4)
c         ranor:debug     
c          write (*,*) stadist(iev,ista), staazi(iev,ista), bazi(iev,ista),
c     1                 stadelt(iev,ista), stalat(iev,ista), stalon(iev,ista)
c          write(*,*) (staamp(iev,ista,ifreq),staph(iev,ista,ifreq),
c     1                 ifreq=1,nfreq)
        enddo
      enddo
c     file11 - phampcor.040
c     file13 - tempd
c     file14 - followit13
c     file18 - stanumlist
      close(unit = 11)
c      close(unit = 12)
c      status表示这个文件要在close后被删除
      close(unit = 13,status='delete')
      close(unit = 14)
      close(unit = 18)
      end
      
c---positive fourier transform
c 
      SUBROUTINE FRT(UP,FR,ARZ,PRZ,NALL,DELT)
      dimension UP(40000)
      DIMENSION W(40000)
      THETA=6.283185*FR*DELT
      C=COS(THETA)*2.0
      NR1=1
      NR2=NALL
      NDR1=NR2-1
      W(1)=UP(NR2)
      W(2)=C*W(1)+UP(NDR1)
      NDATA=NR2-NR1+1
      NTR1=NDATA-1
      NTR2=NDATA-2
      DO 97 I=3,NDATA
      I1=I-1
      I2=I-2
      NDRI=NR2-I+1
      W(I)=C*W(I1)-W(I2)+UP(NDRI)
97    CONTINUE
      ZRE=(W(NDATA)-W(NTR2)+UP(NR1))*DELT/2.0
      ZIM=W(NTR1)*SIN(THETA)*DELT
      CALL PHASE(ZRE,ZIM,PRZ)
c      ARZ=SQRT(ZRE*ZRE+ZIM*ZIM)
      ARZ=SQRT(ZRE*ZRE+ZIM*ZIM)*1.0e-9
      RETURN
      END

      SUBROUTINE PHASE(X,Y,PHI)
      IF(X) 21,20,22
20    IF(Y) 23,24,25
23    PHI=1.5*3.141592
      GO TO 28
24    PHI=0.0
      GO TO 28
25    PHI=0.5*3.141592
      GO TO 28
21    PHI=ATAN(Y/X) +3.141592
      GO TO 28
22    IF(Y) 26,27,27
27    PHI=ATAN(Y/X)
      GO TO 28
26    PHI=ATAN(Y/X)+2.0*3.141592
      GO TO 28
28    CONTINUE
      PHI=PHI/6.283184
      PHI=PHI-AINT(PHI)
      RETURN
      END

      integer function dot(file)
      character file*70
      do 50 i=1,70
      if(file(i:i).ne.'.') goto 50
      dot=i-1
      return
50     continue
      write(1,100) file
100   format(' no dots found in ',a70)
      dot = 0
      return
      end
