RECURSIVE SUBROUTINE miser(func,regn,ndim,npts,dith,ave,var)
USE nrtype; USE nrutil, ONLY : assert_eq
IMPLICIT NONE
INTERFACE
FUNCTION func(x)
USE nrtype
IMPLICIT NONE
REAL(SP) :: func
REAL(SP), DIMENSION(:), INTENT(IN) :: x
END FUNCTION func
END INTERFACE
REAL(SP), DIMENSION(:), INTENT(IN) :: regn
INTEGER(I4B), INTENT(IN) :: ndim,npts
REAL(SP), INTENT(IN) :: dith
REAL(SP), INTENT(OUT) :: ave,var
REAL(SP), PARAMETER :: PFAC=0.1_sp,TINY=1.0e-30_sp,BIG=1.0e30_sp
INTEGER(I4B), PARAMETER :: MNPT=15,MNBS=60
Monte Carlo samples a user-supplied ndim -dimensional function func in a rectangular
volume specified by region , a 2× ndim vector consisting of ndim “lower-left” coordinates
of the region followed by ndim “upper-right” coordinates. The function is sampled a total
of npts times, at locations determined by the method of recursive stratified sampling. The
mean value of the function in the region is returned as ave ; an estimate of the statistical
uncertainty of ave (square of standard deviation) is returned as var . The input parameter
dith should normally be set to zero, but can be set to (e.g.) 0.1 if func ’s active region
falls on the boundary of a power-of-2 subdivision of region .
Parameters: PFAC is the fraction of remaining function evaluations used at each stage to
explore the variance of func . At least MNPT function evaluations are performed in any
terminal subregion; a subregion is further bisected only if at least MNBS function evaluations
are available.
REAL(SP), DIMENSION(:), ALLOCATABLE :: regn_temp
INTEGER(I4B) :: j,jb,n,ndum,npre,nptl,nptr
INTEGER(I4B), SAVE :: iran=0
REAL(SP) :: avel,varl,fracl,fval,rgl,rgm,rgr,&
s,sigl,siglb,sigr,sigrb,sm,sm2,sumb,sumr
REAL(SP), DIMENSION(:), ALLOCATABLE :: fmaxl,fmaxr,fminl,fminr,pt,rmid
ndum=assert_eq(size(regn),2*ndim,’miser’)
allocate(pt(ndim))
if (npts < MNBS) then
Too few points to bisect; do straight Monte
sm=0.0
Carlo.
sm2=0.0
do n=1,npts
call ranpt(pt,regn)
fval=func(pt)
sm=sm+fval
sm2=sm2+fval**2
end do
ave=sm/npts
var=max(TINY,(sm2-sm**2/npts)/npts**2)
else
Do the preliminary (uniform) sampling.
npre=max(int(npts*PFAC),MNPT)
allocate(rmid(ndim),fmaxl(ndim),fmaxr(ndim),fminl(ndim),fminr(ndim))
fminl(:)=BIG
Initialize the left and right bounds for each
fminr(:)=BIG
dimension.
fmaxl(:)=-BIG
fmaxr(:)=-BIG
do j=1,ndim
iran=mod(iran*2661+36979,175000)
s=sign(dith,real(iran-87500,sp))
rmid(j)=(0.5_sp+s)*regn(j)+(0.5_sp-s)*regn(ndim+j)
end do
do n=1,npre
Loop over the points in the sample.
call ranpt(pt,regn)
fval=func(pt)
where (pt <= rmid)
Find the left and right bounds for each di-
fminl=min(fminl,fval)
mension.
fmaxl=max(fmaxl,fval)
elsewhere
fminr=min(fminr,fval)
fmaxr=max(fmaxr,fval)
end where
end do
sumb=BIG
Choose which dimension jb to bisect.
jb=0
siglb=1.0
sigrb=1.0
do j=1,ndim
if (fmaxl(j) > fminl(j) .and. fmaxr(j) > fminr(j)) then
sigl=max(TINY,(fmaxl(j)-fminl(j))**(2.0_sp/3.0_sp))
sigr=max(TINY,(fmaxr(j)-fminr(j))**(2.0_sp/3.0_sp))
sumr=sigl+sigr
Equation (7.8.24); see text.
if (sumr <= sumb) then
sumb=sumr
jb=j
siglb=sigl
sigrb=sigr
end if
end if
end do
deallocate(fminr,fminl,fmaxr,fmaxl)
if (jb == 0) jb=1+(ndim*iran)/175000
MNPT may be too small.
rgl=regn(jb)
Apportion the remaining points between left
rgm=rmid(jb)
and right.
rgr=regn(ndim+jb)
fracl=abs((rgm-rgl)/(rgr-rgl))
nptl=(MNPT+(npts-npre-2*MNPT)*fracl*siglb/ &
Equation (7.8.23).
(fracl*siglb+(1.0_sp-fracl)*sigrb))
nptr=npts-npre-nptl
allocate(regn_temp(2*ndim))
regn_temp(:)=regn(:)
regn_temp(ndim+jb)=rmid(jb)
Set region to left.
call miser(func,regn_temp,ndim,nptl,dith,avel,varl)
Dispatch recursive call; will return back here eventually.
regn_temp(jb)=rmid(jb)
regn_temp(ndim+jb)=regn(ndim+jb)
Set region to right.
call miser(func,regn_temp,ndim,nptr,dith,ave,var)
Dispatch recursive call; will return back here eventually.
deallocate(regn_temp)
ave=fracl*avel+(1-fracl)*ave
Combine left and right regions by equation
var=fracl*fracl*varl+(1-fracl)*(1-fracl)*var
(7.8.11) (1st line).
	deallocate(rmid)
end if
deallocate(pt)

CONTAINS

SUBROUTINE ranpt(pt,region)
USE nr, ONLY : ran1
IMPLICIT NONE
REAL(SP), DIMENSION(:), INTENT(OUT) :: pt
REAL(SP), DIMENSION(:), INTENT(IN) :: region
	!Returns a uniformly random point pt in a rectangular region of dimension d. Used by
	!miser ; calls ran1 for uniform deviates.
INTEGER(I4B) :: n
call ran1(pt)
n=size(pt)
pt(1:n)=region(1:n)+(region(n+1:2*n)-region(1:n))*pt(1:n)
END SUBROUTINE ranpt
END SUBROUTINE miser
