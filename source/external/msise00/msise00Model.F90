MODULE NRLMSISE00_FUNCTIONS

    
        USE MSISE_CONSTANTS
        USE NRLMSISE_TYPE
        USE groops_msise00_kinds
        IMPLICIT NONE
        !PUBLIC :: msise00Density

    CONTAINS

! ======================================================================
! C++: int glatf(const double& lat, double& gv, double& reff)
! F90: SUBROUTINE glatf(lat, gv, reff)
! 
! C++ 版本返回一个 int 0，但主要工作是通过引用修改 gv 和 reff。
! Fortran 的 SUBROUTINE (子程序) 是这种情况的
! 最佳对等体，使用 INTENT(OUT) 声明输出参数。
! ======================================================================    
    SUBROUTINE glatf(lat, gv, reff)
        ! 输入参数
        REAL(RL), INTENT(IN)  :: lat
        ! 输出参数
        REAL(RL), INTENT(OUT) :: gv, reff
        
        ! 局部变量
        REAL(RL) :: c2
        
        c2 = COS(2.0 * d_DGTR * lat)
        gv = 980.616 * (1.0 - 0.0026373 * c2)
        reff = 2.0 * gv / (3.085462e-6 + 2.27e-9 * c2) * 1e-5
        
    END SUBROUTINE glatf

! ======================================================================
! C++: double ccor(const double& alt, const double& r, const double h1, const double zh)
! ======================================================================    
    
    FUNCTION ccor(alt, r, h1, zh) RESULT(ccor_res)
        ! 输入参数
        REAL(RL), INTENT(IN) :: alt, r, h1, zh
        ! 返回值
        REAL(RL) :: ccor_res
        
        ! 局部变量
        REAL(RL) :: e
        
        e = (alt - zh) / h1
        
        IF (e > 70.0) THEN
            ccor_res = 1.0
            RETURN
        END IF
        
        IF (e < -70.0) THEN
            ccor_res = EXP(r)
            RETURN
        END IF
        
        ccor_res = EXP(r / (1.0 + EXP(e)))
        
    END FUNCTION ccor

    
! ======================================================================
! C++: double ccor2(const double& alt, const double& r, const double h1, const double zh, const double h2)
! ======================================================================    
    FUNCTION ccor2(alt, r, h1, zh, h2) RESULT(ccor2_res)
        ! 输入参数
        REAL(RL), INTENT(IN) :: alt, r, h1, zh, h2
        ! 返回值
        REAL(RL) :: ccor2_res
        
        ! 局部变量
        REAL(RL) :: e1, e2
        
        e1 = (alt - zh) / h1
        e2 = (alt - zh) / h2
        
        ! C++: || (or) -> .OR.
        IF ((e1 > 70.0) .OR. (e2 > 70.0)) THEN
            ccor2_res = 1.0  ! C++: exp(0)
            RETURN
        END IF
        
        ! C++: && (and) -> .AND.
        IF ((e1 < -70.0) .AND. (e2 < -70.0)) THEN
            ccor2_res = EXP(r)
            RETURN
        END IF
        
        ccor2_res = EXP(r / (1.0 + 0.5 * (EXP(e1) + EXP(e2))))
        
    END FUNCTION ccor2
    
! ======================================================================
! C++: double scalh(NRLMSISE NRL,const double& alt, const double xm, const double temp)
!
! C++: pow(base, exp) -> F90: base**exp
! C++: NRL.d_gsurf    -> F90: NRL%d_gsurf
! ======================================================================
    FUNCTION scalh(NRL, alt, xm, temp) RESULT(scalh_res)
        ! 输入参数
        TYPE(NRLMSISE), INTENT(IN) :: NRL  ! 派生类型
        REAL(RL), INTENT(IN) :: alt, xm, temp
        ! 返回值
        REAL(RL) :: scalh_res
        
        ! 局部变量
        REAL(RL) :: g
        
        ! 使用 Fortran 的 ** 运算符进行幂运算
        ! 使用 Fortran 的 % 访问派生类型成员
        g = NRL%d_gsurf / ((1.0 + alt / d_re)**2.0)
        
        scalh_res = d_RGAS * temp / (g * xm)
        
    END FUNCTION scalh   
    
! ======================================================================
! C++: double dnet(double& dd, const double& dm, const double& zhm, const double& xmm, const double xm)
!
! 'dd' 是 C++ 中的引用 (double&)，在 F90 中是 INTENT(INOUT)
! C++: printf(...) -> F90: PRINT * 或 WRITE(*,*)
! C++: log() (自然对数) -> F90: LOG() (自然对数)
! C++: std::exp() -> F90: EXP()
! ======================================================================
    FUNCTION dnet(dd, dm, zhm, xmm, xm) RESULT(dnet_res)
        ! TURBOPAUSE CORRECTION FOR MSIS MODELS
        ! DD - diffusive density
        ! DM - full mixed density
        ! ZHM - transition scale length
        ! XMM - full mixed molecular weight
        ! XM - species molecular weight
        ! DNET - combined density
        
        ! 输入/输出参数
        REAL(RL), INTENT(INOUT) :: dd
        ! 输入参数
        REAL(RL), INTENT(IN)    :: dm, zhm, xmm, xm
        ! 返回值
        REAL(RL) :: dnet_res
        
        ! 局部变量
        REAL(RL) :: a, ylog
        
        a = zhm / (xmm - xm)
        
        ! C++: !((dm > 0) && (dd > 0)) -> (dm <= 0) || (dd <= 0)
        IF ((dm <= 0.0) .OR. (dd <= 0.0)) THEN
            ! C++: printf(...)
            PRINT *, 'dnet log error ', dm, dd, xm
            
            IF ((dd == 0.0) .AND. (dm == 0.0)) THEN
                dd = 1.0
            END IF
            
            IF (dm == 0.0) THEN
                dnet_res = dd
                RETURN
            END IF
            
            IF (dd == 0.0) THEN
                dnet_res = dm
                RETURN
            END IF
        END IF
        
        ylog = a * LOG(dm / dd)
        
        IF (ylog < -10.0) THEN
            dnet_res = dd
            RETURN
        END IF
        
        IF (ylog > 10.0) THEN
            dnet_res = dm
            RETURN
        END IF
        
        ! C++: pow(base, 1.0/a) -> base**(1.0/a)
        dnet_res = dd * (1.0 + EXP(ylog))**(1.0 / a)
        
    END FUNCTION dnet
    
    ! ======================================================================
    ! C++: double splini(double *xa, double *ya, double *y2a, int n, double x)
    ! 
    ! 警告：此函数从 C++ 逐字翻译。
    ! C++ 版本似乎混合了 0-based 和 1-based 索引逻辑，
    ! 并且在 C++ 中存在越界访问 (xa[n]) 的高风险。
    ! F90 翻译同样会越界 (xa(n+1))。
    !
    ! Fortran (i) 对应 C++ [i-1]
    ! C++ [klo+1] -> F90 (klo+1+1) -> F90 (klo+2)
    ! C++ [khi+1] -> F90 (khi+1+1) -> F90 (khi+2)
    ! ======================================================================
    FUNCTION splini(xa, ya, y2a, n, x) RESULT(yi)
        ! INTEGRATE CUBIC SPLINE FUNCTION FROM XA(1) TO X
        ! XA, YA: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
        ! Y2A : ARRAY OF SECOND DERIVATIVES
        ! N : SIZE OF ARRAYS XA, YA, Y2A (按 C++ 0-based 习惯)
        ! X : ABSCISSA ENDPOINT FOR INTEGRATION
        ! Y : OUTPUT VALUE
        
        REAL(RL), DIMENSION(*), INTENT(IN) :: xa, ya, y2a
        INTEGER(IT), INTENT(IN) :: n
        REAL(RL), INTENT(IN) :: x
        REAL(RL) :: yi

        INTEGER(IT) :: klo, khi
        REAL(RL) :: xx, h, a, b, a2, b2
            
        yi = 0.0
        klo = 0   ! C++ 0-based 索引 0
        khi = 1   ! C++ 0-based 索引 1
        
        DO WHILE (x > xa(klo + 1) .AND. khi < n)
            xx = x
            IF (khi < n - 1) THEN
                IF (x < xa(khi + 1)) THEN
                    xx = x
                ELSE
                    xx = xa(khi + 1)
                END IF
            END IF
            
          
            h = xa(khi + 1) - xa(klo + 1)
            a = (xa(khi + 1) - xx) / h
            b = (xx - xa(klo + 1)) / h
            a2 = a * a
            b2 = b * b
            
          
            yi = yi + ( ((1.0 - a2) * ya(klo + 1) * 0.5 + &
                       b2 * ya(khi + 1) * 0.5 + &
                       ((-(1.0 + a2 * a2) * 0.25 + a2 * 0.5) * y2a(klo + 1) + &
                       (b2 * b2 * 0.25 - b2 * 0.5) * y2a(khi + 1)) * h * h / 6.0) * h )
            
            klo = klo + 1
            khi = khi + 1
        END DO
        
    END FUNCTION splini
    
    ! ======================================================================
    ! C++: void spline(const double* x, const double* y, int n, ...)
    ! F90: SUBROUTINE spline(x, y, n, yp1, ypn, y2)
    ! ======================================================================
    SUBROUTINE spline(x, y, n, yp1, ypn, y2)
        REAL(RL), DIMENSION(n), INTENT(IN)  :: x, y
        INTEGER(IT), INTENT(IN)  :: n
        REAL(RL), INTENT(IN)  :: yp1, ypn
        REAL(RL), DIMENSION(n), INTENT(OUT) :: y2
        
        REAL(RL), PARAMETER :: BIG = 0.99e30
        
        REAL(RL), ALLOCATABLE :: u(:) ! 临时工作数组
        REAL(RL) :: sig, p, qn, un
        INTEGER(IT) :: i, k
        
        ALLOCATE(u(n))
        
        ! 左边界

        IF (ABS(yp1) > BIG) THEN
            y2(1) = 0.0 
            u(1)  = 0.0 
        ELSE
            y2(1) = -0.5
            u(1)  = (3.0 / (x(2) - x(1))) * &
                    ((y(2) - y(1)) / (x(2) - x(1)) - yp1)
        END IF
        
        ! 三对角方程前向消元
        DO i = 2, n - 1

            sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
            p = sig * y2(i-1) + 2.0
            y2(i) = (sig - 1.0) / p
            u(i) = (6.0 * ((y(i+1) - y(i)) / (x(i+1) - x(i)) - &
                   (y(i) - y(i-1)) / (x(i) - x(i-1))) / &
                   (x(i+1) - x(i-1)) - &
                   sig * u(i-1)) / p
        END DO
        
        ! 右边界
        IF (ABS(ypn) >= BIG) THEN
            qn = 0.0
            un = 0.0
        ELSE
            qn = 0.5
            un = (3.0 / (x(n) - x(n-1))) * &
                 (ypn - (y(n) - y(n-1)) / (x(n) - x(n-1)))
        END IF
        

        y2(n) = (un - qn * u(n-1)) / (qn * y2(n-1) + 1.0)
        
        ! 回代
        DO k = n - 1, 1, -1
            y2(k) = y2(k) * y2(k+1) + u(k)
        END DO
        
        DEALLOCATE(u)
        
    END SUBROUTINE spline

    ! ======================================================================
    ! C++: double splint(const double* xa, const double* ya, const double* y2a, int n, double x)
    ! F90: FUNCTION splint(xa, ya, y2a, n, x)
    ! ======================================================================
    FUNCTION splint(xa, ya, y2a, n, x) RESULT(yi)
        ! CALCULATE CUBIC SPLINE INTERP VALUE
        
        REAL(RL), DIMENSION(n), INTENT(IN) :: xa, ya, y2a
        INTEGER(IT), INTENT(IN) :: n
        REAL(RL), INTENT(IN) :: x
        REAL(RL) :: yi
        
        INTEGER(IT) :: klo, khi, k
        REAL(RL) :: h, a, b
        
        ! C++: if (n < 2) throw ...
        IF (n < 2) THEN
            WRITE(*,*) 'splint: n must >= 2'
            STOP 1
        END IF
        
        klo = 0 
        khi = n-1 
        
        ! 二分查找区间
        ! C++: while (khi - klo > 1)
        DO WHILE (khi - klo > 1)
            k = (khi + klo) / 2 ! Fortran 整数除法
            IF (xa(k+1) > x) THEN
                khi = k
            ELSE
                klo = k
            END IF
        END DO

        h = xa(khi+1) - xa(klo+1)
        IF (h == 0.0) THEN
            WRITE(*,*) 'splint: bad xa input (duplicated x value)'
            STOP 2
        END IF
        
        a = (xa(khi+1) - x) / h
        b = (x - xa(klo+1)) / h
        
        ! C++: yi = a * ya[klo] + b * ya[khi] + ...
        yi = a * ya(klo+1) + b * ya(khi+1) + &
             ((a * a * a - a) * y2a(klo+1) + &
              (b * b * b - b) * y2a(khi+1)) * (h * h) / 6.0
              
    END FUNCTION splint

    ! ======================================================================
    ! C++: double zeta(const double& zz, double zl)
    ! ======================================================================
    FUNCTION zeta(zz, zl) RESULT(zeta_res)
        REAL(RL), INTENT(IN) :: zz, zl
        REAL(RL) :: zeta_res
        
        ! d_re 来自 nrlmsise_data 模块
        zeta_res = ((zz - zl) * (d_re + zl) / (d_re + zz))
        
    END FUNCTION zeta
    
    ! ======================================================================
    ! C++: double densm(...)
    ! ======================================================================
    FUNCTION densm(NRL, alt, d0, xm, tz, &
                   mn3, zn3, tn3, tgn3, &
                   mn2, zn2, tn2, tgn2) RESULT(densm_res)

        ! --- 输入/输出参数 ---
        TYPE(NRLMSISE), INTENT(IN)    :: NRL
        REAL(RL), INTENT(IN)  :: alt, d0, xm
        REAL(RL), INTENT(INOUT) :: tz ! C++ *tz 被读写
        INTEGER(IT), INTENT(IN)           :: mn3, mn2
        REAL(RL), DIMENSION(mn3), INTENT(IN) :: zn3, tn3
        REAL(RL), DIMENSION(2), INTENT(IN)   :: tgn3
        REAL(RL), DIMENSION(mn2), INTENT(IN) :: zn2, tn2
        REAL(RL), DIMENSION(2), INTENT(IN)   :: tgn2
        
        ! --- 返回值 ---
        REAL(RL) :: densm_res
        
        ! --- 局部变量 ---
        REAL(RL) :: z
        REAL(RL) :: z1, z2
        REAL(RL) :: t1, t2, zg, zgdif
        REAL(RL) :: yd1, yd2
        REAL(RL) :: x, y, yi
        REAL(RL) :: expl, gamm, glb
        REAL(RL) :: densm_tmp
        INTEGER(IT) :: mn, k
        
        ! 假设 mn2, mn3 总是 <= 10 (基于 C++ 固定大小数组)
        REAL(RL), DIMENSION(10) :: xs = 0.0
        REAL(RL), DIMENSION(10) :: ys = 0.0
        REAL(RL), DIMENSION(10) :: y2out = 0.0

        densm_tmp = d0

        ! C++: if (alt > zn2[0])
        IF (alt > zn2(1)) THEN
            IF (xm == 0.0) THEN
                ! C++: return *tz; (通过读取 tz)
                densm_res = tz
                RETURN
            ELSE
                ! C++: return d0;
                densm_res = d0
                RETURN
            END IF
        END IF

        ! STRATOSPHERE / MESOSPHERE TEMPERATURE
        ! C++: if (alt > zn2[mn2 - 1])
        IF (alt > zn2(mn2)) THEN
             z = alt
        ELSE
             z = zn2(mn2)
        END IF

        mn = mn2
        z1 = zn2(1)         ! C++ zn2[0]
        z2 = zn2(mn)        ! C++ zn2[mn-1]
        t1 = tn2(1)         ! C++ tn2[0]
        t2 = tn2(mn)        ! C++ tn2[mn-1]
        zg = zeta(z, z1)
        zgdif = zeta(z2, z1)

        ! C++: for (int k = 0; k < mn; k++)
        DO k = 1, mn
            xs(k) = zeta(zn2(k), z1) / zgdif
            ys(k) = 1.0 / tn2(k)
        END DO
        
        yd1 = -tgn2(1) / (t1 * t1) * zgdif
        yd2 = -tgn2(2) / (t2 * t2) * zgdif * (((d_re + z2) / (d_re + z1))**2.0)
        
        ! alculate spline coefficients
        ! 假设 spline 接口为 DIMENSION(*)
        CALL spline(xs, ys, mn, yd1, yd2, y2out)
        
        x = zg / zgdif
        ! 假设 splint 接口为 DIMENSION(*)
        y = splint(xs, ys, y2out, mn, x)

        ! C++: *tz = 1.0 / y; (写入 tz)
        tz = 1.0 / y
        
        IF (xm /= 0.0) THEN
            ! calaculate stratosphere / mesospehere density
            glb = NRL%d_gsurf / ((1.0 + z1 / d_re)**2.0)
            gamm = xm * glb * zgdif / d_RGAS

            ! Integrate temperature profile
            ! 假设 splini 接口为 DIMENSION(*)
            yi = splini(xs, ys, y2out, mn, x)
            expl = gamm * yi
            IF (expl > 50.0) THEN
                expl = 50.0
            END IF
            ! Density at altitude
            densm_tmp = densm_tmp * (t1 / tz) * EXP(-expl)
        END IF

        ! C++: if (alt > zn3[0])
        IF (alt > zn3(1)) THEN
            IF (xm == 0.0) THEN
                densm_res = tz
                RETURN
            ELSE
                densm_res = densm_tmp
                RETURN
            END IF
        END IF

        ! troposhere / stratosphere temperature
        z = alt
        mn = mn3
        z1 = zn3(1)         ! C++ zn3[0]
        z2 = zn3(mn)        ! C++ zn3[mn-1]
        t1 = tn3(1)         ! C++ tn3[0]
        t2 = tn3(mn)        ! C++ tn3[mn-1]
        zg = zeta(z, z1)
        zgdif = zeta(z2, z1)

        ! set up spline nodes
        ! C++: for (k = 0; k < mn; k++)
        DO k = 1, mn
            xs(k) = zeta(zn3(k), z1) / zgdif
            ys(k) = 1.0 / tn3(k)
        END DO
        

        yd1 = -tgn3(1) / (t1 * t1) * zgdif
        yd2 = -tgn3(2) / (t2 * t2) * zgdif * (((d_re + z2) / (d_re + z1))**2.0)

        ! calculate spline coefficients
        CALL spline(xs, ys, mn, yd1, yd2, y2out)
        
        x = zg / zgdif
        y = splint(xs, ys, y2out, mn, x)

        ! temperature at altitude
        ! C++: *tz = 1.0 / y;
        tz = 1.0 / y
        
        IF (xm /= 0.0) THEN
            ! calaculate tropospheric / stratosphere density
            glb = NRL%d_gsurf / ((1.0 + z1 / d_re)**2.0)
            gamm = xm * glb * zgdif / d_RGAS

            ! Integrate temperature profile
            yi = splini(xs, ys, y2out, mn, x)
            expl = gamm * yi
            IF (expl > 50.0) THEN
                expl = 50.0
            END IF

            ! Density at altitude
            densm_tmp = densm_tmp * (t1 / tz) * EXP(-expl)
        END IF

        ! Final return
        IF (xm == 0.0) THEN
            ! C++: return *tz;
            densm_res = tz
        ELSE
            ! C++: return densm_tmp;
            densm_res = densm_tmp
        END IF
        
END FUNCTION densm

! ======================================================================
    ! C++: double  densu(...)
    ! ======================================================================
    FUNCTION densu(alt, dlb, tinf, tlb, xm, alpha, tz, &
                   zlb, s2, mn1, zn1, tn1, tgn1) RESULT(densu_res)

        ! --- 参数 (Arguments) ---
        REAL(RL), INTENT(IN)    :: alt, dlb, tinf, tlb, xm, alpha, zlb, s2
        REAL(RL), INTENT(INOUT) :: tz  ! C++ *tz (读/写)
        INTEGER(IT), INTENT(IN) :: mn1
        REAL(RL), DIMENSION(mn1), INTENT(IN)    :: zn1
        REAL(RL), DIMENSION(mn1), INTENT(INOUT) :: tn1  ! C++ *tn1 (被写入)
        REAL(RL), DIMENSION(2),   INTENT(INOUT) :: tgn1 ! C++ *tgn1 (被写入 tgn1[0])

        ! --- 返回值 ---
        REAL(RL) :: densu_res

        ! --- 局部变量 ---
        REAL(RL) :: yd2, yd1, x, y
        REAL(RL) :: densu_temp
        REAL(RL) :: za, z1, z2
        REAL(RL) :: z, zg2, tt, ta
        REAL(RL) :: dta, t1, t2, zg, zgdif
        INTEGER(IT) :: mn, k
        REAL(RL) :: glb
        REAL(RL) :: expl
        REAL(RL) :: yi
        REAL(RL) :: densa
        REAL(RL) :: gamma, gamm
        
        ! 临时数组 (C++ 中固定大小为 5)
        REAL(RL), DIMENSION(5) :: xs, ys, y2out

        ! C++ 局部变量，遮蔽了全局 d_re 和 d_gsurf
        REAL(RL) :: local_d_gsurf, local_d_re

        ! --- 外部函数 (假定在同一模块中) ---
        !REAL(RL) :: zeta, splint, splini
        
        ! --- 初始化 ---
        x = 0.0_RL
        densu_temp = 1.0_RL
        z1 = 0.0_RL
        t1 = 0.0_RL
        zgdif = 0.0_RL
        mn = 0
        xs = 0.0_RL
        ys = 0.0_RL
        y2out = 0.0_RL
        local_d_gsurf = 0.0_RL
        local_d_re = 0.0_RL

        ! --- 函数体 ---
        
        ! joining altitudes of Bates and spline
        za = zn1(1)
        
        z = MAX(alt, za)

        ! geopotential altitude difference from ZLB
        zg2 = zeta(z, zlb)

        ! Bates temperature
        tt = tinf - (tinf - tlb) * EXP(-s2 * zg2)
        ta = tt
        tz = tt 
        densu_temp = tz 

        IF (alt < za) THEN
            ! calculate temperature below ZA
            ! temperature gradient at ZA from Bates profile
            
            dta = (tinf - ta) * s2 * ((d_re + zlb) / (d_re + za))**2.0_RL
            
            tgn1(1) = dta 
            tn1(1) = ta  
            
            IF (alt > zn1(mn1)) THEN
                z = alt
            ELSE
                z = zn1(mn1)
            END IF
            
            mn = mn1
            z1 = zn1(1)     
            z2 = zn1(mn)   
            t1 = tn1(1)     
            t2 = tn1(mn)   
            
            ! geopotental difference from z1
            zg = zeta(z, z1)
            zgdif = zeta(z2, z1)
            
            ! set up spline nodes
            DO k = 1, mn
                xs(k) = zeta(zn1(k), z1) / zgdif
                ys(k) = 1.0_RL / tn1(k)
            END DO
            
            ! end node derivatives
            yd1 = -tgn1(1) / (t1 * t1) * zgdif
            yd2 = -tgn1(2) / (t2 * t2) * zgdif * ((d_re + z2) / (d_re + z1))**2.0_RL
            
            ! calculate spline coefficients
            CALL spline(xs, ys, mn, yd1, yd2, y2out)
            
            x = zg / zgdif
            y = splint(xs, ys, y2out, mn, x)
            
            ! temperature at altitude
            tz = 1.0_RL / y !
            densu_temp = tz
        END IF

        IF (xm == 0.0_RL) THEN
            densu_res = densu_temp
            RETURN
        END IF

        ! calculate density above za
        
        CALL glatf(0.0_RL, local_d_gsurf, local_d_re)
        
        glb = local_d_gsurf / (1.0_RL + zlb / local_d_re)**2.0_RL
        
        gamma = xm * glb / (s2 * d_RGAS * tinf)
        
        expl = EXP(-s2 * gamma * zg2)
        IF (expl > 50.0_RL) expl = 50.0_RL
        IF (tt <= 0.0_RL) expl = 50.0_RL

        ! density at altitude
        densa = dlb * (tlb / tt)**(1.0_RL + alpha + gamma) * expl
        densu_temp = densa
        
        IF (alt >= za) THEN
            densu_res = densu_temp
            RETURN
        END IF

        ! calculate density below za
        
        glb = local_d_gsurf / (1.0_RL + z1 / local_d_re)**2.0_RL
        
        ! C++: gamm = xm * glb * zgdif / d_RGAS; (使用全局 d_RGAS)
        gamm = xm * glb * zgdif / d_RGAS

        ! integrate spline temperatures
        yi = splini(xs, ys, y2out, mn, x)
        expl = gamm * yi
        IF (expl > 50.0_RL) expl = 50.0_RL
        
        IF (tz <= 0.0_RL) expl = 50.0_RL

        ! density at altitude
        densu_temp = densu_temp * (t1 / tz)**(1.0_RL + alpha) * EXP(-expl)
        
        densu_res = densu_temp
        
    END FUNCTION densu                   
                   
                   
                   
                   
    FUNCTION sumex(ex) RESULT(sumex_res)
        REAL(RL), INTENT(IN) :: ex
        REAL(RL) :: sumex_res
        
        sumex_res = 1.0 + (1.0 - ex**19.0) / (1.0 - ex) * (ex**0.5)
        
    END FUNCTION sumex

    ! ======================================================================
    ! C++: inline double g0(double a, const double *p)
    ! ======================================================================
    FUNCTION g0(a, p) RESULT(g0_res)
        REAL(RL), INTENT(IN) :: a
        REAL(RL), DIMENSION(*), INTENT(IN) :: p
        REAL(RL) :: g0_res
        
        REAL(RL) :: p25_abs
        
        ! C++ p[24] -> F90 p(25)
        ! C++ sqrt(p[24]*p[24]) 是一种计算 ABS(p[24]) 的方法
        p25_abs = SQRT(p(25)**2.0)
        
        ! C++ p[25] -> F90 p(26)
        g0_res = (a - 4.0 + (p(26) - 1.0) * (a - 4.0 + &
                 (EXP(-p25_abs * (a - 4.0)) - 1.0) / p25_abs))
                 
    END FUNCTION g0

    ! ======================================================================
    ! C++: inline double sg0(double ex, const double *p, double *ap)
    ! ======================================================================
    FUNCTION sg0(ex, p, ap) RESULT(sg0_res)
        REAL(RL), INTENT(IN) :: ex
        REAL(RL), DIMENSION(*), INTENT(IN) :: p, ap
        REAL(RL) :: sg0_res
        
        ! 函数 g0 和 sumex 必须在此模块中可用
        !REAL(RL) :: g0, sumex
        
        sg0_res = (g0(ap(2), p) + (g0(ap(3), p) * ex + g0(ap(4), p) * ex**2.0 + &
                   g0(ap(5), p) * ex**3.0 + (g0(ap(6), p) * ex**4.0 + &
                   g0(ap(7), p) * ex**12.0) * (1.0 - ex**8.0) / (1.0 - ex))) &
                   / sumex(ex)
                   
    END FUNCTION sg0
    
! ======================================================================
    ! C++: double globe7(...)
    ! ======================================================================
    FUNCTION globe7(NRL, p, doy, sec, g_lat, g_long, lst, f107A, f107, ap) RESULT(tinf)
        ! --- 函数参数 ---
        TYPE(NRLMSISE), INTENT(INOUT) :: NRL ! INOUT 因为 a_plg 和 d_... 成员被修改
        REAL(RL), DIMENSION(*), INTENT(IN) :: p, ap
        INTEGER(IT), INTENT(IN) :: doy
        REAL(RL), INTENT(IN) :: sec, g_lat, g_long, lst, f107A, f107
        
        ! --- 返回值 ---
        REAL(RL) :: tinf
        
        ! --- 局部变量 ---
        REAL(RL), DIMENSION(15) :: t ! C++ t[15] (0..14) -> F90 t(15) (1..15)
        INTEGER(IT) :: i, j
        REAL(RL) :: apd
        REAL(RL) :: tloc
        REAL(RL) :: c, s, c2, c4, s2
        REAL(RL) :: cd32, cd18, cd14, cd39
        REAL(RL) :: df
        REAL(RL) :: f1, f2
        REAL(RL) :: t71, t72
        REAL(RL) :: t81, t82
        REAL(RL) :: exp1
        REAL(RL) :: p44, p45

        ! --- 外部函数 ---
        !REAL(RL) :: sg0
        
        tloc = lst

        ! C++: for (int j = 0; j < 14; j++) { t[j] = 0; }
        ! C++ 索引 0..13 -> F90 索引 1..14
        t(1:14) = 0.0
        t(15) = 0.0 ! C++ t[14]

        c = SIN(g_lat * d_DGTR)
        s = COS(g_lat * d_DGTR)
        c2 = c * c
        c4 = c2 * c2
        s2 = s * s

        ! C++: NRL.a_plg[0][1] = c; -> F90: NRL%a_plg(1, 2)
        NRL%a_plg(1, 2) = c
        NRL%a_plg(1, 3) = 0.5 * (3.0 * c2 - 1.0)
        NRL%a_plg(1, 4) = 0.5 * (5.0 * c * c2 - 3.0 * c)
        NRL%a_plg(1, 5) = (35.0 * c4 - 30.0 * c2 + 3.0) / 8.0
        NRL%a_plg(1, 6) = (63.0 * c2 * c2 * c - 70.0 * c2 * c + 15.0 * c) / 8.0
        NRL%a_plg(1, 7) = (11.0 * c * NRL%a_plg(1, 6) - 5.0 * NRL%a_plg(1, 5)) / 6.0
        
        NRL%a_plg(2, 2) = s
        NRL%a_plg(2, 3) = 3.0 * c * s
        NRL%a_plg(2, 4) = 1.5 * (5.0 * c2 - 1.0) * s
        NRL%a_plg(2, 5) = 2.5 * (7.0 * c2 * c - 3.0 * c) * s
        NRL%a_plg(2, 6) = 1.875 * (21.0 * c4 - 14.0 * c2 + 1.0) * s
        NRL%a_plg(2, 7) = (11.0 * c * NRL%a_plg(2, 6) - 6.0 * NRL%a_plg(2, 5)) / 5.0
        
        NRL%a_plg(3, 3) = 3.0 * s2
        NRL%a_plg(3, 4) = 15.0 * s2 * c
        NRL%a_plg(3, 5) = 7.5 * (7.0 * c2 - 1.0) * s2
        NRL%a_plg(3, 6) = 3.0 * c * NRL%a_plg(3, 5) - 2.0 * NRL%a_plg(3, 4)
        NRL%a_plg(3, 7) = (11.0 * c * NRL%a_plg(3, 6) - 7.0 * NRL%a_plg(3, 5)) / 4.0
        NRL%a_plg(3, 8) = (13.0 * c * NRL%a_plg(3, 7) - 8.0 * NRL%a_plg(3, 6)) / 5.0
        
        NRL%a_plg(4, 4) = 15.0 * s2 * s
        NRL%a_plg(4, 5) = 105.0 * s2 * s * c
        NRL%a_plg(4, 6) = (9.0 * c * NRL%a_plg(4, 5) - 7.0 * NRL%a_plg(4, 4)) / 2.0
        NRL%a_plg(4, 7) = (11.0 * c * NRL%a_plg(4, 6) - 8.0 * NRL%a_plg(4, 5)) / 3.0

        ! C++: if (!(((NRL.a_sw[7] == 0) && (NRL.a_sw[8] == 0)) && (NRL.a_sw[14] == 0)))
        ! (a_sw(8), a_sw(9), a_sw(15) are non-zero)
        IF ((NRL%a_sw(8) /= 0.0) .OR. (NRL%a_sw(9) /= 0.0) .OR. (NRL%a_sw(15) /= 0.0)) THEN
            NRL%d_stloc  = SIN(d_HR * tloc)
            NRL%d_ctloc  = COS(d_HR * tloc)
            NRL%d_s2tloc = SIN(2.0 * d_HR * tloc)
            NRL%d_c2tloc = COS(2.0 * d_HR * tloc)
            NRL%d_s3tloc = SIN(3.0 * d_HR * tloc)
            NRL%d_c3tloc = COS(3.0 * d_HR * tloc)
        END IF


        cd32 = COS(d_DR * (doy - p(32)))
        cd18 = COS(2.0 * d_DR * (doy - p(18)))
        cd14 = COS(d_DR * (doy - p(14)))
        cd39 = COS(2.0 * d_DR * (doy - p(39)))

        ! F10.7 EFFECT
        df = f107 - f107A
        NRL%d_dfa = f107A - 150.0

        t(1) = p(20) * df * (1.0 + p(60) * NRL%d_dfa) + p(21) * df * df + p(22) * NRL%d_dfa + p(30) * NRL%d_dfa**2.0
        f1 = 1.0 + (p(48) * NRL%d_dfa + p(20) * df + p(21) * df * df) * NRL%a_swc(2)
        f2 = 1.0 + (p(50) * NRL%d_dfa + p(20) * df + p(21) * df * df) * NRL%a_swc(2)

        ! TIME INDEPENDENT
        t(2) = (p(2) * NRL%a_plg(1, 3) + p(3) * NRL%a_plg(1, 5) + p(23) * NRL%a_plg(1, 7)) + &
               (p(15) * NRL%a_plg(1, 3)) * NRL%d_dfa * NRL%a_swc(2) + p(27) * NRL%a_plg(1, 2)

        ! SYMMETRICAL ANNUAL
        t(3) = p(19) * cd32

        ! SYMMETRICAL SEMIANNUAL
        t(4) = (p(16) + p(17) * NRL%a_plg(1, 3)) * cd18

        ! ASYMMETRICAL ANNUAL
        t(5) = f1 * (p(10) * NRL%a_plg(1, 2) + p(11) * NRL%a_plg(1, 4)) * cd14

        ! ASYMMETRICAL SEMIANNUAL
        t(6) = p(38) * NRL%a_plg(1, 2) * cd39
        
        ! DIURNAL
        IF (NRL%a_sw(8) /= 0.0) THEN

            t71 = (p(12) * NRL%a_plg(2, 3)) * cd14 * NRL%a_swc(6)
            t72 = (p(13) * NRL%a_plg(2, 3)) * cd14 * NRL%a_swc(6)
            t(7) = f2 * ((p(4) * NRL%a_plg(2, 2) + p(5) * NRL%a_plg(2, 4) + p(28) * NRL%a_plg(2, 6) + t71) * &
                 NRL%d_ctloc + (p(7) * NRL%a_plg(2, 2) + p(8) * NRL%a_plg(2, 4) + p(29) * NRL%a_plg(2, 6) &
                 + t72) * NRL%d_stloc)
        END IF
        
        ! SEMIDIURNAL
        IF (NRL%a_sw(9) /= 0.0) THEN
            t81 = (p(24) * NRL%a_plg(3, 4) + p(36) * NRL%a_plg(3, 6)) * cd14 * NRL%a_swc(6)
            t82 = (p(34) * NRL%a_plg(3, 4) + p(37) * NRL%a_plg(3, 6)) * cd14 * NRL%a_swc(6)
            t(8) = f2 * ((p(6) * NRL%a_plg(3, 3) + p(42) * NRL%a_plg(3, 5) + t81) * NRL%d_c2tloc + &
                 (p(9) * NRL%a_plg(3, 3) + p(43) * NRL%a_plg(3, 5) + t82) * NRL%d_s2tloc)
        END IF

        ! TERDIURNAL
        IF (NRL%a_sw(15) /= 0.0) THEN
            t(14) = f2 * ((p(40) * NRL%a_plg(4, 4) + (p(94) * NRL%a_plg(4, 5) + p(47) * NRL%a_plg(4, 7)) * cd14 * NRL%a_swc(6)) * &
                 NRL%d_s3tloc + (p(41) * NRL%a_plg(4, 4) + (p(95) * NRL%a_plg(4, 5) + p(49) * NRL%a_plg(4, 7)) * cd14 * NRL%a_swc(6)) * NRL%d_c3tloc)
        END IF

        ! magnetic activity based on daily ap
        IF (NRL%a_sw(10) == -1.0) THEN
            IF (p(52) /= 0.0) THEN
                exp1 = EXP(-10800.0 * SQRT(p(52)**2.0) / &
                       (1.0 + p(139) * (45.0 - SQRT(g_lat**2.0))))
                IF (exp1 > 0.99999) THEN
                    exp1 = 0.99999
                END IF
                NRL%a_apt(1) = sg0(exp1, p, ap)
                
                IF (NRL%a_sw(10) /= 0.0) THEN
                    t(9) = NRL%a_apt(1) * (p(51) + p(97) * NRL%a_plg(1, 3) + p(55) * NRL%a_plg(1, 5) + &
                         (p(126) * NRL%a_plg(1, 2) + p(127) * NRL%a_plg(1, 4) + p(128) * NRL%a_plg(1, 6)) * cd14 * NRL%a_swc(6) + &
                         (p(129) * NRL%a_plg(2, 2) + p(130) * NRL%a_plg(2, 4) + p(131) * NRL%a_plg(2, 6)) * NRL%a_swc(8) * &
                         COS(d_HR * (tloc - p(132))))
                END IF
            END IF
        ELSE
            apd = ap(1) - 4.0
            p44 = p(44)
            p45 = p(45)
            IF (p44 < 0.0) THEN
                p44 = 1.0E-5
            END IF
            NRL%d_apdf = apd + (p45 - 1.0) * (apd + (EXP(-p44 * apd) - 1.0) / p44)
            
            IF (NRL%a_sw(10) /= 0.0) THEN
                t(9) = NRL%d_apdf * (p(33) + p(46) * NRL%a_plg(1, 3) + p(35) * NRL%a_plg(1, 5) + &
                     (p(101) * NRL%a_plg(1, 2) + p(102) * NRL%a_plg(1, 4) + p(103) * NRL%a_plg(1, 6)) * cd14 * NRL%a_swc(6) + &
                     (p(122) * NRL%a_plg(2, 2) + p(123) * NRL%a_plg(2, 4) + p(124) * NRL%a_plg(2, 6)) * NRL%a_swc(8) * &
                     COS(d_HR * (tloc - p(125))))
            END IF
        END IF

        IF ((NRL%a_sw(11) /= 0.0) .AND. (g_long > -1000.0)) THEN
            
            ! longitudinal
            IF (NRL%a_sw(12) /= 0.0) THEN
                t(11) = (1.0 + p(81) * NRL%d_dfa * NRL%a_swc(2)) * &
                      ((p(65) * NRL%a_plg(2, 3) + p(66) * NRL%a_plg(2, 5) + p(67) * NRL%a_plg(2, 7) + &
                      p(104) * NRL%a_plg(2, 2) + p(105) * NRL%a_plg(2, 4) + p(106) * NRL%a_plg(2, 6) + &
                      NRL%a_swc(6) * (p(110) * NRL%a_plg(2, 2) + p(111) * NRL%a_plg(2, 4) + p(112) * NRL%a_plg(2, 6)) * cd14) * &
                      COS(d_DGTR * g_long) + &
                      (p(91) * NRL%a_plg(2, 3) + p(92) * NRL%a_plg(2, 5) + p(93) * NRL%a_plg(2, 7) + &
                      p(107) * NRL%a_plg(2, 2) + p(108) * NRL%a_plg(2, 4) + p(109) * NRL%a_plg(2, 6) + &
                      NRL%a_swc(6) * (p(113) * NRL%a_plg(2, 2) + p(114) * NRL%a_plg(2, 4) + p(115) * NRL%a_plg(2, 6)) * cd14) * &
                      SIN(d_DGTR * g_long))
            END IF

            ! ut and mixed ut, longitude
            IF (NRL%a_sw(13) /= 0.0) THEN
                t(12) = (1.0 + p(96) * NRL%a_plg(1, 2)) * (1.0 + p(82) * NRL%d_dfa * NRL%a_swc(2)) * &
                       (1.0 + p(120) * NRL%a_plg(1, 2) * NRL%a_swc(6) * cd14) * &
                       ((p(69) * NRL%a_plg(1, 2) + p(70) * NRL%a_plg(1, 4) + p(71) * NRL%a_plg(1, 6)) * &
                       COS(d_SR * (sec - p(72))))
                
                t(12) = t(12) + NRL%a_swc(12) * &
                       (p(77) * NRL%a_plg(3, 4) + p(78) * NRL%a_plg(3, 6) + p(79) * NRL%a_plg(3, 8)) * &
                       COS(d_SR * (sec - p(80)) + 2.0 * d_DGTR * g_long) * (1.0 + p(138) * NRL%d_dfa * NRL%a_swc(2))
            END IF

            ! ut, longitude magnetic activity
            IF (NRL%a_sw(14) /= 0.0) THEN
                IF (NRL%a_sw(10) == -1.0) THEN
                    IF (p(52) /= 0.0) THEN
                        t(13) = NRL%a_apt(1) * NRL%a_swc(12) * (1.0 + p(133) * NRL%a_plg(1, 2)) * &
                              ((p(53) * NRL%a_plg(2, 3) + p(99) * NRL%a_plg(2, 5) + p(68) * NRL%a_plg(2, 7)) * &
                              COS(d_DGTR * (g_long - p(98)))) + &
                              NRL%a_apt(1) * NRL%a_swc(12) * NRL%a_swc(6) * &
                              (p(134) * NRL%a_plg(2, 2) + p(135) * NRL%a_plg(2, 4) + p(136) * NRL%a_plg(2, 6)) * &
                              cd14 * COS(d_DGTR * (g_long - p(137))) + &
                              NRL%a_apt(1) * NRL%a_swc(13) * &
                              (p(56) * NRL%a_plg(1, 2) + p(57) * NRL%a_plg(1, 4) + p(58) * NRL%a_plg(1, 6)) * &
                              COS(d_SR * (sec - p(59)))
                    END IF
                ELSE
                    t(13) = NRL%d_apdf * NRL%a_swc(12) * (1.0 + p(121) * NRL%a_plg(1, 2)) * &
                          ((p(61) * NRL%a_plg(2, 3) + p(62) * NRL%a_plg(2, 5) + p(63) * NRL%a_plg(2, 7)) * &
                          COS(d_DGTR * (g_long - p(64)))) + &
                          NRL%d_apdf * NRL%a_swc(12) * NRL%a_swc(6) * &
                          (p(116) * NRL%a_plg(2, 2) + p(117) * NRL%a_plg(2, 4) + p(118) * NRL%a_plg(2, 6)) * &
                          cd14 * COS(d_DGTR * (g_long - p(119))) + &
                          NRL%d_apdf * NRL%a_swc(13) * &
                          (p(84) * NRL%a_plg(1, 2) + p(85) * NRL%a_plg(1, 4) + p(86) * NRL%a_plg(1, 6)) * &
                          COS(d_SR * (sec - p(76)))
                END IF
            END IF
        END IF
        
        tinf = p(31)
        DO i = 1, 14

            tinf = tinf + ABS(NRL%a_sw(i + 1)) * t(i)
        END DO

    END FUNCTION globe7
                   
    FUNCTION glob7s(NRL, p, doy, g_long) RESULT(tt)
        ! --- 函数参数 ---
        TYPE(NRLMSISE), INTENT(IN) :: NRL ! INTENT(IN) 对应 C++ 按值传递
        REAL(RL), DIMENSION(*), INTENT(IN) :: p
        INTEGER(IT), INTENT(IN) :: doy
        REAL(RL), INTENT(IN) :: g_long
        
        ! --- 返回值 ---
        REAL(RL) :: tt

        ! --- 局部变量 ---
        REAL(RL), DIMENSION(14) :: t ! C++ t[14] (0..13)
        REAL(RL) :: cd32, cd18, cd14, cd39
        INTEGER(IT) :: i, j
        REAL(RL) :: t71, t72
        REAL(RL) :: t81, t82

        ! --- 函数体 ---
        
        t(1:14) = 0.0


        cd32 = COS(d_DR * (doy - p(32)))
        cd18 = COS(2.0 * d_DR * (doy - p(18)))
        cd14 = COS(d_DR * (doy - p(14)))
        cd39 = COS(2.0 * d_DR * (doy - p(39)))

        ! F10.7
        t(1) = p(22) * NRL%d_dfa

        ! time independent
        t(2) = p(2) * NRL%a_plg(1, 3) + p(3) * NRL%a_plg(1, 5) + p(23) * NRL%a_plg(1, 7) + &
               p(27) * NRL%a_plg(1, 2) + p(15) * NRL%a_plg(1, 4) + p(60) * NRL%a_plg(1, 6)

        ! SYMMETRICAL ANNUAL
        t(3) = (p(19) + p(48) * NRL%a_plg(1, 3) + p(30) * NRL%a_plg(1, 5)) * cd32

        ! SYMMETRICAL SEMIANNUAL
        t(4) = (p(16) + p(17) * NRL%a_plg(1, 3) + p(31) * NRL%a_plg(1, 5)) * cd18

        ! ASYMMETRICAL ANNUAL
        t(5) = (p(10) * NRL%a_plg(1, 2) + p(11) * NRL%a_plg(1, 4) + p(21) * NRL%a_plg(1, 6)) * cd14

        ! ASYMMETRICAL SEMIANNUAL
        t(6) = (p(38) * NRL%a_plg(1, 2)) * cd39

        ! DIURNAL
        IF (NRL%a_sw(8) /= 0.0) THEN
            t71 = p(12) * NRL%a_plg(2, 3) * cd14 * NRL%a_swc(6)
            t72 = p(13) * NRL%a_plg(2, 3) * cd14 * NRL%a_swc(6)
            t(7) = ((p(4) * NRL%a_plg(2, 2) + p(5) * NRL%a_plg(2, 4) + t71) * NRL%d_ctloc + &
                   (p(7) * NRL%a_plg(2, 2) + p(8) * NRL%a_plg(2, 4) + t72) * NRL%d_stloc)
        END IF

        ! SEMIDIURNAL
        ! C++ if (NRL.a_sw[8])
        IF (NRL%a_sw(9) /= 0.0) THEN
            t81 = (p(24) * NRL%a_plg(3, 4) + p(36) * NRL%a_plg(3, 6)) * cd14 * NRL%a_swc(6)
            t82 = (p(34) * NRL%a_plg(3, 4) + p(37) * NRL%a_plg(3, 6)) * cd14 * NRL%a_swc(6)
            t(8) = ((p(6) * NRL%a_plg(3, 3) + p(42) * NRL%a_plg(3, 5) + t81) * NRL%d_c2tloc + &
                   (p(9) * NRL%a_plg(3, 3) + p(43) * NRL%a_plg(3, 5) + t82) * NRL%d_s2tloc)
        END IF

        ! TERDIURNAL
        IF (NRL%a_sw(15) /= 0.0) THEN
            t(14) = p(40) * NRL%a_plg(4, 4) * NRL%d_s3tloc + p(41) * NRL%a_plg(4, 4) * NRL%d_c3tloc
        END IF

        ! MAGNETIC ACTIVITY
        ! C++ if (NRL.a_sw[9])
        IF (NRL%a_sw(10) /= 0.0) THEN
            IF (NRL%a_sw(10) == 1.0) THEN
                t(9) = NRL%d_apdf * (p(33) + p(46) * NRL%a_plg(1, 3) * NRL%a_swc(3))
            END IF
            ! C++ if (NRL.a_sw[9] == -1)
            IF (NRL%a_sw(10) == -1.0) THEN
                t(9) = (p(51) * NRL%a_apt(1) + p(97) * NRL%a_plg(1, 3) * NRL%a_apt(1) * NRL%a_swc(3))
            END IF
        END IF
        
        ! LONGITUDINAL
        ! C++ if (!((NRL.a_sw[10] == 0) || (NRL.a_sw[11] == 0) || (g_long <= -1000.0)))
        IF ((NRL%a_sw(11) /= 0.0) .AND. (NRL%a_sw(12) /= 0.0) .AND. (g_long > -1000.0)) THEN
            t(11) = (1.0 + NRL%a_plg(1, 2) * (p(81) * NRL%a_swc(6) * COS(d_DR * (doy - p(82))) &
                 + p(86) * NRL%a_swc(7) * COS(2.0 * d_DR * (doy - p(87)))) &
                 + p(84) * NRL%a_swc(4) * COS(d_DR * (doy - p(85))) &
                 + p(88) * NRL%a_swc(5) * COS(2.0 * d_DR * (doy - p(89)))) &
                 * ((p(65) * NRL%a_plg(2, 3) + p(66) * NRL%a_plg(2, 5) + p(67) * NRL%a_plg(2, 7) &
                 + p(75) * NRL%a_plg(2, 2) + p(76) * NRL%a_plg(2, 4) + p(77) * NRL%a_plg(2, 6)) * COS(d_DGTR * g_long) &
                 + (p(91) * NRL%a_plg(2, 3) + p(92) * NRL%a_plg(2, 5) + p(93) * NRL%a_plg(2, 7) &
                 + p(78) * NRL%a_plg(2, 2) + p(79) * NRL%a_plg(2, 4) + p(80) * NRL%a_plg(2, 6)) * SIN(d_DGTR * g_long))
        END IF

        tt = 0.0

        DO i = 1, 14
            tt = tt + ABS(NRL%a_sw(i + 1)) * t(i)
        END DO
        
        ! C++: return tt; (由 RESULT(tt) 自动处理)
        
END FUNCTION glob7s     



SUBROUTINE gts7(NRL, doy, sec, alt, g_lat, g_long, lst, f107A, f107, ap, d, t)
        
        ! --- 参数 (Arguments) ---
        TYPE(NRLMSISE), INTENT(INOUT) :: NRL
        INTEGER(IT), INTENT(IN) :: doy
        REAL(RL), INTENT(IN) :: sec, alt, g_lat, g_long, lst, f107A, f107
        REAL(RL), DIMENSION(*), INTENT(IN) :: ap
        REAL(RL), DIMENSION(*), INTENT(OUT) :: d  ! C++ d[9] (0..8)
        REAL(RL), DIMENSION(*), INTENT(OUT) :: t  ! C++ t[2] (0..1)
        
        ! --- 局部变量 ---
        REAL(RL) :: za
        INTEGER(IT) :: i, j
        REAL(RL) :: z
        REAL(RL), DIMENSION(5) :: zn1 = (/ 120.0, 110.0, 100.0, 90.0, 72.5 /)
        REAL(RL) :: tinf
        INTEGER(IT), PARAMETER :: mn1 = 5
        REAL(RL) :: g0
        REAL(RL) :: tlb
        REAL(RL) :: s
        REAL(RL) :: db01, db04, db14, db16, db28, db32, db40
        REAL(RL) :: zh28, zh04, zh16, zh32, zh40, zh01, zh14
        REAL(RL) :: zhm28, zhm04, zhm16, zhm32, zhm40, zhm01, zhm14
        REAL(RL) :: xmd
        REAL(RL) :: b28, b04, b16, b32, b40, b01, b14
        REAL(RL) :: g28, g4, g16, g32, g40, g1, g14
        REAL(RL) :: zhf, xmm
        REAL(RL) :: zc04, zc16, zc32, zc40, zc01, zc14
        REAL(RL) :: hc04, hc16, hc32, hc40, hc01, hc14
        REAL(RL) :: hcc16, hcc32, hcc01, hcc14
        REAL(RL) :: zcc16, zcc32, zcc01, zcc14
        REAL(RL) :: rc16, rc32, rc01, rc14
        REAL(RL) :: rll
        REAL(RL) :: g16h, db16h, zsho
        REAL(RL) :: tho, zsht, zmho
        REAL(RL) :: dd
        REAL(RL) :: hc216, hcc232
        

        REAL(RL), DIMENSION(9), PARAMETER :: alpha = &
            (/ -0.38, 0.0, 0.0, 0.0, 0.17, 0.0, -0.38, 0.0, 0.0 /)

        REAL(RL), DIMENSION(8), PARAMETER :: altl = &
            (/ 200.0, 300.0, 160.0, 250.0, 240.0, 450.0, 320.0, 450.0 /)

        ! --- 外部函数 ---
        !REAL(RL) :: globe7, glob7s, densu, dnet, ccor, ccor2, scalh, zeta
        
        ! --- 函数体 ---
        
        za = pdl(2, 16)
        zn1(1) = za
        
        d(1:9) = 0.0

        ! TINF VARIATIONS NOT IMPORTANT BELOW ZA OR ZN1(1)
        IF (alt > zn1(1)) THEN
            tinf = ptm(1) * pt(1) * (1.0 + NRL%a_sw(17) * &
                   globe7(NRL, pt, doy, sec, g_lat, g_long, lst, f107A, f107, ap))
        ELSE
            tinf = ptm(1) * pt(1)
        END IF
        
        t(1) = tinf

        ! GRADIENT VARIATIONS NOT IMPORTANT BELOW ZN1(5)
        IF (alt > zn1(5)) THEN
            g0 = ptm(4) * ps(1) * (1.0 + NRL%a_sw(20) * &
                 globe7(NRL, ps, doy, sec, g_lat, g_long, lst, f107A, f107, ap))
        ELSE
            g0 = ptm(4) * ps(1)
        END IF
        
        ! C++: tlb = ptm[1] * (1.0 + NRL.a_sw[17] * globe7(NRL , pd[3], ...)) * pd[3][0];
        tlb = ptm(2) * (1.0 + NRL%a_sw(18) * &
              globe7(NRL, pd(4, :), doy, sec, g_lat, g_long, lst, f107A, f107, ap)) * pd(4, 1)
        
        s = g0 / (tinf - tlb)

        IF (alt < 300.0) THEN
            NRL%a_meso_tn1(2) = ptm(7) * ptl(1, 1) / (1.0 - NRL%a_sw(19) * glob7s(NRL, ptl(1, :), doy, g_long))
            NRL%a_meso_tn1(3) = ptm(3) * ptl(2, 1) / (1.0 - NRL%a_sw(19) * glob7s(NRL, ptl(2, :), doy, g_long))
            NRL%a_meso_tn1(4) = ptm(8) * ptl(3, 1) / (1.0 - NRL%a_sw(19) * glob7s(NRL, ptl(3, :), doy, g_long))
            NRL%a_meso_tn1(5) = ptm(5) * ptl(4, 1) / (1.0 - NRL%a_sw(19) * NRL%a_sw(21) * glob7s(NRL, ptl(4, :), doy, g_long))
            NRL%a_meso_tgn1(2) = ptm(9) * pma(9, 1) * (1.0 + NRL%a_sw(19) * NRL%a_sw(21) * &
                                 glob7s(NRL, pma(9, :), doy, g_long)) * NRL%a_meso_tn1(5) * NRL%a_meso_tn1(5) / ((ptm(5) * ptl(4, 1))**2.0)
        ELSE
            NRL%a_meso_tn1(2) = ptm(7) * ptl(1, 1)
            NRL%a_meso_tn1(3) = ptm(3) * ptl(2, 1)
            NRL%a_meso_tn1(4) = ptm(8) * ptl(3, 1)
            NRL%a_meso_tn1(5) = ptm(5) * ptl(4, 1)
            NRL%a_meso_tgn1(2) = ptm(9) * pma(9, 1) * NRL%a_meso_tn1(5) * NRL%a_meso_tn1(5) / ((ptm(5) * ptl(4, 1))**2.0)
        END IF

        ! N2 variation factor at Zlb
        g28 = NRL%a_sw(22) * globe7(NRL, pd(3, :), doy, sec, g_lat, g_long, lst, f107A, f107, ap)

        ! VARIATION OF TURBOPAUSE HEIGHT
        zhf = pdl(2, 25) * (1.0 + NRL%a_sw(6) * pdl(1, 25) * SIN(d_DGTR * g_lat) * COS(d_DR * (doy - pt(14))))
        
        ! C++: t[0] = tinf; (t(1) 已在上面设置)
        xmm = pdm(3, 5)
        z = alt

        !**** N2 DENSITY ****
        db28 = pdm(3, 1) * EXP(g28) * pd(3, 1)
        d(3) = densu(z, db28, tinf, tlb, 28.0_RL, alpha(3), t(2), ptm(6), s, mn1, zn1, NRL%a_meso_tn1, NRL%a_meso_tgn1)
        dd = d(3)
        ! Turbopause
        zh28 = pdm(3, 3) * zhf
        zhm28 = pdm(3, 4) * pdl(2, 6)
        xmd = 28.0 - xmm
        b28 = densu(zh28, db28, tinf, tlb, xmd, (alpha(3) - 1.0), t(2), ptm(6), s, mn1, zn1, NRL%a_meso_tn1, NRL%a_meso_tgn1)
        IF ((NRL%a_sw(16) /= 0.0) .AND. (z <= altl(3))) THEN
            NRL%d_dm28 = densu(z, b28, tinf, tlb, xmm, alpha(3), t(2), ptm(6), s, mn1, zn1, NRL%a_meso_tn1, NRL%a_meso_tgn1)
            d(3) = dnet(d(3), NRL%d_dm28, zhm28, xmm, 28.0_RL)
        END IF

        !**** HE DENSITY ****
        g4 = NRL%a_sw(22) * globe7(NRL, pd(1, :), doy, sec, g_lat, g_long, lst, f107A, f107, ap)
        db04 = pdm(1, 1) * EXP(g4) * pd(1, 1)
        d(1) = densu(z, db04, tinf, tlb, 4.0_RL, alpha(1), t(2), ptm(6), s, mn1, zn1, NRL%a_meso_tn1, NRL%a_meso_tgn1)
        dd = d(1)

        IF ((NRL%a_sw(16) /= 0.0) .AND. (z < altl(1))) THEN
            zh04 = pdm(1, 3)
            b04 = densu(zh04, db04, tinf, tlb, 4.0 - xmm, alpha(1) - 1.0, t(2), ptm(6), s, mn1, zn1, NRL%a_meso_tn1, NRL%a_meso_tgn1)
            NRL%d_dm04 = densu(z, b04, tinf, tlb, xmm, 0.0_RL, t(2), ptm(6), s, mn1, zn1, NRL%a_meso_tn1, NRL%a_meso_tgn1)
            zhm04 = zhm28
            d(1) = dnet(d(1), NRL%d_dm04, zhm04, xmm, 4.0_RL)
            rll = LOG(b28 * pdm(1, 2) / b04)
            zc04 = pdm(1, 5) * pdl(2, 1)
            hc04 = pdm(1, 6) * pdl(2, 2)
            d(1) = d(1) * ccor(z, rll, hc04, zc04)
        END IF

        !**** O DENSITY ****
        g16 = NRL%a_sw(22) * globe7(NRL, pd(2, :), doy, sec, g_lat, g_long, lst, f107A, f107, ap)
        db16 = pdm(2, 1) * EXP(g16) * pd(2, 1)
        d(2) = densu(z, db16, tinf, tlb, 16.0_RL, alpha(2), t(2), ptm(6), s, mn1, zn1, NRL%a_meso_tn1, NRL%a_meso_tgn1)
        dd = d(2)

        IF ((NRL%a_sw(16) /= 0.0) .AND. (z <= altl(2))) THEN
            zh16 = pdm(2, 3)
            b16 = densu(zh16, db16, tinf, tlb, 16.0 - xmm, (alpha(2) - 1.0), t(2), ptm(6), s, mn1, zn1, NRL%a_meso_tn1, NRL%a_meso_tgn1)
            NRL%d_dm16 = densu(z, b16, tinf, tlb, xmm, 0.0_RL, t(2), ptm(6), s, mn1, zn1, NRL%a_meso_tn1, NRL%a_meso_tgn1)
            zhm16 = zhm28
            d(2) = dnet(d(2), NRL%d_dm16, zhm16, xmm, 16.0_RL)
            rll = pdm(2, 2) * pdl(2, 17) * (1.0 + NRL%a_sw(2) * pdl(1, 24) * (f107A - 150.0))
            hc16 = pdm(2, 6) * pdl(2, 4)
            zc16 = pdm(2, 5) * pdl(2, 3)
            hc216 = pdm(2, 6) * pdl(2, 5)
            d(2) = d(2) * ccor2(z, rll, hc16, zc16, hc216)
            hcc16 = pdm(2, 8) * pdl(2, 14)
            zcc16 = pdm(2, 7) * pdl(2, 13)
            rc16 = pdm(2, 4) * pdl(2, 15)
            d(2) = d(2) * ccor(z, rc16, hcc16, zcc16)
        END IF

        !**** O2 DENSITY ****
        g32 = NRL%a_sw(22) * globe7(NRL, pd(5, :), doy, sec, g_lat, g_long, lst, f107A, f107, ap)
        db32 = pdm(4, 1) * EXP(g32) * pd(5, 1)
        d(4) = densu(z, db32, tinf, tlb, 32.0_RL, alpha(4), t(2), ptm(6), s, mn1, zn1, NRL%a_meso_tn1, NRL%a_meso_tgn1)
        dd = d(4)
        IF (NRL%a_sw(16) /= 0.0) THEN
            IF (z <= altl(4)) THEN
                zh32 = pdm(4, 3)
                b32 = densu(zh32, db32, tinf, tlb, 32.0 - xmm, alpha(4) - 1.0, t(2), ptm(6), s, mn1, zn1, NRL%a_meso_tn1, NRL%a_meso_tgn1)
                NRL%d_dm32 = densu(z, b32, tinf, tlb, xmm, 0.0_RL, t(2), ptm(6), s, mn1, zn1, NRL%a_meso_tn1, NRL%a_meso_tgn1)
                zhm32 = zhm28
                d(4) = dnet(d(4), NRL%d_dm32, zhm32, xmm, 32.0_RL)
                rll = LOG(b28 * pdm(4, 2) / b32)
                hc32 = pdm(4, 6) * pdl(2, 8)
                zc32 = pdm(4, 5) * pdl(2, 7)
                d(4) = d(4) * ccor(z, rll, hc32, zc32)
            END IF
            hcc32 = pdm(4, 8) * pdl(2, 23)
            hcc232 = pdm(4, 8) * pdl(1, 23)
            zcc32 = pdm(4, 7) * pdl(2, 22)
            rc32 = pdm(4, 4) * pdl(2, 24) * (1.0 + NRL%a_sw(2) * pdl(1, 24) * (f107A - 150.0))
            d(4) = d(4) * ccor2(z, rc32, hcc32, zcc32, hcc232)
        END IF

        !**** AR DENSITY ****
        g40 = NRL%a_sw(22) * globe7(NRL, pd(6, :), doy, sec, g_lat, g_long, lst, f107A, f107, ap)
        db40 = pdm(5, 1) * EXP(g40) * pd(6, 1)
        d(5) = densu(z, db40, tinf, tlb, 40.0_RL, alpha(5), t(2), ptm(6), s, mn1, zn1, NRL%a_meso_tn1, NRL%a_meso_tgn1)
        dd = d(5)
        
        IF ((NRL%a_sw(16) /= 0.0) .AND. (z <= altl(5))) THEN
            zh40 = pdm(5, 3)
            b40 = densu(zh40, db40, tinf, tlb, 40.0 - xmm, alpha(5) - 1.0, t(2), ptm(6), s, mn1, zn1, NRL%a_meso_tn1, NRL%a_meso_tgn1)
            NRL%d_dm40 = densu(z, b40, tinf, tlb, xmm, 0.0_RL, t(2), ptm(6), s, mn1, zn1, NRL%a_meso_tn1, NRL%a_meso_tgn1)
            zhm40 = zhm28
            d(5) = dnet(d(5), NRL%d_dm40, zhm40, xmm, 40.0_RL)
            rll = LOG(b28 * pdm(5, 2) / b40)
            hc40 = pdm(5, 6) * pdl(2, 10)
            zc40 = pdm(5, 5) * pdl(2, 9)
            d(5) = d(5) * ccor(z, rll, hc40, zc40)
        END IF

        !**** HYDROGEN DENSITY ****
        g1 = NRL%a_sw(22) * globe7(NRL, pd(7, :), doy, sec, g_lat, g_long, lst, f107A, f107, ap)
        db01 = pdm(6, 1) * EXP(g1) * pd(7, 1)
        d(7) = densu(z, db01, tinf, tlb, 1.0_RL, alpha(7), t(2), ptm(6), s, mn1, zn1, NRL%a_meso_tn1, NRL%a_meso_tgn1)
        dd = d(7)
        IF ((NRL%a_sw(16) /= 0.0) .AND. (z <= altl(7))) THEN
            zh01 = pdm(6, 3)
            b01 = densu(zh01, db01, tinf, tlb, 1.0 - xmm, alpha(7) - 1.0, t(2), ptm(6), s, mn1, zn1, NRL%a_meso_tn1, NRL%a_meso_tgn1)
            NRL%d_dm01 = densu(z, b01, tinf, tlb, xmm, 0.0_RL, t(2), ptm(6), s, mn1, zn1, NRL%a_meso_tn1, NRL%a_meso_tgn1)
            zhm01 = zhm28
            d(7) = dnet(d(7), NRL%d_dm01, zhm01, xmm, 1.0_RL)
            rll = LOG(b28 * pdm(6, 2) * SQRT(pdl(2, 18)**2.0) / b01)
            hc01 = pdm(6, 6) * pdl(2, 12)
            zc01 = pdm(6, 5) * pdl(2, 11)
            d(7) = d(7) * ccor(z, rll, hc01, zc01)
            hcc01 = pdm(6, 8) * pdl(2, 20)
            zcc01 = pdm(6, 7) * pdl(2, 19)
            rc01 = pdm(6, 4) * pdl(2, 21)
            d(7) = d(7) * ccor(z, rc01, hcc01, zcc01)
        END IF

        !**** ATOMIC NITROGEN DENSITY ****
        g14 = NRL%a_sw(22) * globe7(NRL, pd(8, :), doy, sec, g_lat, g_long, lst, f107A, f107, ap)
        db14 = pdm(7, 1) * EXP(g14) * pd(8, 1)
        d(8) = densu(z, db14, tinf, tlb, 14.0_RL, alpha(8), t(2), ptm(6), s, mn1, zn1, NRL%a_meso_tn1, NRL%a_meso_tgn1)
        dd = d(8)
        IF ((NRL%a_sw(16) /= 0.0) .AND. (z <= altl(8))) THEN
            zh14 = pdm(7, 3)
            b14 = densu(zh14, db14, tinf, tlb, 14.0 - xmm, alpha(8) - 1.0, t(2), ptm(6), s, mn1, zn1, NRL%a_meso_tn1, NRL%a_meso_tgn1)
            NRL%d_dm14 = densu(z, b14, tinf, tlb, xmm, 0.0_RL, t(2), ptm(6), s, mn1, zn1, NRL%a_meso_tn1, NRL%a_meso_tgn1)
            zhm14 = zhm28
            d(8) = dnet(d(8), NRL%d_dm14, zhm14, xmm, 14.0_RL)
            rll = LOG(b28 * pdm(7, 2) * SQRT(pdl(1, 3)**2.0) / b14)
            hc14 = pdm(7, 6) * pdl(1, 2)
            zc14 = pdm(7, 5) * pdl(1, 1)
            d(8) = d(8) * ccor(z, rll, hc14, zc14)
            hcc14 = pdm(7, 8) * pdl(1, 5)
            zcc14 = pdm(7, 7) * pdl(1, 4)
            rc14 = pdm(7, 4) * pdl(1, 6)
            d(8) = d(8) * ccor(z, rc14, hcc14, zcc14)
        END IF

        !**** Anomalous OXYGEN DENSITY ****
        g16h = NRL%a_sw(22) * globe7(NRL, pd(9, :), doy, sec, g_lat, g_long, lst, f107A, f107, ap)
        db16h = pdm(8, 1) * EXP(g16h) * pd(9, 1)
        tho = pdm(8, 10) * pdl(1, 7)
        dd = densu(z, db16h, tho, tho, 16.0_RL, alpha(9), t(2), ptm(6), s, mn1, zn1, NRL%a_meso_tn1, NRL%a_meso_tgn1)
        zsht = pdm(8, 6)
        zmho = pdm(8, 5)
        zsho = scalh(NRL, zmho, 16.0_RL, tho)
        d(9) = dd * EXP(-zsht / zsho * (EXP(-(z - zmho) / zsht) - 1.0))

        ! total mass density
        d(6) = 1.66E-24 * (4.0 * d(1) + 16.0 * d(2) + 28.0 * d(3) + &
                                 32.0 * d(4) + 40.0 * d(5) + d(7) + 14.0 * d(8))

        ! temperature
        z = SQRT(alt**2.0)
        dd = densu(z, 1.0_RL, tinf, tlb, 0.0_RL, 0.0_RL, t(2), ptm(6), s, mn1, zn1, NRL%a_meso_tn1, NRL%a_meso_tgn1)

        ! C++: if (NRL.a_sw[0])
        IF (NRL%a_sw(1) /= 0.0) THEN
            d(1:9) = d(1:9) * 1.0E6
            d(6) = d(6) / 1000.0
        END IF

END SUBROUTINE gts7


! ======================================================================
! C++: void gtd7(NRLMSISE &NRL, int doy, double sec, double &alt, ...)
! ======================================================================
    SUBROUTINE gtd7(NRL, doy, sec, alt, g_lat, g_long, lst, f107A, f107, ap, d, t)
        
        ! --- 参数 (Arguments) ---
        TYPE(NRLMSISE), INTENT(INOUT) :: NRL
        INTEGER(IT), INTENT(IN)           :: doy
        REAL(RL), INTENT(IN)  :: sec, alt, g_lat, g_long, lst, f107A, f107
        REAL(RL), DIMENSION(*), INTENT(IN) :: ap
        REAL(RL), DIMENSION(*), INTENT(OUT) :: d  ! C++ d[9] (0..8)
        REAL(RL), DIMENSION(*), INTENT(OUT) :: t  ! C++ t[2] (0..1)
        
        ! --- 局部变量 ---
        REAL(RL) :: xlat, xmm
        INTEGER(IT), PARAMETER :: mn3 = 5
        REAL(RL), DIMENSION(5) :: zn3 = (/ 32.5, 20.0, 15.0, 10.0, 0.0 /)
        INTEGER(IT), PARAMETER :: mn2 = 4
        REAL(RL), DIMENSION(4) :: zn2 = (/ 72.5, 55.0, 45.0, 32.5 /)
        REAL(RL) :: altt
        REAL(RL), PARAMETER :: zmix = 62.5
        REAL(RL) :: dm28m, tz, dmc, dmr, dz28
        
        ! 临时 (scratch) 数组
        REAL(RL), DIMENSION(9) :: sdens
        REAL(RL), DIMENSION(2) :: stemp
        
        INTEGER(IT) :: i
        
        ! C++ 局部变量，遮蔽了全局变量
        REAL(RL) :: local_d_gsurf = 0.0
        REAL(RL) :: local_d_re = 0.0
        
        ! --- 外部函数/子程序 ---
        !REAL(RL) :: glob7s, densm
        
        ! --- 函数体 ---
        xlat = g_lat
        IF (NRL%a_sw(3) == 0.0) xlat = 45.0
        
        CALL glatf(xlat, local_d_gsurf, local_d_re)
        
        xmm = pdm(3, 5)

        IF (alt > zn2(1)) THEN
            altt = alt
        ELSE
            altt = zn2(1)
        END IF
        
        ! 首先调用 gts7 (热层)
        CALL gts7(NRL, doy, sec, altt, g_lat, g_long, lst, f107A, f107, ap, sdens, stemp)
        
        altt = alt
        
        IF (NRL%a_sw(1) /= 0.0) THEN
            dm28m = NRL%d_dm28 * 1.0E6
        ELSE
            dm28m = NRL%d_dm28
        END IF
        

        t(1) = stemp(1)
        t(2) = stemp(2)
        
        IF (alt >= zn2(1)) THEN
            d(1:9) = sdens(1:9)
            RETURN
        END IF
        
        ! --- 高度低于 zn2(1) [72.5km]，需要计算中层 ---

        ! LOWER MESOSPHERE/UPPER STRATOSPHERE (between zn3[0] and zn2[0])
        NRL%a_meso_tgn2(1) = NRL%a_meso_tgn1(2)
        NRL%a_meso_tn2(1) = NRL%a_meso_tn1(5)
        
        NRL%a_meso_tn2(2) = pma(1, 1) * pavgm(1) / (1.0 - NRL%a_sw(21) * glob7s(NRL, pma(1, :), doy, g_long))
        NRL%a_meso_tn2(3) = pma(2, 1) * pavgm(2) / (1.0 - NRL%a_sw(21) * glob7s(NRL, pma(2, :), doy, g_long))
        NRL%a_meso_tn2(4) = pma(3, 1) * pavgm(3) / (1.0 - NRL%a_sw(21) * NRL%a_sw(23) * glob7s(NRL, pma(3, :), doy, g_long))
        
        NRL%a_meso_tgn2(2) = pavgm(9) * pma(10, 1) * (1.0 + NRL%a_sw(21) * NRL%a_sw(23) * &
                             glob7s(NRL, pma(10, :), doy, g_long)) * NRL%a_meso_tn2(4) * NRL%a_meso_tn2(4) / ((pma(3, 1) * pavgm(3))**2.0)
        

        NRL%a_meso_tn3(1) = NRL%a_meso_tn2(4)

        ! C++: if (alt <= zn3[0])
        IF (alt <= zn3(1)) THEN
            ! LOWER STRATOSPHERE AND TROPOSPHERE (below zn3[0])
            NRL%a_meso_tgn3(1) = NRL%a_meso_tgn2(2)
            NRL%a_meso_tn3(2) = pma(4, 1) * pavgm(4) / (1.0 - NRL%a_sw(23) * glob7s(NRL, pma(4, :), doy, g_long))
            NRL%a_meso_tn3(3) = pma(5, 1) * pavgm(5) / (1.0 - NRL%a_sw(23) * glob7s(NRL, pma(5, :), doy, g_long))
            NRL%a_meso_tn3(4) = pma(6, 1) * pavgm(6) / (1.0 - NRL%a_sw(23) * glob7s(NRL, pma(6, :), doy, g_long))
            NRL%a_meso_tn3(5) = pma(7, 1) * pavgm(7) / (1.0 - NRL%a_sw(23) * glob7s(NRL, pma(7, :), doy, g_long))
            
            NRL%a_meso_tgn3(2) = pma(8, 1) * pavgm(8) * (1.0 + NRL%a_sw(23) * &
                                 glob7s(NRL, pma(8, :), doy, g_long)) * NRL%a_meso_tn3(5) * NRL%a_meso_tn3(5) / ((pma(7, 1) * pavgm(7))**2.0)
        END IF

        ! LINEAR TRANSITION TO FULL MIXING BELOW zn2[0]
        dmc = 0.0
        IF (alt > zmix) THEN
            dmc = 1.0 - (zn2(1) - alt) / (zn2(1) - zmix)
        END IF
        dz28 = sdens(3)

        !**** N2 density ****
        dmr = sdens(3) / dm28m - 1.0
        d(3) = densm(NRL, alt, dm28m, xmm, tz, mn3, zn3, NRL%a_meso_tn3, NRL%a_meso_tgn3, mn2, zn2, NRL%a_meso_tn2, NRL%a_meso_tgn2)
        d(3) = d(3) * (1.0 + dmr * dmc)

        !**** HE density ****
        dmr = sdens(1) / (dz28 * pdm(1, 2)) - 1.0
        d(1) = d(3) * pdm(1, 2) * (1.0 + dmr * dmc)

        !**** O density ****
        d(2) = 0.0
        d(9) = 0.0

        !**** O2 density ****
        dmr = sdens(4) / (dz28 * pdm(4, 2)) - 1.0
        d(4) = d(3) * pdm(4, 2) * (1.0 + dmr * dmc)

        !**** AR density ***
        dmr = sdens(5) / (dz28 * pdm(5, 2)) - 1.0
        d(5) = d(3) * pdm(5, 2) * (1.0 + dmr * dmc)

        !**** Hydrogen density ****
        d(7) = 0.0

        !**** Atomic nitrogen density ****
        d(8) = 0.0

        !**** Total mass density */
        d(6) = 1.66E-24 * (4.0 * d(1) + 16.0 * d(2) + 28.0 * d(3) &
                               + 32.0 * d(4) + 40.0 * d(5) + d(7) + 14.0 * d(8))

        ! C++: if (NRL.a_sw[0])
        IF (NRL%a_sw(1) /= 0.0) THEN
            d(6) = d(6) / 1000.0
        END IF

        !**** temperature at altitude ****
        NRL%d_dd = densm(NRL, alt, 1.0_RL, 0.0_RL, tz, mn3, zn3, NRL%a_meso_tn3, NRL%a_meso_tgn3, mn2, zn2, NRL%a_meso_tn2, NRL%a_meso_tgn2)
        t(2) = tz

    END SUBROUTINE gtd7
    
! ======================================================================
! C++: void gtd7d(NRLMSISE &NRL, int doy, double sec, double &alt, ...)
!      (GTD7 WITH ANOMALOUS OXYGEN INCLUDED IN TOTAL DENSITY)
! ======================================================================   
    
    SUBROUTINE gtd7d(NRL, doy, sec, alt, g_lat, g_long, lst, f107A, f107, ap, d, t)
        
        ! --- 参数 (Arguments) ---
        TYPE(NRLMSISE), INTENT(INOUT) :: NRL
        INTEGER(IT), INTENT(IN)           :: doy
        REAL(RL), INTENT(IN)  :: sec, alt, g_lat, g_long, lst, f107A, f107
        REAL(RL), DIMENSION(*), INTENT(IN) :: ap
        REAL(RL), DIMENSION(*), INTENT(OUT) :: d  ! C++ d[9] (0..8) -> F90 d(9) (1..9)
        REAL(RL), DIMENSION(*), INTENT(OUT) :: t  ! C++ t[2] (0..1) -> F90 t(2) (1..2)
        
        ! --- 函数体 ---
        
        ! 1. 首先，像 gtd7 一样计算所有值
        CALL gtd7(NRL, doy, sec, alt, g_lat, g_long, lst, f107A, f107, ap, d, t)
        
        ! 2. 覆盖总质量密度 d(6) (C++ d[5])

        d(6) = 1.66E-24 * (4.0 * d(1) + 16.0 * d(2) + 28.0 * d(3) + &
                                 32.0 * d(4) + 40.0 * d(5) + d(7) + 14.0 * d(8) + &
                                 16.0 * d(9)) ! 包含异常氧 d(9)
        
        ! 3. 应用公制单位转换
        IF (NRL%a_sw(1) /= 0.0) THEN
            d(6) = d(6) / 1000.0
        END IF
        
    END SUBROUTINE gtd7d


    ! C interoperable wrapper used by GROOPS. Inputs follow the conventional
    ! NRLMSISE-00 units: altitude [km], angles [deg], time [s], flux [sfu].
    SUBROUTINE msise00CalcWrapper(doy, sec, alt, lat, lon, lst, f107A, f107, ap, &
                                  density, temperature, exosphericTemperature) &
                                  BIND(C, name="msise00CalcWrapper")
        USE, INTRINSIC :: ISO_C_BINDING
        IMPLICIT NONE
        INTEGER(C_INT), INTENT(IN) :: doy
        REAL(C_DOUBLE), INTENT(IN) :: sec, alt, lat, lon, lst, f107A, f107
        REAL(C_DOUBLE), INTENT(IN) :: ap(7)
        REAL(C_DOUBLE), INTENT(OUT) :: density, temperature, exosphericTemperature
        TYPE(NRLMSISE) :: NRL
        REAL(RL) :: d(9), t(2)

        CALL gts7(NRL, doy, sec, alt, lat, lon, lst, f107A, f107, ap, d, t)
        density = d(6)
        exosphericTemperature = t(1)
        temperature = t(2)
    END SUBROUTINE msise00CalcWrapper
                
        
END MODULE NRLMSISE00_FUNCTIONS
    
