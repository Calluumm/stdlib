program test_distribution_gamma
    use stdlib_kinds, only : sp, dp, xdp
    use stdlib_error, only : check
    use stdlib_random, only : random_seed
    use stdlib_stats_distribution_gamma, only : rgamma => rvs_gamma,           &
                              gamma_pdf => pdf_gamma, gamma_cdf => cdf_gamma

    implicit none
    real(sp), parameter :: sptol = 1000 * epsilon(1.0_sp)
    real(dp), parameter :: dptol = 1000 * epsilon(1.0_dp)
    logical ::  warn = .true.
    integer :: put, get

    put = 12345678
    call random_seed(put, get)

    call test_gamma_random_generator

    call test_gamma_rvs_rsp
    call test_gamma_rvs_rdp
    call test_gamma_rvs_csp
    call test_gamma_rvs_cdp

    call test_gamma_pdf_rsp
    call test_gamma_pdf_rdp
    call test_gamma_pdf_csp
    call test_gamma_pdf_cdp

    call test_gamma_cdf_rsp
    call test_gamma_cdf_rdp
    call test_gamma_cdf_csp
    call test_gamma_cdf_cdp

contains

    subroutine test_gamma_random_generator
        integer, parameter :: num = 10000000, array_size = 1000
        integer :: i, j, freq(0:array_size-1)
        real(dp) :: chisq, expct

        print *, ""
        print *, "Test gamma random generator with chi-squared"

        freq = 0
        do i = 1, num
            j = min(array_size - 1, int(array_size * gamma_cdf(rgamma(2.0_dp, 1.5_dp), 2.0_dp, 1.5_dp)))
            freq(j) = freq(j) + 1
        end do

        chisq = 0.0_dp
        expct = num / array_size

        do i = 0, array_size - 1
           chisq = chisq + (freq(i) - expct) ** 2 / expct
        end do

        write(*,*) "The critical values for chi-squared with 1000 d. of f. is" &
            //" 1143.92"
        write(*,*) "Chi-squared for gamma random generator is : ", chisq
        call check((chisq < 1143.9), &
            msg="gamma randomness failed chi-squared test", warn=warn)

    end subroutine test_gamma_random_generator

    subroutine test_gamma_rvs_rsp
        integer, parameter :: k = 5, n = 10
        real(sp) :: res(n), gshape, rate
        integer :: i
        integer :: seed, get
        real(sp) :: ans(n) = [0.85758907497718884_sp,                        &
                            1.0206623865526090_sp,                         &
                            0.99753931024198650_sp,                        &
                            0.97653359790345839_sp,                        &
                            0.41853482638322043_sp,                        &
                            2.2012288073086310_sp,                         &
                            2.0639542613306592_sp,                         &
                            3.1794669730880192_sp,                         &
                            1.9329744662223280_sp,                         &
                            1.0257959670932111_sp]

        print *, "Test gamma_distribution_rvs_rsp"
        seed = 639741825
        call random_seed(seed, get)

        gshape = 2.0_sp; rate = 1.0_sp

        do i = 1, k
            res(i) = rgamma(gshape, rate)
        end do

        res(k + 1 : n) = rgamma(gshape, rate, k)

        do i = 1, n
            call check(abs(res(i) - ans(i)) < sptol,                      &
                       msg="gamma_distribution_rvs_rsp failed", &
                       warn=warn)
        end do
    end subroutine test_gamma_rvs_rsp

    subroutine test_gamma_rvs_rdp
        integer, parameter :: k = 5, n = 10
        real(dp) :: res(n), gshape, rate
        integer :: i
        integer :: seed, get
        real(dp) :: ans(n) = [0.85758907497718884_dp,                        &
                            1.0206623865526090_dp,                         &
                            0.99753931024198650_dp,                        &
                            0.97653359790345839_dp,                        &
                            0.41853482638322043_dp,                        &
                            2.2012288073086310_dp,                         &
                            2.0639542613306592_dp,                         &
                            3.1794669730880192_dp,                         &
                            1.9329744662223280_dp,                         &
                            1.0257959670932111_dp]

        print *, "Test gamma_distribution_rvs_rdp"
        seed = 639741825
        call random_seed(seed, get)

        gshape = 2.0_dp; rate = 1.0_dp

        do i = 1, k
            res(i) = rgamma(gshape, rate)
        end do

        res(k + 1 : n) = rgamma(gshape, rate, k)

        do i = 1, n
            call check(abs(res(i) - ans(i)) < dptol,                      &
                       msg="gamma_distribution_rvs_rdp failed", &
                       warn=warn)
        end do
    end subroutine test_gamma_rvs_rdp

    subroutine test_gamma_rvs_csp
        integer, parameter :: k = 5, n = 10
        complex(sp) :: res(n), gshape, rate
        integer :: i
        integer :: seed, get
        complex(sp) :: ans(n) =                                                     &
                     [(1.0719863437214860_sp, 0.46775532101393819_sp), &
                     (0.42382516926807201_sp, 0.96340496644915230_sp), &
                      (2.7515360091357888_sp, 0.14837198853150388_sp), &
                      (1.4536367104245524_sp, 0.56852736336951559_sp), &
                (0.34559143458416125_sp, 4.96217685362488267E-002_sp), &
                      (1.9657884897696516_sp,  3.1124314799641013_sp), &
                 (3.4155160623540453_sp, 5.04948933894018709E-002_sp), &
                     (0.94594398345216302_sp, 0.45691588305890624_sp), &
                      (1.1493158751025965_sp, 0.12944763723941669_sp), &
                      (2.9691469633592282_sp, 1.1617408197125874_sp)]

        print *, "Test gamma_distribution_rvs_csp"
        seed = 639741825
        call random_seed(seed, get)

        gshape = (2.0_sp, 0.7_sp); rate = (0.8_sp, 1.2_sp)

        do i = 1, k
            res(i) = rgamma(gshape, rate)
        end do

        res(k + 1 : n) = rgamma(gshape, rate, k)

        do i = 1, n
            call check(abs(res(i) - ans(i)) < sptol,                      &
                       msg="gamma_distribution_rvs_csp failed", &
                       warn=warn)
        end do
    end subroutine test_gamma_rvs_csp

    subroutine test_gamma_rvs_cdp
        integer, parameter :: k = 5, n = 10
        complex(dp) :: res(n), gshape, rate
        integer :: i
        integer :: seed, get
        complex(dp) :: ans(n) =                                                     &
                     [(1.0719863437214860_dp, 0.46775532101393819_dp), &
                     (0.42382516926807201_dp, 0.96340496644915230_dp), &
                      (2.7515360091357888_dp, 0.14837198853150388_dp), &
                      (1.4536367104245524_dp, 0.56852736336951559_dp), &
                (0.34559143458416125_dp, 4.96217685362488267E-002_dp), &
                      (1.9657884897696516_dp,  3.1124314799641013_dp), &
                 (3.4155160623540453_dp, 5.04948933894018709E-002_dp), &
                     (0.94594398345216302_dp, 0.45691588305890624_dp), &
                      (1.1493158751025965_dp, 0.12944763723941669_dp), &
                      (2.9691469633592282_dp, 1.1617408197125874_dp)]

        print *, "Test gamma_distribution_rvs_cdp"
        seed = 639741825
        call random_seed(seed, get)

        gshape = (2.0_dp, 0.7_dp); rate = (0.8_dp, 1.2_dp)

        do i = 1, k
            res(i) = rgamma(gshape, rate)
        end do

        res(k + 1 : n) = rgamma(gshape, rate, k)

        do i = 1, n
            call check(abs(res(i) - ans(i)) < dptol,                      &
                       msg="gamma_distribution_rvs_cdp failed", &
                       warn=warn)
        end do
    end subroutine test_gamma_rvs_cdp


    subroutine test_gamma_pdf_rsp
        real(sp) :: x1, x2(3,4), gshape, rate
        integer :: i
        integer :: seed, get
        real(sp) :: res(15)
        real(sp), parameter :: ans(15) =                                   &
                                   [3.4495412572168718E-002_sp,            &
                                    3.4495412572168718E-002_sp,            &
                                    3.4495412572168718E-002_sp,            &
                                    0.29116634347089576_sp,                &
                                    0.28338290850731412_sp,                &
                                    0.27922270935613586_sp,                &
                                    0.36440665523348270_sp,                &
                                    0.24379209619143699_sp,                &
                                    6.3815638087140858E-002_sp,            &
                                    0.25844600948718588_sp,                &
                                    0.17268118913523497_sp,                &
                                    0.31181223194308200_sp,                &
                                    0.24027095040543087_sp,                &
                                    0.36765502365831570_sp,                &
                                    9.9011714088769673E-002_sp]

        print *, "Test gamma_distribution_pdf_rsp"
        seed = 345987126
        call random_seed(seed, get)
        gshape = 2.0_sp; rate = 1.0_sp

        x1 = rgamma(gshape, rate)
        x2 = reshape(rgamma(gshape, rate, 12), [3,4])

        res(1:3) = gamma_pdf(x1, gshape, rate)
        res(4:15) = reshape(gamma_pdf(x2, gshape, rate), [12])

        do i = 1, 15
            call check(abs(res(i) - ans(i)) < sptol,                      &
                       msg="gamma_distribution_pdf_rsp failed", &
                       warn=warn)
        end do
    end subroutine test_gamma_pdf_rsp

    subroutine test_gamma_pdf_rdp
        real(dp) :: x1, x2(3,4), gshape, rate
        integer :: i
        integer :: seed, get
        real(dp) :: res(15)
        real(dp), parameter :: ans(15) =                                   &
                                   [3.4495412572168718E-002_dp,            &
                                    3.4495412572168718E-002_dp,            &
                                    3.4495412572168718E-002_dp,            &
                                    0.29116634347089576_dp,                &
                                    0.28338290850731412_dp,                &
                                    0.27922270935613586_dp,                &
                                    0.36440665523348270_dp,                &
                                    0.24379209619143699_dp,                &
                                    6.3815638087140858E-002_dp,            &
                                    0.25844600948718588_dp,                &
                                    0.17268118913523497_dp,                &
                                    0.31181223194308200_dp,                &
                                    0.24027095040543087_dp,                &
                                    0.36765502365831570_dp,                &
                                    9.9011714088769673E-002_dp]

        print *, "Test gamma_distribution_pdf_rdp"
        seed = 345987126
        call random_seed(seed, get)
        gshape = 2.0_dp; rate = 1.0_dp

        x1 = rgamma(gshape, rate)
        x2 = reshape(rgamma(gshape, rate, 12), [3,4])

        res(1:3) = gamma_pdf(x1, gshape, rate)
        res(4:15) = reshape(gamma_pdf(x2, gshape, rate), [12])

        do i = 1, 15
            call check(abs(res(i) - ans(i)) < dptol,                      &
                       msg="gamma_distribution_pdf_rdp failed", &
                       warn=warn)
        end do
    end subroutine test_gamma_pdf_rdp

    subroutine test_gamma_pdf_csp
        complex(sp) :: x1, x2(3,4), gshape, rate
        integer :: i
        integer :: seed, get
        real(sp) :: res(15)
        real(sp), parameter :: ans(15) =                                   &
                                   [0.11554282574059289_sp,                &
                                    0.11554282574059289_sp,                &
                                    0.11554282574059289_sp,                &
                                    9.2682318951901529E-002_sp,            &
                                    0.40166849087286088_sp,                &
                                    0.37468980496232701_sp,                &
                                    0.14712363446345342_sp,                &
                                    0.22561628567985184_sp,                &
                                    0.12765403024301181_sp,                &
                                    3.9182498867847360E-002_sp,            &
                                    2.5873533461032859E-003_sp,            &
                                    0.10105832622792968_sp,                &
                                    0.24044091896609490_sp,                &
                                    4.9885356046115948E-003_sp,            &
                                    0.11085827028639164_sp]

        print *, "Test gamma_distribution_pdf_csp"
        seed = 345987126
        call random_seed(seed, get)
        gshape = (2.0_sp, 0.7_sp); rate = (0.8_sp, 1.2_sp)

        x1 = rgamma(gshape, rate)
        x2 = reshape(rgamma(gshape, rate, 12), [3,4])

        res(1:3) = gamma_pdf(x1, gshape, rate)
        res(4:15) = reshape(gamma_pdf(x2, gshape, rate), [12])

        do i = 1, 15
            call check(abs(res(i) - ans(i)) < sptol,                      &
                       msg="gamma_distribution_pdf_csp failed", &
                       warn=warn)
        end do
    end subroutine test_gamma_pdf_csp

    subroutine test_gamma_pdf_cdp
        complex(dp) :: x1, x2(3,4), gshape, rate
        integer :: i
        integer :: seed, get
        real(dp) :: res(15)
        real(dp), parameter :: ans(15) =                                   &
                                   [0.11554282574059289_dp,                &
                                    0.11554282574059289_dp,                &
                                    0.11554282574059289_dp,                &
                                    9.2682318951901529E-002_dp,            &
                                    0.40166849087286088_dp,                &
                                    0.37468980496232701_dp,                &
                                    0.14712363446345342_dp,                &
                                    0.22561628567985184_dp,                &
                                    0.12765403024301181_dp,                &
                                    3.9182498867847360E-002_dp,            &
                                    2.5873533461032859E-003_dp,            &
                                    0.10105832622792968_dp,                &
                                    0.24044091896609490_dp,                &
                                    4.9885356046115948E-003_dp,            &
                                    0.11085827028639164_dp]

        print *, "Test gamma_distribution_pdf_cdp"
        seed = 345987126
        call random_seed(seed, get)
        gshape = (2.0_dp, 0.7_dp); rate = (0.8_dp, 1.2_dp)

        x1 = rgamma(gshape, rate)
        x2 = reshape(rgamma(gshape, rate, 12), [3,4])

        res(1:3) = gamma_pdf(x1, gshape, rate)
        res(4:15) = reshape(gamma_pdf(x2, gshape, rate), [12])

        do i = 1, 15
            call check(abs(res(i) - ans(i)) < dptol,                      &
                       msg="gamma_distribution_pdf_cdp failed", &
                       warn=warn)
        end do
    end subroutine test_gamma_pdf_cdp


    subroutine test_gamma_cdf_rsp
        real(sp) :: x1, x2(3,4), gshape, rate
        integer :: i
        integer :: seed, get
        real(sp) :: res(15)
        real(sp), parameter :: ans(15) =                                   &
                                   [0.14616451332312537_sp,                &
                                    0.14616451332312537_sp,                &
                                    0.14616451332312537_sp,                &
                                    0.65659267488407371_sp,                &
                                    0.60964815134398197_sp,                &
                                    0.59750757512835858_sp,                &
                                    0.69018498665337554_sp,                &
                                    0.52070621772130649_sp,                &
                                    0.23208261199896811_sp,                &
                                    0.57395848656344925_sp,                &
                                    0.50556401081889642_sp,                &
                                    0.60050726865093853_sp,                &
                                    0.52521835876447772_sp,                &
                                    0.74358834194807097_sp,                &
                                    0.35588223028193658_sp]

        print *, "Test gamma_distribution_cdf_rsp"
        seed = 345987126
        call random_seed(seed, get)
        gshape = 2.0_sp; rate = 1.0_sp

        x1 = rgamma(gshape, rate)
        x2 = reshape(rgamma(gshape, rate, 12), [3,4])

        res(1:3) = gamma_cdf(x1, gshape, rate)
        res(4:15) = reshape(gamma_cdf(x2, gshape, rate), [12])

        do i = 1, 15
            call check(abs(res(i) - ans(i)) < sptol,                      &
                       msg="gamma_distribution_cdf_rsp failed", &
                       warn=warn)
        end do
    end subroutine test_gamma_cdf_rsp

    subroutine test_gamma_cdf_rdp
        real(dp) :: x1, x2(3,4), gshape, rate
        integer :: i
        integer :: seed, get
        real(dp) :: res(15)
        real(dp), parameter :: ans(15) =                                   &
                                   [0.14616451332312537_dp,                &
                                    0.14616451332312537_dp,                &
                                    0.14616451332312537_dp,                &
                                    0.65659267488407371_dp,                &
                                    0.60964815134398197_dp,                &
                                    0.59750757512835858_dp,                &
                                    0.69018498665337554_dp,                &
                                    0.52070621772130649_dp,                &
                                    0.23208261199896811_dp,                &
                                    0.57395848656344925_dp,                &
                                    0.50556401081889642_dp,                &
                                    0.60050726865093853_dp,                &
                                    0.52521835876447772_dp,                &
                                    0.74358834194807097_dp,                &
                                    0.35588223028193658_dp]

        print *, "Test gamma_distribution_cdf_rdp"
        seed = 345987126
        call random_seed(seed, get)
        gshape = 2.0_dp; rate = 1.0_dp

        x1 = rgamma(gshape, rate)
        x2 = reshape(rgamma(gshape, rate, 12), [3,4])

        res(1:3) = gamma_cdf(x1, gshape, rate)
        res(4:15) = reshape(gamma_cdf(x2, gshape, rate), [12])

        do i = 1, 15
            call check(abs(res(i) - ans(i)) < dptol,                      &
                       msg="gamma_distribution_cdf_rdp failed", &
                       warn=warn)
        end do
    end subroutine test_gamma_cdf_rdp

    subroutine test_gamma_cdf_csp
        complex(sp) :: x1, x2(3,4), gshape, rate
        integer :: i
        integer :: seed, get
        real(sp) :: res(15)
        real(sp), parameter :: ans(15) =                                   &
                                   [0.70448339487949867_sp,                &
                                    0.70448339487949867_sp,                &
                                    0.70448339487949867_sp,                &
                                    0.58326189602977894_sp,                &
                                    0.64528366310108972_sp,                &
                                    0.64940551894231044_sp,                &
                                    0.90124666284311420_sp,                &
                                    0.84345913462074950_sp,                &
                                    0.85089208702762014_sp,                &
                                    0.77705700206310725_sp,                &
                                    0.80494892292159898_sp,                &
                                    0.87781094911086321_sp,                &
                                    0.87456078103386559_sp,                &
                                    0.88580959491539949_sp,                &
                                    0.94240247435575605_sp]

        print *, "Test gamma_distribution_cdf_csp"
        seed = 345987126
        call random_seed(seed, get)
        gshape = (2.0_sp, 0.7_sp); rate = (0.8_sp, 1.2_sp)

        x1 = rgamma(gshape, rate)
        x2 = reshape(rgamma(gshape, rate, 12), [3,4])

        res(1:3) = gamma_cdf(x1, gshape, rate)
        res(4:15) = reshape(gamma_cdf(x2, gshape, rate), [12])

        do i = 1, 15
            call check(abs(res(i) - ans(i)) < sptol,                      &
                       msg="gamma_distribution_cdf_csp failed", &
                       warn=warn)
        end do
    end subroutine test_gamma_cdf_csp

    subroutine test_gamma_cdf_cdp
        complex(dp) :: x1, x2(3,4), gshape, rate
        integer :: i
        integer :: seed, get
        real(dp) :: res(15)
        real(dp), parameter :: ans(15) =                                   &
                                   [0.70448339487949867_dp,                &
                                    0.70448339487949867_dp,                &
                                    0.70448339487949867_dp,                &
                                    0.58326189602977894_dp,                &
                                    0.64528366310108972_dp,                &
                                    0.64940551894231044_dp,                &
                                    0.90124666284311420_dp,                &
                                    0.84345913462074950_dp,                &
                                    0.85089208702762014_dp,                &
                                    0.77705700206310725_dp,                &
                                    0.80494892292159898_dp,                &
                                    0.87781094911086321_dp,                &
                                    0.87456078103386559_dp,                &
                                    0.88580959491539949_dp,                &
                                    0.94240247435575605_dp]

        print *, "Test gamma_distribution_cdf_cdp"
        seed = 345987126
        call random_seed(seed, get)
        gshape = (2.0_dp, 0.7_dp); rate = (0.8_dp, 1.2_dp)

        x1 = rgamma(gshape, rate)
        x2 = reshape(rgamma(gshape, rate, 12), [3,4])

        res(1:3) = gamma_cdf(x1, gshape, rate)
        res(4:15) = reshape(gamma_cdf(x2, gshape, rate), [12])

        do i = 1, 15
            call check(abs(res(i) - ans(i)) < dptol,                      &
                       msg="gamma_distribution_cdf_cdp failed", &
                       warn=warn)
        end do
    end subroutine test_gamma_cdf_cdp


end program test_distribution_gamma
