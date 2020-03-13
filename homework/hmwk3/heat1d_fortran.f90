subroutine init_q_fort(N,mbc,a,b,init,q)
    implicit none

    external init 

    integer N, mbc
    double precision a,b
    double precision q(-mbc:N+mbc)
    double precision init

    integer i
    double precision x, dx

    dx = (b-a)/N

    do i = -mbc, N+mbc
        x = a + i*dx
        q(i) = init(x)
    end do

end subroutine init_q_fort

subroutine update_q_fort(N,mbc,dt, dx, q, qp)
    implicit none

    integer N, mbc
    double precision dt, dx
    double precision  q(-mbc:N+mbc)
    double precision qp(-mbc:N+mbc)

    integer i, m
    double precision dx2

    dx2 = dx*dx

    !! Set boundary conditions
    do m = 1,mbc
        q(-m) = q(m)
        q(N+m) = q(N-m)
    end do

    !! Update qp
    do i = 0,N
        qp(i) = q(i) + dt*(q(i-1) - 2*q(i) + q(i+1))/dx2
    end do

end subroutine update_q_fort