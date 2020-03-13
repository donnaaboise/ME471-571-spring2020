subroutine init_q_fort(N,M,mbc,ax,bx,ay,by,init,q)
    implicit none

    external init 

    integer N, M, mbc
    double precision ax,bx,ay, by
    double precision q(-mbc:N+mbc, -mbc:M+mbc)
    double precision init

    integer i, j
    double precision x, dx, y, dy

    dx = (bx-ax)/N
    dy = (by-ay)/M

    do i = -mbc, N+mbc
        x = ax + i*dx
        do j = -mbc,M+mbc
            y = ay + j*dy
            q(i,j) = init(x,y)
        end do
    end do

end subroutine init_q_fort

subroutine update_q_fort(N,M,mbc,dt, dx, dy, q, qp)
    implicit none

    integer N, M, mbc
    double precision dt, dx, dy
    double precision  q(-mbc:N+mbc,-mbc:M+mbc)
    double precision qp(-mbc:N+mbc,-mbc:M+mbc)

    integer i, j
    double precision dx2, dy2, qxx, qyy

    dx2 = dx*dx
    dy2 = dy*dy

    !! Set boundary conditions
    do m = 1,mbc
        do j = 0,M
            q(-m,j) = q(m,j)
            q(N+m,j) = q(N-m,j)
        end do
        do i = 0,N
            q(i,-m) = q(i,m)
            q(i,M+m) = q(i,M-m)
        end do
    end do


    !! Update qp
    do j = 0,M
        do i = 0,N        
            qxx = (q(i-1,j) - 2*q(i,j) + q(i+1,j))/dx2
            qyy = (q(i,j-1) - 2*q(i,j) + q(i,j+1))/dy2
            qp(i,j) = q(i,j) + dt*(qxx + qyy)
        end do
    end do

end subroutine update_q_fort