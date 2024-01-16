
module My_func
    implicit none 

    contains

    !@cuf attributes(host, device) & 
    ! real(8) pure function  linear(x1, t1, x2, t2, x3, t3, y)
    !     ! Главное значение с параметрами 2
    !     ! Строим линии между 1 и 2,  2 и 3, потом находим минмодом значение в y
    !     implicit none
    !     real(8), intent(in) :: x1, x2, x3, y, t1, t2, t3
    !     real(8) :: d

    !     d = minmod((t1 - t2) / (x1 - x2), (t2 - t3) / (x2 - x3))
    !     linear =  (d * (y - x2) + t2)
    !     return
        
    ! end function linear

    !@cuf attributes(host, device) & 
    real(8) pure function  linear1(x1, t1, x2, t2, y)
        implicit none
        real(8), intent(in) :: x1, x2, y, t1, t2
        real(8) :: d

        d = (t1 - t2) / (x1 - x2)
        linear1 =  (d * (y - x2) + t2)
        return
        
    end function linear1

end module My_func