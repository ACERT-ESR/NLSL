	integer i
        character*20 c
        data c/'abcdefghijklmnopqrst'/
 10	read(*,*)i
        write(*,*)'c(i)i = ','test'//c(i:i)
        go to 10

        stop
        end
