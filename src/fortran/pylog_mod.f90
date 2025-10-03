      module pylog_mod
      use iso_c_binding
      implicit none
      private
      public pylog_open, pylog_emit, pylog_close, log_enabled, log_buffer,
     &       ensure_log_buffer, flush_log_buffer

      logical, save :: log_enabled = .false.
      character(len=:), allocatable, save :: log_buffer(:)

      character(len=1), parameter :: sentinel = '~'

      interface
         integer(c_int) function nlsl_pylog_open(filename) bind(C,
     &        name="nlsl_pylog_open")
           import :: c_char, c_int
           character(kind=c_char), dimension(*), intent(in) :: filename
         end function nlsl_pylog_open

         integer(c_int) function nlsl_pylog_debug(message) bind(C,
     &        name="nlsl_pylog_debug")
           import :: c_char, c_int
           character(kind=c_char), dimension(*), intent(in) :: message
         end function nlsl_pylog_debug

         integer(c_int) function nlsl_pylog_close() bind(C,
     &        name="nlsl_pylog_close")
           import :: c_int
         end function nlsl_pylog_close
      end interface

      contains

      subroutine ensure_log_buffer(buffer, nlines, linelen)
        character(len=:), allocatable, intent(inout) :: buffer(:)
        integer, intent(in), optional :: nlines, linelen
        integer :: lines, length
        character(len=:), allocatable :: tmp(:)
        integer :: i, copy_lines, old_len, new_len

        lines = 8
        if (present(nlines)) lines = max(lines, nlines)
        length = 256
        if (present(linelen)) length = max(length, linelen)

        if (.not. allocated(buffer)) then
           allocate(character(len=length) :: buffer(lines))
        else if (size(buffer) < lines .or. len(buffer(1)) < length) then
           new_len = max(length, len(buffer(1)))
           allocate(character(len=new_len) :: tmp(max(lines,
     &        size(buffer))))
           tmp = sentinel
           copy_lines = min(size(buffer), size(tmp))
           old_len = min(len(buffer(1)), len(tmp(1)))
           do i = 1, copy_lines
              tmp(i)(:old_len) = buffer(i)(:old_len)
           end do
           call move_alloc(tmp, buffer)
        end if
        buffer = sentinel
      end subroutine ensure_log_buffer

      pure function to_c_string(text) result(buf)
        character(*), intent(in) :: text
        character(kind=c_char), allocatable :: buf(:)
        integer :: n, i, code

        n = len_trim(text)
        if (n < 0) n = 0
        allocate(buf(n + 1))
        do i = 1, n
           code = iachar(text(i:i))
           if (code < 0) code = 0
           buf(i) = char(code, kind=c_char)
        end do
        buf(n + 1) = c_null_char
      end function to_c_string

      logical function has_log_suffix(name)
        character(*), intent(in) :: name
        integer :: n
        character(len=4) :: suffix

        n = len_trim(name)
        if (n < 4) then
           has_log_suffix = .false.
           return
        end if
        suffix = name(n-3:n)
        suffix = uppercase(suffix)
        has_log_suffix = (suffix == '.LOG')
      end function has_log_suffix

      pure function uppercase(text) result(out)
        character(*), intent(in) :: text
        character(len=len(text)) :: out
        integer :: i

        do i = 1, len(text)
           select case (text(i:i))
           case('a':'z')
              out(i:i) = char(iachar(text(i:i)) - 32)
           case default
              out(i:i) = text(i:i)
           end select
        end do
      end function uppercase

      subroutine pylog_open(file_id)
        character(*), intent(in) :: file_id
        character(:), allocatable :: target
        character(kind=c_char), allocatable :: c_filename(:)
        integer(c_int) :: status

        call pylog_close()

        if (len_trim(file_id) == 0) then
           log_enabled = .false.
           return
        end if

        target = adjustl(trim(file_id))
        if (.not. has_log_suffix(target)) then
           target = target // '.log'
        end if

        c_filename = to_c_string(target)
        status = nlsl_pylog_open(c_filename)
        if (status == 0) then
           log_enabled = .true.
        else
           log_enabled = .false.
        end if
      end subroutine pylog_open

      subroutine pylog_close()
        integer(c_int) :: status

        status = nlsl_pylog_close()
        log_enabled = .false.
      end subroutine pylog_close

      subroutine pylog_emit(message)
        character(*), intent(in) :: message
        character(:), allocatable :: text
        character(kind=c_char), allocatable :: c_message(:)
        integer(c_int) :: status

        if (.not. log_enabled) return

        if (len_trim(message) == 0) then
           text = ''
        else
           text = trim(message)
        end if

        c_message = to_c_string(text)
        status = nlsl_pylog_debug(c_message)
        if (status /= 0) then
           log_enabled = .false.
        end if
      end subroutine pylog_emit

      subroutine flush_log_buffer()
        integer :: i
        if (.not. log_enabled) return
        if (.not. allocated(log_buffer)) return

        do i = 1, size(log_buffer)
           if (log_buffer(i)(1:1) == sentinel) exit
           call pylog_emit(log_buffer(i))
        end do
      end subroutine flush_log_buffer

      end module pylog_mod
