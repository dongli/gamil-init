program test_str_utils

    use utils

    character(30) :: string = "abc+de+fghi"
    character(1) :: substr = "+"
    integer c

    print *, count_substr(string, substr)
    print *, get_substr(string, substr, 1)
    print *, get_substr(string, substr, 2)
    print *, get_substr(string, substr, 3)

end program test_str_utils
