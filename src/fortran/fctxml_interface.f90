! This file is part of Cantera. See License.txt in the top-level directory or
! at https://cantera.org/license.txt for license and copyright information.

module fctxml
interface
    integer function fxml_new(name)
        character*(*), intent(in) :: name
    end function fxml_new

    integer function fxml_get_xml_file(file)
        character*(*), intent(in) :: file
    end function fxml_get_xml_file

    integer function fxml_clear()
    end function fxml_clear

    integer function fxml_del(i)
        integer, intent(in) :: i
    end function fxml_del

    integer function fxml_removechild(i, j)
        integer, intent(in) :: i
        integer, intent(in) :: j
    end function fxml_removechild

    integer function fxml_copy(i)
        integer, intent(in) :: i
    end function fxml_copy

    integer function fxml_preprocess_and_build(i, file)
        integer, intent(in) :: i
        character*(*), intent(in) :: file
    end function fxml_preprocess_and_build

    integer function fxml_attrib(i, key, value)
        integer, intent(in) :: i
        character*(*), intent(in) :: key
        character*(*), intent(out) :: value
    end function fxml_attrib

    integer function fxml_addattrib(i, key, value)
        integer, intent(in) :: i
        character*(*), intent(in) :: key
        character*(*), intent(in) :: value
    end function fxml_addattrib

    integer function fxml_addcomment(i, comment)
        integer, intent(in) :: i
        character*(*), intent(in) :: comment
    end function fxml_addcomment

    integer function fxml_tag(i, tag)
        integer, intent(in) :: i
        character*(*), intent(out) :: tag
    end function fxml_tag

    integer function fxml_value(i, value)
        integer, intent(in) :: i
        character*(*), intent(out) :: value
    end function fxml_value

    integer function fxml_child(i, loc)
        integer, intent(in) :: i
        character*(*), intent(in) :: loc
    end function fxml_child

    integer function fxml_child_bynumber(i, m)
        integer, intent(in) :: i
        integer, intent(in) :: m
    end function fxml_child_bynumber

    integer function fxml_findid(i, id)
        integer, intent(in) :: i
        character*(*), intent(in) :: id
    end function fxml_findid

    integer function fxml_findbyname(i, nm)
        integer, intent(in) :: i
        character*(*), intent(in) :: nm
    end function fxml_findbyname

    integer function fxml_nchildren(i)
        integer, intent(in) :: i
    end function fxml_nchildren

    integer function fxml_addchild(i, name, value)
        integer, intent(in) :: i
        character*(*), intent(in) :: name
        character*(*), intent(in) :: value
    end function fxml_addchild

    integer function fxml_addchildnode(i, j)
        integer, intent(in) :: i
        integer, intent(in) :: j
    end function fxml_addchildnode

    integer function fxml_write(i, file)
        integer, intent(in) :: i
        character*(*), intent(in) :: file
    end function fxml_write

    integer function ctml_getfloatarray(i, n, data, iconvert)
        integer, intent(in) :: i
        integer, intent(in) :: n
        integer, intent(in) :: iconvert
    end function ctml_getfloatarray

end interface
end module fctxml
