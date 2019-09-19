! This file is part of Cantera. See License.txt in the top-level directory or
! at https://cantera.org/license.txt for license and copyright information.

module cantera_xml

  use fctxml  ! interface for C functions

  type XML_Node
    integer :: xml_id
    integer :: wrapper
    integer :: err
  end type XML_Node

contains

    type(XML_Node) function new_XML_Node(src, name, wrap)
      implicit none
      character*(*), optional, intent(in) :: name
      character*(*), optional, intent(in) :: src
      integer, optional, intent(in) :: wrap
      type(XML_Node) self
      self%err = 0

      if (present(wrap)) then
         ! create a wrapper around an existing XML_Node
         self%xml_id = wrap
         self%wrapper = 1

      else if (present(src)) then
         !  read in an XML tree
         self%xml_id = fxml_get_xml_file(src)
         self%wrapper = 0

      else
         ! create an empty node
         self%xml_id = fxml_new(name)
         self%wrapper = 0
      end if
      new_XML_Node = self
      return
    end function new_XML_Node


    subroutine ctxml_clear()
      implicit none
      integer i
      i = fxml_clear()
    end subroutine ctxml_clear

    subroutine ctxml_getAttrib(self, key, value)
      implicit none
      type(XML_Node), intent(inout) :: self
      character*(*), intent(in) :: key
      character*(*), intent(out) :: value
      self%err = fxml_attrib(self%xml_id, key, value)
    end subroutine ctxml_getAttrib

    subroutine ctxml_addAttrib(self, key, value)
      implicit none
      type(XML_Node), intent(inout) :: self
      character*(*), intent(in) :: key
      character*(*), intent(in) :: value
      self%err = fxml_addattrib(self%xml_id, key, value)
    end subroutine ctxml_addAttrib

    subroutine ctxml_addComment(self, comment)
      implicit none
      type(XML_Node), intent(inout) :: self
      character*(*), intent(in) :: comment
      self%err = fxml_addcomment(self%xml_id, comment)
    end subroutine ctxml_addComment

    subroutine ctxml_getTag(self, tag)
      implicit none
      type(XML_Node), intent(inout) :: self
      character*(*), intent(out) :: tag
      self%err = fxml_tag(self%xml_id, tag)
    end subroutine ctxml_getTag

    subroutine ctxml_getValue(self, value)
      implicit none
      type(XML_Node), intent(inout) :: self
      character*(*), intent(out) :: value
      self%err = fxml_value(self%xml_id, value)
    end subroutine ctxml_getValue

    type(XML_Node) function ctxml_child(self, loc, id, name)
      implicit none
      type(XML_Node), intent(in) :: self
      character*(*), optional, intent(in) :: loc
      character*(*), optional, intent(in) :: id
      character*(*), optional, intent(in) :: name
      integer ichild

      if (present(loc)) then
         ichild = fxml_child(self%xml_id, loc)
      else if (present(id)) then
         ichild = fxml_findid(self%xml_id, id)
      else if (present(name)) then
         ichild = fxml_findbyname(self%xml_id, name)
      end if
      ctxml_child = new_XML_Node(wrap = ichild)
    end function ctxml_child

    integer function ctxml_nChildren(self)
      implicit none
      type(XML_Node), intent(in) :: self
      ctxml_nChildren = fxml_nchildren(self%xml_id)
    end function ctxml_nChildren

    subroutine ctxml_addChild(self, name, value)
      implicit none
      type(XML_Node), intent(inout) :: self
      character*(*), intent(in) :: name
      character*(*), intent(in) :: value
      integer ichild
      ichild = fxml_addchild(self%xml_id, name, value)
    end subroutine ctxml_addChild

    subroutine ctxml_write(self, file)
      implicit none
      type(XML_Node), intent(inout) :: self
      character*(*), intent(in) :: file
      self%err = fxml_write(self%xml_id, file)
    end subroutine ctxml_write

end module cantera_xml
