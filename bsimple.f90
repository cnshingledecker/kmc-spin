MODULE bsimple
  !
  !  Purpose:
  !    To define the derived data type used as a node in the
  !    binary tree, and to define the operations >, <. and ==
  !    for this data type.  This module also contains the
  !    subroutines to add a node to the tree, write out the
  !    values in the tree, and find a value in the tree.
  !
  !  Record of revisions:
  !      Date       Programmer          Description of change
  !      ====       ==========          =====================
  !    12/24/06    S. J. Chapman        Original code
  !
  USE typedefs
  USE parameters
  !  USE subroutines
  IMPLICIT NONE

  ! Restrict access to module contents.
  PRIVATE
  PUBLIC :: node, OPERATOR(>), OPERATOR(<), OPERATOR(==)
  PUBLIC :: add_node, write_node, find_node, find_min, delete_node, find_previous

  INTERFACE OPERATOR (>)
     MODULE PROCEDURE greater_than
  END INTERFACE OPERATOR (>)

  INTERFACE OPERATOR (<)
     MODULE PROCEDURE less_than
  END INTERFACE OPERATOR (<)

  INTERFACE OPERATOR (==)
     MODULE PROCEDURE equal_to
  END INTERFACE OPERATOR (==)

CONTAINS
  RECURSIVE SUBROUTINE add_node (ptr, new_node)
    !
    !  Purpose:
    !    To add a new node to the binary tree structure.
    !
    TYPE (node), POINTER :: ptr  ! Pointer to current pos. in tree
    TYPE (node), POINTER :: new_node ! Pointer to new node

    IF ( .NOT. ASSOCIATED(ptr) ) THEN
       ! There is no tree yet.  Add the node right here.
       ptr => new_node
    ELSE IF ( new_node < ptr ) THEN
       !      IF ( (new_node%iy) .EQ. 2 .AND. (new_node%iz .EQ. 2) ) PRINT *, new_node%tau,"<",ptr%tau
       IF ( ASSOCIATED(ptr%before) ) THEN
          CALL add_node ( ptr%before, new_node )
       ELSE
          ptr%before => new_node
          new_node%parent => ptr
          new_node%leftRight = 0
       END IF
    ELSE
       !      IF ( (new_node%iy) .EQ. 2 .AND. (new_node%iz .EQ. 2) ) PRINT *, new_node%tau,">",ptr%tau
       IF ( ASSOCIATED(ptr%after) ) THEN
          CALL add_node ( ptr%after, new_node )
       ELSE
          ptr%after => new_node
          new_node%parent => ptr
          new_node%leftRight = 1
       END IF
    END IF
  END SUBROUTINE add_node

  RECURSIVE SUBROUTINE write_node (ptr)
    !
    !  Purpose:
    !    To write out the contents of the binary tree
    !    structure in order.
    !
    TYPE (node), POINTER :: ptr  ! Pointer to current pos. in tree

    ! Write contents of previous node.
    IF ( ASSOCIATED(ptr%before) ) THEN
       CALL write_node ( ptr%before )
    END IF

    ! Write contents of current node.
    !WRITE (*,*) ptr%tau, ptr%ix, ptr%iy, ptr%iz

    ! Write contents of next node.
    IF ( ASSOCIATED(ptr%after) ) THEN
       CALL write_node ( ptr%after )
    END IF
  END SUBROUTINE write_node

  RECURSIVE SUBROUTINE find_node (ptr, search, error)
    !
    !  Purpose:
    !    To find a particular node in the binary tree structure.
    !    "Search" is a pointer to the name to find, and will
    !    also contain the results when the subroutine finishes
    !    if the node is found.
    !
    TYPE (node), POINTER :: ptr    ! Pointer to curr pos. in tree
    TYPE (node), POINTER :: search ! Pointer to value to find.
    INTEGER :: error               ! Error: 0 = ok, 1 = not found

    IF ( search < ptr ) THEN
       IF ( ASSOCIATED(ptr%before) ) THEN
          CALL find_node (ptr%before, search, error)
       ELSE
          error = 1
       END IF
    ELSE IF ( search == ptr ) THEN
       search = ptr
       error = 0
    ELSE
       IF ( ASSOCIATED(ptr%after) ) THEN
          CALL find_node (ptr%after, search, error)
       ELSE
          error = 1
       END IF
    END IF
  END SUBROUTINE find_node

  RECURSIVE SUBROUTINE find_previous (ptr, search, prevNode, leftRight, error)
    !
    !  Purpose:
    !    To find the node just prior to the one of interest
    !
    TYPE (node), POINTER :: ptr    ! Pointer to curr pos. in tree
    TYPE (node), POINTER :: search ! Pointer to value to find.
    TYPE (node), POINTER :: prevNode ! Pointer to value to find.
    INTEGER :: error               ! Error: 0 = ok, 1 = not found
    INTEGER :: leftRight
    INTEGER :: assocstat

    ! PRINT *, "Now in find previous"
    ! PRINT *, "ptr at",ptr%ix,ptr%iy,ptr%iz,ptr%tau
    ! PRINT *, "search at",search%ix,search%iy,search%iz,search%tau
    leftRight = 99
    assocstat = 0

    ! Find association status of the pointers
    IF ( ASSOCIATED(ptr%before) .AND. ASSOCIATED(ptr%after) ) THEN
       assocstat = 1
    ELSE IF ( ASSOCIATED(ptr%before) ) THEN !.AND. (.NOT. ASSOCIATED(ptr%after)))) THEN
       assocstat = 2
    ELSE IF ( ASSOCIATED(ptr%after) ) THEN !.AND. (.NOT. ASSOCIATED(ptr%before))) THEN
       assocstat = 3
    END IF

    SELECT CASE (assocstat)
    CASE(0)
       PRINT *, "Both pointers null! Exiting!!"
       PRINT *, "search=",search%tau,search%omega,search%ix,search%iy,search%iz
       PRINT *, "ptr=",ptr%tau,ptr%omega,ptr%ix,ptr%iy,ptr%iz
       PRINT *, "Now writing ptr"
       CALL write_node(ptr)
       PRINT *, "Now writing search:"
       CALL write_node(search)
       error = 1
       CALL EXIT()
    CASE(1)
       IF ( search == ptr%before ) THEN
          prevNode => ptr
          error = 0
          leftRight = 0
       ELSE IF ( search == ptr%after ) THEN
          !      PRINT *, search%tau,"=after",ptr%tau
          prevNode => ptr
          error = 0
          leftRight = 1
       ELSE IF ( search < ptr ) THEN
          !      PRINT *, search%tau," is less than ptr%before"
          CALL find_previous (ptr%before, search, prevNode, leftRight, error)
       ELSE
          CALL find_previous (ptr%after, search, prevNode, leftRight, error)
       END IF
    CASE(2)
       IF ( search == ptr%before ) THEN
          prevNode => ptr
          error = 0
          leftRight = 0
       ELSE
          CALL find_previous (ptr%before, search, prevNode, leftRight, error)
       END IF
    CASE(3)
       IF ( search == ptr%after ) THEN
          !      PRINT *, search%tau,"=after",ptr%tau
          prevNode => ptr
          error = 0
          leftRight = 1
       ELSE
          CALL find_previous (ptr%after, search, prevNode, leftRight, error)
       END IF
    END SELECT
    !    PRINT *, "End find_previous"
  END SUBROUTINE find_previous

  RECURSIVE SUBROUTINE find_min (ptr, search)
    IMPLICIT NONE
    TYPE (node), POINTER :: ptr    ! Pointer to curr pos. in tree
    TYPE (node), POINTER :: search ! Pointer to minimum value

    IF ( ASSOCIATED(ptr%before) ) THEN
       CALL find_min ( ptr%before, search )
    ELSE
       search => ptr
    END IF
  END SUBROUTINE find_min

  RECURSIVE SUBROUTINE find_max (ptr, search)
    IMPLICIT NONE
    TYPE (node), POINTER :: ptr    ! Pointer to curr pos. in tree
    TYPE (node), POINTER :: search ! Pointer to minimum value

    IF ( ASSOCIATED(ptr%after) ) THEN
       CALL find_min ( ptr%after, search )
    ELSE
       search => ptr
    END IF
  END SUBROUTINE find_max

  RECURSIVE SUBROUTINE delete_node (root, toDelete, prevNode, nextNode, error)
    !
    !  Purpose:
    !    To find a particular node in the binary tree structure.
    !    "Search" is a pointer to the name to find, and will
    !    also contain the results when the subroutine finishes
    !    if the node is found.
    !
    TYPE (node), POINTER :: toDelete ! Node to delete
    TYPE (node), POINTER :: root     ! Root node
    TYPE (node), POINTER :: prevNode, nextNode
    INTEGER :: error                 ! Error: 0 = ok, 1 = not found
    INTEGER :: dType                 ! Deletion type: 0 = leaf node, 1 = one child, 2 = two children
    INTEGER :: d1,d2,d3              ! Coordinates 1, 2, and 3
    INTEGER :: onetype
    INTEGER :: twotype

    d1 = toDelete%ix
    d2 = toDelete%iy
    d3 = toDelete%iz

    CALL node_type(toDelete,dType,onetype,twotype,0)
    IF ( VERBOSE .EQV. .TRUE. ) THEN
       IF ( dType .EQ. 1 ) PRINT *, "dType=",dType,"onetype=",onetype
       IF (dType .EQ. 2 ) PRINT *, "dType=",dType,"twotype=",twotype
    END IF

    SELECT CASE (dType)
    CASE (0)
       ! Since toDelete corresponds to a persistant lattice site, we just nullify
       ! the pointers and will re-add it to the tree later
       IF ( toDelete == root ) THEN
          NULLIFY(root)
          RETURN
       ELSE
          prevNode => toDelete%parent
          IF ( toDelete%leftRight .EQ. 0 ) THEN
             NULLIFY(prevNode%before,toDelete%parent)
          ELSE IF ( toDelete%leftRight .EQ. 1 ) THEN
             NULLIFY(prevNode%after,toDelete%parent)
          ELSE
             IF ( toDelete == root ) THEN
                CONTINUE
             ELSE
                PRINT *, "leftRight=",toDelete%leftRight
                !            CALL write_matrix_dataframe(root%ix,root%iy,root%iz)
                !            CALL write_matrix_dataframe(toDelete%ix,toDelete%iy,toDelete%iz)
                CALL EXIT()
             END IF
          END IF
       END IF
       toDelete%leftRight = -1
    CASE (1)
       SELECT CASE (onetype)
       CASE (0)
          PRINT *, "Error in onetype case: onetype=0"
          CALL EXIT()
       CASE (1)
          IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "onetype=",onetype
          nextNode => toDelete%before
          prevNode => toDelete%parent
          nextNode%parent => prevNode        ! 1
          nextNode%leftRight = 0             ! 2
          prevNode%before => nextNode        ! 3
       CASE (2)
          IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "onetype=",onetype
          nextNode => toDelete%before
          prevNode => toDelete%parent
          nextNode%parent => prevNode        ! 1
          prevNode%after =>  nextNode        ! 2
          nextNode%leftRight = 1             ! 3
       CASE (3)
          IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "onetype=",onetype
          nextNode => toDelete%after
          prevNode => toDelete%parent
          nextNode%parent => prevNode        ! 1
          prevNode%before => nextNode        ! 2
          nextNode%leftRight = 0             ! 3
       CASE (4)
          IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "onetype=",onetype
          nextNode => toDelete%after
          prevNode => toDelete%parent
          nextNode%parent => prevNode        ! 1
          prevNode%after => nextNode         ! 2
          nextNode%leftRight = 1             ! 3
       CASE (5)
          IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "onetype=",onetype
          nextNode => toDelete%before
          NULLIFY(nextNode%parent)           ! 1
          nextNode%leftRight = -1            ! 2
          root => nextNode                   ! 3
       CASE (6)
          IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "onetype=",onetype
          nextNode => toDelete%after
          NULLIFY(nextNode%parent)           ! 1
          nextNode%leftRight = -1            ! 2
          root => nextNode                   ! 3
       END SELECT
       IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "Nullifying pointers"
       NULLIFY(matrix(d1,d2,d3)%before,matrix(d1,d2,d3)%after,matrix(d1,d2,d3)%parent)
       matrix(d1,d2,d3)%leftRight = -1
       IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "Pointers nullified"
    CASE (2)
       ! First step is to delete the min of toDelete%after
       CALL find_min(toDelete%after,nextNode)
       toDelete => nextNode
       CALL delete_node(root,toDelete,prevNode,nextNode,error)
       SELECT CASE (twotype)
       CASE (1)
          IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "twotype=",twotype
          prevNode => matrix(d1,d2,d3)%before
          prevNode%parent => toDelete ! 1
          toDelete%before => prevNode ! 2
          root => toDelete
       CASE (2)
          IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "twotype=",twotype
          prevNode => matrix(d1,d2,d3)%before
          nextNode => matrix(d1,d2,d3)%after
          prevNode%parent => toDelete ! 1
          nextNode%parent => toDelete ! 2
          toDelete%before => prevNode ! 3
          toDelete%after  => nextNode ! 4
          root => toDelete
       CASE (3)
          IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "twotype=",twotype
          prevNode => matrix(d1,d2,d3)%before
          prevNode%parent => toDelete ! 1
          toDelete%before => prevNode ! 2
          root => toDelete
       CASE (4)
          IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "twotype=",twotype
          prevNode => matrix(d1,d2,d3)%before
          nextNode => matrix(d1,d2,d3)%after
          prevNode%parent => toDelete ! 1
          nextNode%parent => toDelete ! 2
          toDelete%before => prevNode ! 3
          toDelete%after  => nextNode ! 4
          root => toDelete
       CASE (5)
          IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "twotype=",twotype
          prevNode => matrix(d1,d2,d3)%before
          prevNode%parent => toDelete                     ! 1
          toDelete%before => prevNode                     ! 2
          toDelete%parent => matrix(d1,d2,d3)%parent      ! 3
          toDelete%leftRight = matrix(d1,d2,d3)%leftRight ! 4
          prevNode => matrix(d1,d2,d3)%parent
          IF ( toDelete%leftRight .EQ. 0 ) THEN           ! 5
             prevNode%before => toDelete
          ELSE
             prevNode%after  => toDelete
          END IF
       CASE (6)
          IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "twotype=",twotype
          prevNode => matrix(d1,d2,d3)%before
          nextNode => matrix(d1,d2,d3)%after
          prevNode%parent => toDelete                     ! 1
          nextNode%parent => toDelete                     ! 2
          toDelete%before => prevNode                     ! 3
          toDelete%after  => nextNode                     ! 4
          toDelete%parent => matrix(d1,d2,d3)%parent      ! 5
          toDelete%leftRight = matrix(d1,d2,d3)%leftRight ! 6
          prevNode => matrix(d1,d2,d3)%parent
          IF ( toDelete%leftRight .EQ. 0 ) THEN           ! 7
             prevNode%before => toDelete
          ELSE
             prevNode%after => toDelete
          END IF
       CASE (7)
          IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "twotype=",twotype
          prevNode => matrix(d1,d2,d3)%before
          prevNode%parent => toDelete                     ! 1
          toDelete%before => prevNode                     ! 2
          toDelete%parent => matrix(d1,d2,d3)%parent      ! 3
          toDelete%leftRight = matrix(d1,d2,d3)%leftRight ! 4
          prevNode => matrix(d1,d2,d3)%parent
          IF ( toDelete%leftRight .EQ. 0 ) THEN           ! 5
             prevNode%before => toDelete
          ELSE
             prevNode%after  => toDelete
          END IF
       CASE (8)
          IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "twotype=",twotype
          prevNode => matrix(d1,d2,d3)%before
          nextNode => matrix(d1,d2,d3)%after
          prevNode%parent => toDelete                     ! 1
          nextNode%parent => toDelete                     ! 2
          toDelete%before => prevNode                     ! 3
          toDelete%after  => nextNode                     ! 4
          toDelete%parent => matrix(d1,d2,d3)%parent      ! 5
          toDelete%leftRight = matrix(d1,d2,d3)%leftRight ! 6
          prevNode => matrix(d1,d2,d3)%parent
          IF ( toDelete%leftRight .EQ. 0 ) THEN           ! 7
             prevNode%before => toDelete
          ELSE
             prevNode%after => toDelete
          END IF
       END SELECT
       NULLIFY(matrix(d1,d2,d3)%before,matrix(d1,d2,d3)%after,matrix(d1,d2,d3)%parent)
       matrix(d1,d2,d3)%leftRight = -1
    END SELECT
  END SUBROUTINE delete_node

  LOGICAL FUNCTION greater_than (op1, op2)
    !
    !  Purpose:
    !    To test to see if operand 1 is > operand 2
    !    in alphabetical order.
    !
    TYPE (node), INTENT(IN) :: op1, op2

    IF (op1%tau > op2%tau) THEN
       greater_than = .TRUE.
    ELSE IF ( op1%tau == op2%tau ) THEN
       IF ( op1%ix > op2%ix ) THEN
          greater_than = .TRUE.
       ELSE IF ( op1%ix == op2%ix ) THEN
          IF ( op1%iy > op2%iy ) THEN
             greater_than = .TRUE.
          ELSE IF ( op1%iy == op2%iy ) THEN
             IF ( op1%iz > op2%iz ) THEN
                greater_than = .TRUE.
             ELSE
                greater_than = .FALSE.
             END IF
          ELSE
             greater_than = .FALSE.
          END IF
       ELSE
          greater_than = .FALSE.
       END IF
    ELSE
       greater_than = .FALSE.
    END IF
  END FUNCTION greater_than

  LOGICAL FUNCTION less_than (op1, op2)
    !
    !  Purpose:
    !    To test to see if operand 1 is < operand 2
    !    in alphabetical order.
    !
    TYPE (node), INTENT(IN) :: op1, op2


    IF (op1%tau < op2%tau) THEN
       less_than = .TRUE.
    ELSE IF ( op1%tau == op2%tau ) THEN
       IF ( op1%ix < op2%ix ) THEN
          less_than = .TRUE.
       ELSE IF ( op1%ix == op2%ix ) THEN
          IF ( op1%iy < op2%iy ) THEN
             less_than = .TRUE.
          ELSE IF ( op1%iy == op2%iy ) THEN
             IF ( op1%iz < op2%iz ) THEN
                less_than = .TRUE.
             ELSE
                less_than = .FALSE.
             END IF
          ELSE
             less_than = .FALSE.
          END IF
       ELSE
          less_than = .FALSE.
       END IF
    ELSE
       less_than = .FALSE.
    END IF
  END FUNCTION less_than

  LOGICAL FUNCTION equal_to (op1, op2)
    !
    !  Purpose:
    !    To test to see if operand 1 is equal to operand 2
    !    alphabetically.
    !
    TYPE (node), INTENT(IN) :: op1, op2

    IF ( (op1%tau == op2%tau ) .AND. &
         (op1%ix == op2%ix) .AND. &
         (op1%iy == op2%iy) .AND. &
         (op1%iz == op2%iz) ) THEN
       equal_to = .TRUE.
    ELSE
       equal_to = .FALSE.
    END IF
  END FUNCTION equal_to

  RECURSIVE SUBROUTINE node_type (temp, dType, onetype, twotype, recursive )
    !
    !  Purpose:
    !    To find a particular node in the binary tree structure.
    !    "Search" is a pointer to the name to find, and will
    !    also contain the results when the subroutine finishes
    !    if the node is found.
    !
    IMPLICIT NONE
    TYPE (node), POINTER :: temp ! Node to delete
    INTEGER :: dType, dType1, dType2                 ! Deletion type: 0 = leaf node, 1 = one child, 2 = two children
    INTEGER :: onetype, twotype
    LOGICAL :: leftBranch
    LOGICAL :: rightBranch
    !    LOGICAL :: recursive
    INTEGER :: recursive

    dType       = 0
    dType1      = 0
    dType2      = 0
    onetype     = 0
    twotype     = 0
    leftBranch  = .FALSE.
    rightBranch = .FALSE.

    IF ( ASSOCIATED(temp%before) ) THEN
       ! PRINT *, "Before node="
       CALL write_node(temp%before)
       dType = dType + 1
       leftBranch = .TRUE.
    END IF

    IF ( ASSOCIATED(temp%after)) THEN
       ! PRINT *, "After node="
       CALL write_node(temp%after)
       dType = dType + 1
       rightBranch = .TRUE.
    END IF

    ! Determine Type 1 case
    IF ( dType .EQ. 1 ) THEN
       IF ( (temp%leftRight .EQ. 0) .AND. (leftBranch .EQV. .TRUE. ) ) THEN
          onetype = 1
       ELSE IF ( (temp%leftRight .EQ. 1) .AND. (leftBranch .EQV. .TRUE.) ) THEN
          onetype = 2
       ELSE IF ( (temp%leftRight .EQ. 0) .AND. (rightBranch .EQV. .TRUE. ) ) THEN
          onetype = 3
       ELSE IF ( (temp%leftRight .EQ. 1) .AND. (rightBranch .EQV. .TRUE.) ) THEN
          onetype = 4
       ELSE IF ( (temp%leftRight .EQ. -1) .AND. (leftBranch .EQV. .TRUE.) ) THEN
          onetype = 5
       ELSE IF ( (temp%leftRight .EQ. -1) .AND. (rightBranch .EQV. .TRUE.) ) THEN
          onetype = 6
       END IF
    END IF

    IF ( dType .EQ. 2 ) THEN
       !      IF ( recursive .EQV. .TRUE. ) THEN
       IF ( recursive == 0 ) THEN
          !        recursive = .FALSE.
          CALL node_type(temp%before,dType1,onetype,twotype,1)
          CALL node_type(temp%after,dType2,onetype,twotype,1)
          IF ( .NOT. ASSOCIATED(temp%parent) ) THEN
             ! Root node
             IF ( (dType1 .EQ. 0) .AND. (dType2 .EQ. 0 ) ) THEN
                twotype = 1
             ELSE IF ( (dType1 .EQ. 0) .AND. (dType2 .NE. 0 ) ) THEN
                twotype = 2
             ELSE IF ( (dType1 .NE. 0) .AND. (dType2 .EQ. 0 ) ) THEN
                twotype = 3
             ELSE IF ( (dType1 .NE. 0) .AND. (dType2 .NE. 0 ) ) THEN
                twotype = 4
             END IF
          ELSE
             ! Not root node
             IF ( (dType1 .EQ. 0) .AND. (dType2 .EQ. 0 ) ) THEN
                twotype = 5
             ELSE IF ( (dType1 .EQ. 0) .AND. (dType2 .NE. 0 ) ) THEN
                twotype = 6
             ELSE IF ( (dType1 .NE. 0) .AND. (dType2 .EQ. 0 ) ) THEN
                twotype = 7
             ELSE IF ( (dType1 .NE. 0) .AND. (dType2 .NE. 0 ) ) THEN
                twotype = 8
             END IF
          END IF
       ELSE
          RETURN
       END IF
    END IF
  END SUBROUTINE node_type
END MODULE bsimple
