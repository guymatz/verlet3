module kinds
  implicit none
  integer, PARAMETER, PUBLIC :: sp = SELECTED_REAL_KIND (p=6, r=37)
  integer, PARAMETER, PUBLIC :: dp = SELECTED_REAL_KIND (p=13, r=300)
  !integer, PARAMETER, PUBLIC :: dp = SELECTED_REAL_KIND (p=13, r=310)
end module kinds
