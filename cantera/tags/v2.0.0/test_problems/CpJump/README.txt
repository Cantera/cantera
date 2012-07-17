

#
#
# This test checks a new algorithm used in setState_HP(). Basically, we are trying
# to make the convergence algorithm fault tolerant of (small) jumps in the value of H or Cp
# at temperature boundaries. The new algorithm, which is a root finder, achieves this
# while the old algorithm, a bare Newton's method, diverges on this sample problem.
